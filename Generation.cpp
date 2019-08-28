#include "Generation.h"

#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

using namespace std;

const double EPSILON = numeric_limits<double>::epsilon();

/**
 *  Given lambda, the expected number of ancestral children per member of the population,
 *  calculate the probability of having 'ic' ancestral children on the assumption that
 *  the parent has one or more ancestral children. This is the probability of 'ic' lineages
 *  coalescing into one lineage. These probabilities are calculated up to 'max_coalescent'
 *  or until they are numerically negligible.
 *
 *  The 0th index of the array contains a special case: the probability of having no children.
 *  This quantity is useful for calculating the expected ancestral parent population.
 */


vector<double> getCoalescenceProbabilities(double n, double p, int max_coalescent) {
    const double lambda = n * p;
    const double eml = exp(-lambda);
    vector<double> coal_probs; coal_probs.reserve(max_coalescent+1);
    coal_probs.push_back(eml);         // special case: p(n=0)
    double coal_p = eml / (1.0 - eml); // p(n given n!=0)
    for (long ic=1; ic<=max_coalescent; ++ic) {
        coal_p *= lambda / ic; // probability that one ancestor has 'ic' ancestral children
        if (coal_p < EPSILON && lambda < ic) break;
        coal_probs.push_back(coal_p);
    }
    return coal_probs;
}
/*
vector<double> getCoalescenceProbabilities(double n, double p, int max_coalescent) {
    const double lambda = n * p;
    vector<double> coal_probs; coal_probs.reserve(max_coalescent+1);
    const double p0 = pow(1.0-p, n);
    coal_probs.push_back(p0);         // special case: p(n=0)
    double pre = 1.0 / (1.0 - p0); // p(n given n!=0)
    double logp = log(p);
    double logq = log(1.0-p);
    for (long ic=1; ic<=max_coalescent; ++ic) {
        double coal_p = pre*exp(lgamma(n+1) -lgamma(ic+1)-lgamma(n-ic+1) + ic*logp +(n-ic)*logq); // probability that one ancestor has 'ic' ancestral children
        if (coal_p < EPSILON && lambda < ic) break;
        coal_probs.push_back(coal_p);
    }
    return coal_probs;
}
*/

/**
 *  Calculate coalescence for one generation.
 */
void Generation::coalesceOneGeneration(double new_total_popn) {
    const auto& weights = this->number_of_descendants_distribution;
    const long sample_popn = weights.size()-1;
    const double n = this->ancestral_population;
    const double p = 1.0 / new_total_popn;

    vector<double> coal_probs = getCoalescenceProbabilities(n, p, sample_popn);

    this->total_population = new_total_popn;
    this->ancestral_population = new_total_popn * (1.0 - coal_probs[0]);

    this->coalescence_probabilities = move(coal_probs);
    updateSample();
}

/** Determines whether or not probabilities warrant an update of the sample */
inline bool maxed(vector<double>& probs, int imax) {
    while (!probs.empty() && probs.back() < EPSILON) probs.pop_back();
    return (probs.size() > (size_t)imax);
}

/**
 *  Calculate coalescence for one generation,
 *  but don't update whole sample unless the probabilities warrant it.
 */
void Generation::coalesceLiteOneGeneration(double new_total_popn, int max_coalescent){
    const auto& weights = this->number_of_descendants_distribution;
    const long sample_popn = weights.size();
    const double n = this->ancestral_population;
    const double p = 1.0 / new_total_popn;

    vector<double>& accum_probs = this->coalescence_probabilities;

    vector<double> coal_probs = getCoalescenceProbabilities(n, p, sample_popn);

    this->total_population = new_total_popn;
    this->ancestral_population = new_total_popn * (1.0 - coal_probs[0]);

    // If lots of coalescence this gen, update the sample to flush any existing probabilities.
    if (maxed(coal_probs, max_coalescent)) {
        updateSample(); // update to flush existing probabilities first
        accum_probs = move(coal_probs);
        updateSample(); // update with new probabilities
        return;
    }
    // If there are no accumulated probabilities (and current probabilities are not maxed),
    // the initial assignment is simple.
    if (accum_probs.empty()) {
        accum_probs = move(coal_probs);
        return;
    }

    // Compare this section with the section for coalescing the whole sample
    coal_probs.resize(max_coalescent+1);
    accum_probs.resize(max_coalescent+1);
    vector<double> new_accum_probs(max_coalescent+1);
    vector<double> uvec = accum_probs;
    vector<double> new_uvec(uvec.size());

    for (long i=1; i<=max_coalescent; ++i) new_accum_probs[i] = coal_probs[1] * uvec[i];
    for (long ic=2; ic<=max_coalescent; ++ic) {
        fill(new_uvec.begin(), new_uvec.end(), 0.0);
        for(long i=ic; i<=max_coalescent; ++i) {
            long pic = ic-1;
            for(long j=1; j<=i-pic; ++j) {
                new_uvec[i] += accum_probs[j] * uvec[i-j];
            }
            if (new_uvec[i] < EPSILON && i > 0 && new_uvec[i-1] >= new_uvec[i]) break;
        }
        uvec.swap(new_uvec);
        for(long i=ic;i<=max_coalescent;++i) {
            new_accum_probs[i] += coal_probs[ic] * uvec[i];
        }
    }

    this->coalescence_probabilities = move(new_accum_probs);
    if (maxed(accum_probs, max_coalescent)) updateSample();

}

/**
 *  Update the sample (number_of_descendants_distribution) using
 *  the stored coalescence probabilities,
 *  which might be just from this generation,
 *  or may have been accumulated over several generations.
 */
void Generation::updateSample() {

    const auto& weights = this->number_of_descendants_distribution;
    const long sample_popn = weights.size();

    vector<double>& coal_probs = this->coalescence_probabilities;
    while (!coal_probs.empty() && coal_probs.back()<EPSILON) coal_probs.pop_back();

    const long max_ic = coal_probs.size()-1;
    if (max_ic < 2) return;

    vector<double> new_weights(sample_popn, 0.0);
    // uvec_n is a helper vector that can be thought of as weights^n
    // uvec_1[i] is just weights[i]
    // uvec_2[i] is sum_j weights[j] * weights[i-j]
    // uvec_3[i] is sum_jk weights[j] * weights[k] * weights[i-j-k]
    vector<double> uvec = weights;
    vector<double> new_uvec(uvec.size());
    for(long i=1;i<sample_popn;++i) new_weights[i] += coal_probs[1] * uvec[i];

    // For each number of ancestral children ic >= 2
    for (long ic=2; ic<=max_ic; ++ic) {
        // Update uvec_ic -> uvec_(ic+1);
        fill(new_uvec.begin(), new_uvec.end(), 0.0);
        for(long i=ic; i<sample_popn; ++i) {
            long pic = ic-1;
            for(long j=1; j<=i-pic; ++j) {
                new_uvec[i] += weights[j] * uvec[i-j];
            }
            if (new_uvec[i] < EPSILON && i > 0 && new_uvec[i-1] >= new_uvec[i]) break;
        }
        uvec.swap(new_uvec);
        // Add a contribution to the weights
        for(long i=ic;i<sample_popn;++i) {
            new_weights[i] += coal_probs[ic] * uvec[i];
        }
    }

    this->number_of_descendants_distribution = move(new_weights);
    this->coalescence_probabilities.clear();

}

/**
 *  Add mutations at the given rate, incrementing the allele frequency spectrum using
 *  the number_of_descendants_distribution at this generation.
 */
void Generation::incrementAlleleFrequencySpectrumMutations(
        vector<double>& afs, const double mutation_rate) const {
    const double prefactor = this->ancestral_population * mutation_rate;
    const auto& weight = this->number_of_descendants_distribution;
    for (size_t i=1; i<weight.size(); ++i) afs[i] += prefactor * weight[i];
}

/** Combinatorial function, "N choose m" or C^n_m = n!/(m!(n-m)!) */
double n_choose_m(int ni, int mi) {
    double n = (double)ni;
    double m = (double)min(mi, ni-mi);
    double numerator = 1.0;
    for (int i=0; i<m; ++i) numerator *= (double)(n-i);
    double denominator = 1.0;
    for (int i=2; i<=m; ++i) denominator *= (double)i;
    return numerator / denominator;
}

/**
 *  Add genetic diversity from 1,2,3,... founders, incrementing the allele frequency spectrum
 *  using the number_of_descendants_distribution for that number of founders.
 */
void Generation::incrementAlleleFrequencySpectrumFoundingDiversity(
        vector<double>& afs, const vector<double>& founding_diversity) const {

    const vector<double>& weights = this->number_of_descendants_distribution;
    const long sample_popn = weights.size();
    const long max_ic = founding_diversity.size() - 1;
    if (max_ic < 1) return;

    // First transform founding diversity into ancestral founding diversity.
    vector<double> ancestral_diversity(max_ic+1);
    double pnonzero = this->ancestral_population / this->total_population;
    double pzero = 1.0 - pnonzero;
    for (int i=1;i<=max_ic;++i) {
        for (int j=i;j<=max_ic;++j) {
            double pij = pow(pnonzero,i) * pow(pzero,j-i) * n_choose_m(j,i);
            ancestral_diversity[i] += pij * founding_diversity[j];
        }
    }

    // Add contributions to the Allele Frequency Spectrum
    vector<double> uvec = weights;
    vector<double> new_uvec(uvec.size());
    //   For single founders
    for(long i=1;i<sample_popn;++i) afs[i] += ancestral_diversity[1] * uvec[i];
    //   For multiple founders
    for (int ic=2; ic<=max_ic; ++ic) {
        // Update uvec_ic -> uvec_(ic+1);
        fill(new_uvec.begin(), new_uvec.end(), 0.0);
        for(long i=ic; i<sample_popn; ++i) {
            long pic = ic-1;
            for(long j=1; j<=i-pic; ++j) {
                new_uvec[i] += weights[j] * uvec[i-j];
            }
            if (new_uvec[i] < EPSILON && i > 0 && new_uvec[i-1] >= new_uvec[i]) break;
        }
        uvec.swap(new_uvec);
        // Add a contribution to the weights
        for(long i=ic;i<sample_popn;++i) {
            afs[i] += ancestral_diversity[ic] * uvec[i];
        }
    }
}


/** Calculate AFS using the population history, mutation rate history, and any founding diversity. */
vector<double> Generation::calculateAFS(const double sample_popn,
                                        const vector<double>& popn_history,
                                        const vector<double>& mutn_history,
                                        const vector<double>& founding_diversity) {
    const int max_coalescent = 10;

    Generation generation (popn_history[0], sample_popn);

    vector<double> afs(sample_popn+1, 0.0);

    const long ngens = popn_history.size();
    for (long igen = 0; igen < ngens; ++igen) {
        generation.coalesceLiteOneGeneration(popn_history[igen], max_coalescent);
        generation.incrementAlleleFrequencySpectrumMutations(afs, mutn_history[igen]);
    }
    generation.incrementAlleleFrequencySpectrumFoundingDiversity(afs, founding_diversity);

    return afs;
}
