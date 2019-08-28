#include "StocGeneration.h"
#include <algorithm>
#include <random>

using namespace std;

/** Simulate AFS using the population history, mutation rate history, and any founding diversity. */
vector<double> StocGeneration::calculateAFS(const double sample_popn,
                                            const vector<double>& popn_history,
                                            const vector<double>& mutn_history,
                                            const vector<double>& founding_diversity,
                                            const double selection) {

    mt19937_64 rng;
    vector<long> variant_frequencies;

    // Add founding diversity
    const long nic = founding_diversity.size() - 1;
    for (long ic=1; ic<=nic; ++ic) {
        for (long i=0; i<founding_diversity[ic]; ++i) variant_frequencies.push_back(ic);
    }

    const long ngens = popn_history.size();
    const long firstgen = ngens-1;
    const long lastgen = 0L;
    const double exps = exp(selection);

    for (long igen=firstgen; igen>lastgen; --igen) {
        const double popn = popn_history[igen];
        // Add new variants by mutation.
        double expected_mutations = popn * mutn_history[igen];
        long actual_mutations = poisson_distribution<long>(expected_mutations)(rng);
        for(long i=0; i<actual_mutations; ++i) variant_frequencies.push_back(1L);
        // Get population fluctuations.
        const long nextgen = igen-1;
        const long new_popn = nextgen==lastgen ? sample_popn : popn_history[nextgen];
        const size_t n = variant_frequencies.size();
        #pragma omp parallel for
        for(size_t i=0; i<n; ++i) {
            double p = variant_frequencies[i]/popn;
            if (selection) {
                double pexps = p * exps;
                p = pexps / (1 + pexps - p);
            }
            double new_freq = binomial_distribution<int>(new_popn, p)(rng);
            if (new_freq == new_popn) new_freq = 0;
            variant_frequencies[i] = new_freq;
        }
        // Remove extinct variants.
        auto& vec = variant_frequencies; // use shorter name
        vec.erase(remove(vec.begin(), vec.end(), 0L), vec.end());
    }

    // Collate Results
    vector<double> bins(sample_popn+1);
    for (long variant_frequency : variant_frequencies) {
        bins[variant_frequency]++;
    }

    return bins;

}
