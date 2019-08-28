#include "ForeGeneration.h"

#include <cmath>
#include <limits>
#include <numeric>
#include <vector>
#include <cstdio>

using namespace std;

const double EPSILON = numeric_limits<double>::epsilon();

inline double logfact(double x) { return lgamma(x+1.0); }

void ForeGeneration::mutateOneGeneration(double mutation_rate) {
    vector<double>& afs = allele_frequency_spectrum;
    const long popn = afs.size()-1;
    afs[1] += popn * mutation_rate;
}

void ForeGeneration::generateOneGeneration(const double new_popn, const double selection) {

    const vector<double>& old_afs = allele_frequency_spectrum;
    const long old_n = old_afs.size()-1;

    vector<double> new_afs((long)new_popn+1);
    const long new_n = new_afs.size()-1;

    const double exps = exp(selection);

    new_afs[0] = old_afs[0];
    new_afs[new_n] = old_afs[old_n];
    for (long old_i=1; old_i<old_n; ++old_i) {
        double op = (double)old_i / (double)old_n;
        if (selection) {
            double opexps = op * exps;
            op = opexps / (1 + opexps - op);
        }
        long mid_j = (long)(new_popn*op);
        double mid_p = old_afs[old_i] * exp(logfact(new_n)-logfact(mid_j)-logfact(new_n-mid_j) + mid_j*log(op) + (new_n-mid_j)*log(1.0-op));
        double pj;
        pj = mid_p;
        for (long j=mid_j-1; j>=0; --j) {
            pj *= ((j+1) * (1.0-op)) / ((new_n-j) * op);
            if (pj < EPSILON) break;
            new_afs[j] += pj;
        }
        pj = mid_p;
        for (long j=mid_j; j<=new_n; ++j) {
            new_afs[j] += pj;
            pj *= ((new_n-j) * op) / ((j+1) * (1.0-op));
            if (pj < EPSILON) break;
        }
    }

    allele_frequency_spectrum = move(new_afs);
}

/** Calculate AFS using the population history, mutation rate history, and any founding diversity. */
vector<double> ForeGeneration::calculateAFS(const double sample_popn,
                                            const vector<double>& popn_history,
                                            const vector<double>& mutn_history,
                                            const vector<double>& founding_diversity,
                                            const double selection) {

    const long ngens = popn_history.size();

    ForeGeneration generation (popn_history[ngens-1], founding_diversity);

    long firstgen = ngens-1;
    long lastgen = 0;

    for (long igen = firstgen; igen > lastgen; --igen) {
        if (igen%100==0) fprintf(stderr,"# %ld\n",igen);
        double mutation_rate = mutn_history[igen];
        generation.mutateOneGeneration(mutation_rate);

        long nextgen = igen-1;
        double new_popn = nextgen>0 ? popn_history[nextgen] : sample_popn;
        generation.generateOneGeneration(new_popn, selection);
    }

    return generation.allele_frequency_spectrum;
}
