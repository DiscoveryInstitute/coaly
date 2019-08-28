#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include "Inputs.h"

using namespace std;

using LongType = int;

int main(int argc, char** argv) {

    // Input population size history
    const vector<double> popn_history = Inputs::createPopnHistoryFromString(argv[1]);
    const long ngens = popn_history.size();

    // Input mutation rate
    const double mutation_rate = stod(argv[2]);

    // Input founding diversity
    istringstream ss("0 " + ((argc-1)==3 ? string(argv[3]) : ""));
    vector<long> founding_diversity{istream_iterator<double>(ss),istream_iterator<double>()};

    // Run simulation
    mt19937_64 rng;
    vector<long> variant_frequencies;

    // Add founding diversity
    const long nic = founding_diversity.size() - 1;
    for (long ic=1; ic<=nic; ++ic) {
        for (long i=0; i<founding_diversity[ic]; ++i) variant_frequencies.push_back(ic);
    }

    const long firstgen = ngens-1;
    const long lastgen = 0L;

    for (long igen=firstgen; igen>=lastgen; --igen) {
        const double popn = popn_history[igen];
        // Add new variants by mutation.
        double expected_mutations = popn * mutation_rate;
        long actual_mutations = poisson_distribution<long>(expected_mutations)(rng);
        for(long i=0; i<actual_mutations; ++i) variant_frequencies.push_back(1L);
        // Special case: taking the final sample constitutes the first coalescent
        //               so skip the population fluctuations.
        if (igen==lastgen) break;
        // Get population fluctuations.
        const long new_popn = popn_history[igen-1];
        const size_t n = variant_frequencies.size();
        #pragma omp parallel for
        for(size_t i=0; i<n; ++i) {
            double new_freq = binomial_distribution<int>(new_popn, variant_frequencies[i]/popn)(rng);
            if (new_freq >= new_popn) new_freq = 0;
            variant_frequencies[i] = new_freq;
        }
        // Remove extinct variants.
        auto& vec = variant_frequencies; // use shorter name
        vec.erase(remove(vec.begin(), vec.end(), 0L), vec.end());
        // Log
        fprintf(stderr,"#gen %ld %ld %ld %ld\n", igen, (long)popn, actual_mutations, vec.size());
    }

    // Collate Results
    const int sample_popn = 1000;
    const long total_popn = popn_history[lastgen];
    vector<long> bins(sample_popn+1);
    for (long variant_frequency : variant_frequencies) {
        // This sampling constitutes the final fluctuation step.
        int i = binomial_distribution<int>(sample_popn, variant_frequency / (double)total_popn)(rng);
        bins[i]++;
    }

    // Print results
    double pre_i = 1.0/sample_popn;
    long pre_bin = sample_popn;
    for (int i=1; i<sample_popn; ++i) {
        printf("%lf %ld\n", pre_i*i, pre_bin*bins[i]);
    }

}

