#pragma once

#include <vector>

class Generation {

public:

    double total_population;
    double ancestral_population;
    std::vector<double> coalescence_probabilities;
    std::vector<double> number_of_descendants_distribution;
    std::vector<std::vector<double>> founding_allele_descendants_distribution;

    Generation(double total_popn, double sample_popn)
        : total_population(total_popn), ancestral_population(sample_popn),
          number_of_descendants_distribution(sample_popn)
        { number_of_descendants_distribution[1] = 1.0; }

    void coalesceOneGeneration(double total_popn);
    void coalesceLiteOneGeneration(double total_popn, int max_coalescent);
    void updateSample();

    void incrementAlleleFrequencySpectrumMutations(
            std::vector<double>& afs, const double mutation_rate) const;
    void incrementAlleleFrequencySpectrumFoundingDiversity(
            std::vector<double>& afs, const std::vector<double>& diversity) const;

    static std::vector<double> calculateAFS(const double sample_popn,
                                            const std::vector<double>& popn_history,
                                            const std::vector<double>& mutn_history,
                                            const std::vector<double>& founding_diversity);

};
