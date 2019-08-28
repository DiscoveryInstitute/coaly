#pragma once

#include <vector>

class ForeGeneration {

public:

    std::vector<double> allele_frequency_spectrum;

    ForeGeneration(double initial_popn,
                   const std::vector<double>& diversity = std::vector<double>{})
        : allele_frequency_spectrum(initial_popn+1)
        { std::copy(diversity.begin(), diversity.end(), allele_frequency_spectrum.begin()); }

    void generateOneGeneration(const double new_popn, const double selection);
    void mutateOneGeneration(const double mutation_rate);

    static std::vector<double> calculateAFS(const double sample_popn,
                                            const std::vector<double>& popn_history,
                                            const std::vector<double>& mutn_history,
                                            const std::vector<double>& founding_diversity,
                                            const double selection);

};
