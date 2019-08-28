#pragma once

#include <vector>

class StocGeneration {

public:

    static std::vector<double> calculateAFS(const double sample_popn,
                                            const std::vector<double>& popn_history,
                                            const std::vector<double>& mutn_history,
                                            const std::vector<double>& founding_diversity,
                                            const double selection);

};

