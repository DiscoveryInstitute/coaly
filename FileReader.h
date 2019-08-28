#pragma once

#include <string>
#include <vector>

class FileReader {
public:
    static std::vector<double> readPopulationHistory(const std::string& input_filename);
    static std::vector<double> readMutationRateHistory(const std::string& input_filename);
    static std::vector<double> readPrimordialFrequencies(const std::string& input_filename);
};
