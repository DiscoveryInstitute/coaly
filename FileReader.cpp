#include "FileReader.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

inline long nearestLong(double x) { return (long)(x>0 ? x+0.5 : x-0.5); }


vector<vector<double>> readFile(string input_filename) {

    ifstream ifs(input_filename, ifstream::in);
    if (ifs.fail()) throw runtime_error("Cannot open " + input_filename);

    vector<vector<double>> rows;
    for (string line; getline(ifs,line); ) {
        if (line.length()==0) continue;
        istringstream tokens_ss(line);
        vector<double> row;
        for (string token; getline(tokens_ss, token, ' '); ) {
            if (!token.empty()) row.push_back(stod(token));
        }
        rows.push_back(move(row));
    }
    return rows;
}


vector<double> FileReader::readPopulationHistory(const string& input_filename) {
    const auto rows = readFile(input_filename);
    for(const auto& row : rows) if (row.size()!=2) {
        throw runtime_error("Need two values per line");
    }
    double prevGen = nearestLong(rows[0][0])+1;
    double prevPopn= rows[0][1];
    vector<double> history(prevGen);
    for(const auto& row : rows) {
        double nextGen = nearestLong(row[0]);
        double nextPopn = row[1];
        if (nextGen >= prevGen) throw runtime_error("Need decreasing generation numbers. Got" + to_string(prevGen) + " then " + to_string(nextGen));
        double norm = 1.0 / (nextGen - prevGen);
        for (long gen = prevGen-1; gen >= nextGen; --gen) {
            double t = (gen - prevGen) * norm;
            history[gen] = nearestLong((1-t)* prevPopn + t*nextPopn);
        }
        prevGen = nextGen;
        prevPopn = nextPopn;
    }
    return history;
}


vector<double> FileReader::readMutationRateHistory(const string& input_filename) {
    const auto rows = readFile(input_filename);
    for(const auto& row : rows) if (row.size()!=2) throw runtime_error("Need two values per line");

    double prevGen = nearestLong(rows[0][0])+1;
    double prevRate = rows[0][1];
    vector<double> history(prevGen);
    for(const auto& row : rows) {
        double nextGen = nearestLong(row[0]);
        double nextRate = row[1];
        if (nextGen >= prevGen) throw runtime_error("Need decreasing generation numbers. Got" + to_string(prevGen) + " then " + to_string(nextGen));
        double norm = 1.0 / (nextGen - prevGen);
        for (long gen = prevGen-1; gen >= nextGen; --gen) {
            double t = (gen - prevGen) * norm;
            history[gen] = (1-t)* prevRate + t*nextRate;
        }
        prevGen = nextGen;
        prevRate = nextRate;
    }
    return history;
}


vector<double> FileReader::readPrimordialFrequencies(const string& input_filename) {
    const auto rows = readFile(input_filename);

    vector<double> nus;
    for(const auto& row : rows) {
        size_t freq = nearestLong(row[0]);
        double density = row[1];
        if (freq+1 > nus.size()) nus.resize(freq+1, 0.0);
        nus[freq] = density;
    }
    return nus;
}



