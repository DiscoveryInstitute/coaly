#include <cmath>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include "CoalGeneration.h"
#include "FileReader.h"
#include "ForeGeneration.h"
#include "StocGeneration.h"

using namespace std;

void error(const string& msg) { fprintf(stderr, "%s\n", msg.c_str()); exit(0); }


int main(int argc, char** argv) {
try {
    // Inputs
    const string method = string(argv[1]);
    const vector<double> popn_history = FileReader::readPopulationHistory(argv[2]);
    const vector<double> mutn_rate_history = FileReader::readMutationRateHistory(argv[3]);
    const vector<double> founding_diversity = (argc-1)>=4 ? FileReader::readPrimordialFrequencies(argv[4]) : vector<double>();
    const double selection = (argc-1)>=5 ? stod(argv[5]) : 0.0;

    // Run simulation
    const size_t sample_popn = min(1000.0, popn_history[0]);

    vector<double> afs;
    if (method=="stochastic")
        afs = StocGeneration::calculateAFS(sample_popn, popn_history, mutn_rate_history, founding_diversity, selection);
    if (method=="forward")
        afs = ForeGeneration::calculateAFS(sample_popn, popn_history, mutn_rate_history, founding_diversity, selection);
    if (method=="backward")
        afs = CoalGeneration::calculateAFS(sample_popn, popn_history, mutn_rate_history, founding_diversity, selection);

    printf("# %s\n", argv[1]);
    double pre_i = 1.0 / sample_popn;
    double pre_afs = sample_popn;
    for (size_t i=1; i<afs.size(); ++i) {
        printf("%lf %lf\n", pre_i * i, pre_afs * afs[i]);
    }


}
catch (const exception& e) { error(e.what()); }
}

