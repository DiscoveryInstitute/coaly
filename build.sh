rm -f Coaly Coaly.exe
g++ --std=c++11 -Wall -Wno-unknown-pragmas -g -O3 -o Coaly \
    Coaly.cpp FileReader.cpp CoalGeneration.cpp \
    StocGeneration.cpp ForeGeneration.cpp

