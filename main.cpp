#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include "Portfolio.h"
#include "GeneticAlgorithm.h"
#include "DataProcessing.h"

using namespace std;





int main() {
    DataProcessing Data;
    Data.read_file("Data.csv");
    Data.CalculateIndicators(0.01);
    Portfolio Portfolio(Data.get_assets());
    Portfolio.display();
    GeneticAlgorithm ga(&Portfolio, 500, 100, 0.8, 0.1, 0.2,0.01);
    ga.optimiseWeights();

    return 0;
}




