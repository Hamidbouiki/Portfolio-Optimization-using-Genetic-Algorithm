#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include "Portfolio.h"

using namespace std;


class GeneticAlgorithm {
public:
    GeneticAlgorithm(
            Portfolio *data,
            int population,
            int generations,
            double crossover_rate,
            double mutation_rate,
            double elite_rate ,
            double riskFreeRate
    ) :
            data(*data),
            population(population),
            generations(generations),
            crossover_rate(crossover_rate),
            mutation_rate(mutation_rate),
            elite_rate(elite_rate),
            riskFreeRate(riskFreeRate)
            {
        n_assets = data->get_assets().size(); //data[0].size();
        port_return = vector<double>(n_assets, 0.0);
        port_risk = vector<double>(n_assets, 0.0);
        port_beta = vector<double>(n_assets, 0.0);
    }

    void optimiseWeights() {
        CalculatePortfolioIndicators();
        generateWeights();

        for (int i = 0; i < generations; ++i) {
            fitnessFunc();
            elitism();
            selection();
            crossover();
            mutation();
            avgGenResult();
            cout << "Generation " << i << ": Average fitness score of " <<avg_result << " from "
                 << weights.size() << " chromosomes" << endl;
        }

        fitnessFunc();
        optimalSolution();
    }

private:
    Portfolio data;
    int n_assets;
    int population;
    int generations;
    double crossover_rate;
    double mutation_rate;
    double elite_rate;
    double riskFreeRate ;
    vector<double> port_risk;
    vector<double> port_return;
    vector<double> excessReturn;
    vector<double> port_beta;
    vector<vector<double>> port_cov;
    vector<vector<double>> weights;
    vector<double> exp_ret;
    vector<double> portfolio_sharp;
    vector<double> port_volatilite;
    vector<double> treynorRatio;
    vector<double> fitness;
    vector<double> exp_beta;
    vector<double> sd;
    vector<int> non_elite_index;
    vector<int> crossover_index;
    vector<double> acc_fitness;
    double avg_result;
    vector<double> optimal_weight;
    double optimal_sharpe;

    void generateWeights() {
        if (population <= 0 || n_assets <= 0) {
            cout << "Invalid population or n_assets values." << endl;
        }
        weights.resize(population, vector<double>(n_assets, 0.0));

        for (int i = 0; i < population; ++i) {
            for (int j = 0; j < n_assets; ++j) {
                weights[i][j] = static_cast<double>(rand()) / RAND_MAX;
            }
            weights[i] = normalizeArray(weights[i]);
        }
    }

    void fitnessFunc() {
        exp_ret.resize(weights.size());
        exp_beta.resize(weights.size());
        treynorRatio.resize(weights.size());
        excessReturn.resize(weights.size());
        portfolio_sharp.resize(weights.size());
        fitness.resize(weights.size());
        sd.resize(weights.size());
        port_volatilite.resize(weights.size());


        for (size_t i = 0; i < weights.size(); ++i) {
            exp_ret[i] = 0.0;
            sd[i] = 0.0;
            portfolio_sharp[i] = 0.0;
            port_volatilite[i] = 0.0 ;
            exp_beta[i] = 0.0 ;
            treynorRatio[i] = 0.0;
            excessReturn[i]= 0.0  ;
            for (int j = 0; j < n_assets; ++j) {
                exp_ret[i] += weights[i][j] * port_return[j];
                exp_beta[i]+= weights[i][j] * port_beta[j];
                port_volatilite[i] += std::pow(weights[i][j] * data.get_assets()[j].get_volatility(), 2);
                portfolio_sharp[i] += weights[i][j] * data.get_assets()[j].get_sharp();

            }
            excessReturn[i] = exp_ret[i] - riskFreeRate;

                if (exp_ret[i]!= 0.0) {
                    treynorRatio[i] = exp_ret[i]/exp_beta[i];
                } else {
                    // Handle the case where portfolioBeta is 0 to avoid division by zero
                    treynorRatio[i] = 0.0;
                }

            port_volatilite[i] = std::sqrt(port_volatilite[i]);
            fitness[i] = exp_ret[i] - 0.5 * portfolio_sharp[i] - 0.25 * exp_beta[i]- 0.25 *port_volatilite[i]+0.25 *treynorRatio[i];


        }

    }

    void elitism() {
        /*
        Perform elitism step by finding n highest sharpe ratios.
        */

        int n_elite = static_cast<int>(fitness.size() * elite_rate);
        vector<int> temp(fitness.size());
        // Sort individuals in descending order of fitness during the selection process.
        iota(temp.begin(), temp.end(), 0);

        sort(temp.begin(), temp.end(), [this](int i, int j) { return fitness[i] > fitness[j]; });

        non_elite_index.resize(fitness.size() - n_elite);
        copy(temp.begin() + n_elite, temp.end(), non_elite_index.begin());
        // The highest fitness are directly transferred to the next generation without changes.
     std::vector<std::vector<double>> nextGeneration(temp.begin(), temp.begin() + n_elite);
    }

void selection() {
    // Compute Non-Elite Fitness
    vector<double> non_elite_fitness;
    for (int index : non_elite_index) {
        non_elite_fitness.push_back(fitness[index]);
    }

    // The number of pairs selected of individuals for crossover.
    int n_selections = non_elite_fitness.size() / 2;

    // Clear Previous Crossover Indices
    crossover_index.clear();

    // Normalize and calculate cumulative fitness
    acc_fitness = normalizeCumsum(non_elite_fitness);

    // Roulette wheel selection
    for (int _ = 0; _ < n_selections; ++_) {
        // Random probability between 0 and 1
        double rw_prob = static_cast<double>(rand()) / RAND_MAX;

        // Find the first index where acc_fitness is greater than or equal to rw_prob
        int index = 0;
        while (index < acc_fitness.size() && acc_fitness[index] < rw_prob) {
            ++index;
        }

        crossover_index.push_back(index);
    }
}
    // Crossover
    /*The purpose of the random crossover point is to introduce diversity in the offspring portfolios.
    By randomly selecting different crossover points during the evolutionary process,
    the algorithm explores different combinations of genetic material from the parents,
    aiming to discover potentially better solutions in the solution space.*/

    void crossover() {
    for (size_t i = 1; i < crossover_index.size(); i += 2) {
        int gen1 = crossover_index[i - 1];
        int gen2 = crossover_index[i];
        // Generates a random probability  between 0 and 1
        if (gen1 < weights.size() && gen2 < weights.size()) {
            for (int assetIndex = 0; assetIndex < n_assets; ++assetIndex) {
                double crossoverProb = static_cast<double>(rand()) / RAND_MAX;

                // Perform uniform crossover
                if (crossoverProb > crossover_rate) {
                    swap(weights[gen1][assetIndex], weights[gen2][assetIndex]);
                    //The use of swap simplifies the code for exchanging values between two variables
                }
            }

            // Normalize weights for both individuals
            weights[gen1] = normalizeArray(weights[gen1]);
            weights[gen2] = normalizeArray(weights[gen2]);
        }
    }
}

// Mutation
// Mutations are random modifications applied to individuals in the population
// to introduce genetic diversity and explore potential new solutions.
// Mutation helps prevent premature convergence towards suboptimal solutions.

void mutation() {
    // Calculate the total number of weights to be mutated
    //example : n_assets = 3, mutation_rate = 0.1, and crossover_index = {0, 1, 2, 3}
    int weight_n = crossover_index.size() * n_assets; // weight_n = 4 * 3 = 12

    // Calculate the number of weights to be mutated based on the mutation rate
    int mutate_gens = static_cast<int>(weight_n * mutation_rate); // mutate_gens = 12 * 0.1 = 1

    // Check if mutation_rate is non-zero to perform mutation
    if (mutation_rate != 0) {
        // Iterate over the number of weights to be mutated
        for (int _ = 0; _ < mutate_gens; ++_) {
            // Select a random index from crossover_index
            int rand_index = rand() % crossover_index.size(); // rand_index = random number between 0 and 3

            // Get the vector of weights at the selected index
            std::vector<double>& generation = weights[crossover_index[rand_index]];

            // Select a random asset index
            //rand_asset = 2 (randomly selected asset index)
            int rand_asset = rand() % n_assets; // rand_asset = random number between 0 and  2

            // Get the current value of the selected asset
            double mu_gen = generation[rand_asset];
            // Assume mutated_ind = 0.7 * mu_gen (random mutation factor)
            // Generate a random mutation factor between 0 and 1
            double mutated_ind = mu_gen * static_cast<double>(rand()) / RAND_MAX; // Update the asset with the mutated value

            // Update the asset with the mutated value
            generation[rand_asset] = abs(mutated_ind);

            // Normalize the weights to ensure the sum is 1.0
            generation = normalizeArray(generation);
        }
    }
}



    void CalculatePortfolioIndicators(){
        for (size_t i = 0; i < n_assets; ++i) {
            port_return[i] = data.get_assets()[i].get_return();
            port_risk[i]= data.get_assets()[i].get_volatility();
            port_beta[i] = data.get_assets()[i].get_beta();
        }
    }

    vector<double> normalizeArray(const vector<double> &arr) {
        vector<double> normal_arr(arr.size(), 0.0);
        double sum = 0.0;

        for (double val: arr) {
            sum += val;
        }

        for (size_t i = 0; i < arr.size(); ++i) {
            normal_arr[i] = arr[i] / sum;
        }

        return normal_arr;
    }

    vector<double> normalizeCumsum(const vector<double> &arr) {
        vector<double> normal_arr(arr.size(), 0.0);
        double sum = 0.0;

        for (double val: arr) {
            sum += val;
        }

        for (size_t i = 0; i < arr.size(); ++i) {
            normal_arr[i] = arr[i] / sum;
        }

        for (size_t i = 1; i < normal_arr.size(); ++i) {
            normal_arr[i] += normal_arr[i - 1];
        }

        return normal_arr;
    }

    void optimalSolution() {
        int optimal_index = distance(fitness.begin(), max_element(fitness.begin(), fitness.end()));
        optimal_weight = weights[optimal_index];
        optimal_sharpe = fitness[optimal_index];
        cout << "the optimal weights are :"<< endl;
        for(size_t i=0; i<optimal_weight.size(); i++){
            cout <<data.get_assets()[i].get_name()<<" : " <<optimal_weight[i]<<endl;
        }
    }


void avgGenResult() {
    avg_result = accumulate(portfolio_sharp.begin(), portfolio_sharp.end(), 0.0,
        [this](double acc, double sharpe) {
            size_t i = &sharpe - &portfolio_sharp[0]; // Get the index of the current element
            return acc + (exp_ret[i] - 0.5 * sharpe - 0.25 * exp_beta[i] - 0.25 * port_volatilite[i] + 0.25 * treynorRatio[i]);
        });

    avg_result = round(avg_result / portfolio_sharp.size() * 100) / 100.0;
}

};
