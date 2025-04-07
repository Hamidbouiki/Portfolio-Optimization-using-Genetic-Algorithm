#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Asset.h"

using namespace std;

class DataProcessing  {
public:
    vector<Asset> assets;
    vector<vector<string>> content;
    string name;

public:
// This function will read the data from the csv File
    void read_file(string file_name){
        vector<string> row;
        string line, word;

        fstream file (file_name, ios::in);
        if(file.is_open())
        {
            while(getline(file, line))
            {
                row.clear();

                stringstream str(line);

                while(getline(str, word, ','))
                    row.push_back(word);
                content.push_back(row);

            }
        }
        else{
            cout<<"Could not open the file\n";
        }

    }

    //The CalculateIndicators() function will compute all the indicators for each asset, including Return, Volatility, Sharpe Ratio, and Beta.
    void CalculateIndicators(double risk_free) {
        double n_assets = content[0].size();
        vector<double> port_return(n_assets, 0.0);
        vector<double> port_risk(n_assets, 0.0);
        vector<double> port_beta(n_assets, 0.0);
        vector<double> mean_price(n_assets, 0.0);
        vector<double> sharpe_ratio(n_assets, 0.0);
        vector<vector<double>> port_cov(n_assets, vector<double>(n_assets, 0.0));
        //vector<double> assets_return;
        vector<vector<double>> assets_return;
        vector<double> row(n_assets, 0.0);
        for (int i = 0; i < content.size(); i++) {
            assets_return.push_back(row);
        }

        for (size_t j = 0; j < n_assets; ++j) {
            for (size_t i = 2; i < content.size(); ++i) {
                assets_return[i][j] = log(stod(content[i][j]) / stod(content[i - 1][j]));
                port_return[j] += log(stod(content[i][j]) / stod(content[i - 1][j]));
                mean_price[j] += stod(content[i][j])/(content.size() - 2) ;
            }
            port_return[j] = port_return[j] / (content.size() - 2) * 252 * 100;

        }

        for (size_t i = 2; i < content.size(); ++i) {
            for (size_t j = 0; j < n_assets; ++j) {
                //port_risk[j] += pow(assets_return[i][j] - port_return[j], 2);
                port_risk[j] += pow(stod(content[i][j]) - mean_price[j], 2);
                for (size_t k = 0; k < n_assets; ++k) {
                    port_cov[j][k] += (assets_return[i][j] - port_return[j]) * (assets_return[i][k] - port_return[k]);
                }
            }
        }
        for (size_t i = 0; i < n_assets; ++i) {
            port_risk[i] = sqrt(port_risk[i] / (content.size() - 2));
            for (size_t j = 0; j < n_assets; ++j) {
                port_cov[i][j] = port_cov[i][j] / (content.size() - 2);
            }

        }
        for (size_t i = 1; i < n_assets; ++i) {
        port_beta[i] = port_cov[0][i]/port_risk[0];
        sharpe_ratio[i] = (port_return[i] - risk_free * 100)/ port_risk[i];
        Asset newasset(content[0][i], port_return[i], port_risk[i], sharpe_ratio[i], port_beta[i]);
        assets.push_back(newasset);
        }
    }

    void display(){
        for(int i=0;i<this->content.size();i++)
        {
            for(int j=0;j<this->content[i].size();j++)
            {
                cout<<this->content[i][j]<<" ";
            }
            cout<<"\n";
        }
    }

    vector<Asset> get_assets(){
        return this->assets;
    }



};

