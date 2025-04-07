#pragma once
#include <fstream>
#include <string>
#include <vector>

#include "Asset.h"

using namespace std;

class Portfolio{
private:
    vector<Asset> assets;

public:
    Portfolio(vector<Asset> assets):
        assets(assets){}


    void display(){
        cout<< "Name"<<"|\t"<<"Return"<<"|\t"<<"Volatility"<<"|\t" <<"Sharpe_ratio"<<"|\t"<<"Beta"<<"|\t";
        cout<<endl;

        for(int i=0;i<this->assets.size();i++)
        {
              this-> assets[i].display();
              cout<<" ";

            cout<<"\n";
        }
    }

    vector<Asset> get_assets(){
        return this->assets;
    }



};