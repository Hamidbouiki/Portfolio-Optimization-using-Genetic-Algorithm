#pragma once
#include <iostream>
#include <string>

using namespace std;


class Asset {
private:
    string name;
    double ret;
    double volatility;
    double sharpe_ratio;
    double beta;

public:
    Asset(string name, double ret, double volatility, double sharpe_ratio, double beta) :
            name(name), ret(ret), volatility(volatility), sharpe_ratio(sharpe_ratio), beta(beta){}

    void display(){
        cout << name << "|\t" << ret  << "|\t" << volatility  << "|\t" << sharpe_ratio  << "|\t" << beta <<"|\t" << endl;
    }
    string get_name(){
        return this->name;
    }


    double get_return(){
        return this->ret;
    }

    double get_sharp(){
        return  this->sharpe_ratio;
    }

    double get_beta(){
        return  this->beta;
    }
    double get_volatility(){
        return  this->volatility;
    }
};