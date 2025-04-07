#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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
        cout << name << "\t" << ret  << "\t" << volatility  << "\t" << sharpe_ratio  << "\t" << beta << endl;
    }
};

class Portifolio{
private:
    vector<Asset> assets;

public:
    void read_file(string file_name){
        ifstream myfile (file_name);

        string name;
        double ret;
        double volatility;
        double sharpe_ratio;
        double beta;

        if (myfile.is_open()){
            getline(myfile, name); // Discarding the first line 
            while(!myfile.eof()) {
                myfile >> name >> ret >> volatility >> sharpe_ratio >> beta;
                Asset newasset(name, ret, volatility, sharpe_ratio, beta);
                assets.push_back(newasset);
            }
        } else {
            cout << "Unable to open file";
        }
    }

    void display(){
        for(int i = 0; i < assets.size() ; i++){
            assets[i].display();
        }

        // // Alternative:
        // for(Asset asset : assets){
        //     asset.display();
        // }
    }
};

int main() {
    Portifolio portifolio;
    portifolio.read_file("file.txt");
    portifolio.display();
    return 0;
}