#include <iostream>
#include <fstream>
#include<string>
#include <vector>
#include <sstream>

using namespace std;

class Dimensions {
public:
    vector<vector<int>> daneX;
    vector<vector<int>> daneY;

};

ostream &operator<<(ostream& o,vector<vector<int>> tab){
    for(auto row : tab){
        o << "{ ";
        for(auto item : row){
            o << item << " " ;
        }
        o << " }";
    }
    return o;
};

Dimensions getDataFromFile(string fileName);

int main() {
    vector<int> potentialNeighbour_x;
    vector<int> potentialNeighbour_y;

    Dimensions dimensions = getDataFromFile("../input.txt");


    cout << dimensions.daneX << endl;
    cout << dimensions.daneY << endl;


    return 0;
}




vector<int> simple_tokenizer(string s) {
    vector<int> vec;
    stringstream ss(s);
    string word;
    while (ss >> word) {
        vec.push_back(stoi(word));
    }
    return vec;
}

Dimensions getDataFromFile(string fileName) {
    vector<vector<int>> vectorX, vectorY;

    std::ifstream file(fileName);
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            //poprawic
            if (simple_tokenizer(line).empty()) {
                break;
            }
            vectorX.push_back(simple_tokenizer(line));
        }

        while (std::getline(file, line)) {
            vectorY.push_back(simple_tokenizer(line));
        }
        file.close();
    }
    return {vectorX,vectorY};
}


vector<vector<string>> createRandomNonogram(int x, int y, int elements){

}





