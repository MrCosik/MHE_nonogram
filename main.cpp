#include <iostream>
#include <fstream>
#include<string>
#include <vector>
#include <sstream>
#include <random>

using namespace std;


bool equals(vector<int> a, vector<int> b) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < (int) a.size(); i++) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

struct Dimensions {

public:
    vector<vector<int>> daneX;
    vector<vector<int>> daneY;
    int dimension;


    vector<vector<int>> matrix;
    long dimensionX = daneX.size();
    long dimensionY = daneY.size();


    void showMatrix() {
        cout << "Resulting grid is: " << endl;
        for (int i = 0; i < dimensionX; i++) {
            for (int j = 0; j < dimensionY; j++) {
                if (j != 0) cout << " ";
                cout << ((matrix[i][j] == 1) ? 'X' : 'O');
            }
            cout << endl;
        }
    }



    bool updateLine(int idx) {
        int newVal, pos;
        bool hasChanged = false, go;

        vector<int> aux;
        vector<vector<int> > auxLines;

        for (int i = 1; i < (1 << dimensionX); i++) {
            go = true;
            aux.clear();
            for (int j = 0; j < dimensionX; j++) {
                if ((i & (1 << j)) != 0) newVal = 1;
                else newVal = 0;

                if (matrix[idx][j] != -1 and matrix[idx][j] != newVal) go = false;
                aux.push_back(newVal);
            }
            if (go) auxLines.push_back(aux);
        }

        for (int i = auxLines.size() - 1; i >= 0; i--) {
            aux.clear();
            newVal = pos = 0;
            while (pos < dimensionX) {
                if (auxLines[i][pos] == 0) {
                    if (newVal != 0) aux.push_back(newVal), newVal = 0;
                } else newVal++;
                pos++;
            }
            if (newVal != 0) aux.push_back(newVal);
            if (not equals(aux, daneX[idx])) auxLines.erase(auxLines.begin() + i);
        }

        if (auxLines.size() > 0) {
            for (int j = 0; j < dimensionX; j++) {
                if (matrix[idx][j] != -1) continue;
                go = true;
                newVal = auxLines[0][j];
                for (int i = 1; i < (int) auxLines.size(); i++) {
                    if (newVal != auxLines[i][j]) go = false;
                }
                if (go) matrix[idx][j] = newVal, hasChanged = true;
            }
        }

        return hasChanged;
    }

    bool updateColumn(int idx) {
        int newVal, pos;
        bool hasChanged = false, go;


        vector<int> aux;
        vector<vector<int> > auxLines;

        for (int i = 1; i < (1 << dimensionY); i++) {
            go = true;
            aux.clear();
            for (int j = 0; j < dimensionY; j++) {
                if ((i & (1 << j)) != 0) newVal = 1;
                else newVal = 0;

                if (matrix[j][idx] != -1 and matrix[j][idx] != newVal) go = false;
                aux.push_back(newVal);
            }
            if (go) auxLines.push_back(aux);
        }

        for (int i = auxLines.size() - 1; i >= 0; i--) {
            aux.clear();
            newVal = pos = 0;
            while (pos < dimensionY) {
                if (auxLines[i][pos] == 0) {
                    if (newVal != 0) aux.push_back(newVal), newVal = 0;
                } else newVal++;
                pos++;
            }
            if (newVal != 0) aux.push_back(newVal);
            if (not equals(aux, daneY[idx])) auxLines.erase(auxLines.begin() + i);
        }

        if (auxLines.size() > 0) {
            for (int j = 0; j < dimensionY; j++) {
                if (matrix[j][idx] != -1) continue;
                go = true;
                newVal = auxLines[0][j];
                for (int i = 1; i < (int) auxLines.size(); i++) {
                    if (newVal != auxLines[i][j]) go = false;
                }
                if (go) matrix[j][idx] = newVal, hasChanged = true;
            }
        }

        return hasChanged;
    }

    void solve() {
        matrix.resize(dimension);
        for (int i = 0; i < dimension; i++) matrix[i].resize(dimension, -1);

        bool finished = false;
        while (not finished) {
            finished = true;
            for (int i = 0; i < dimension; i++) {
                if (updateLine(i)) finished = false;
                if (updateColumn(i)) finished = false;
            }
        }
    }
};

int randomNumber(int to) {
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uni(0, to);

    int random_integer = uni(rng);

    return random_integer;
}

struct Nonogram {
public:
    vector<vector<int>> nonogram;
};


Nonogram generatePuzzle(int height, int width, int allElements) {


    int maxHeight;
    int maxWidth;

    vector<vector<int> > nonogram(
            height,
            vector<int>(width, 0));


    for (auto i = 0; i < allElements; i++) {
        maxHeight = randomNumber(height - 1);
        maxWidth = randomNumber(width - 1);
        if(nonogram[maxHeight][maxWidth] != 1) {
            nonogram[maxHeight][maxWidth] = 1;
        }else{
            i--;
        }
    }

    return {nonogram};
}

Dimensions createDimensionsFromNonogram(Nonogram nonogram, int height) {
    vector<vector<int>> daneX;
    vector<vector<int>> daneY;

    vector<int> xLine;
    vector<int> yLine;

    int x = 0;
    int y = 0;

    for (int i = 0; i < height ; i++) {
        for (int j = 0; j < nonogram.nonogram[0].size() ; ++j) {
            if (nonogram.nonogram[i][j] == 1) {
                x++;
            } else if(x != 0) {
                xLine.push_back(x);
                x = 0;
            }
        }
        if(x != 0){
            xLine.push_back(x);
        }
        daneX.push_back(xLine);
        xLine.clear();
        x = 0;
    }


    for (int j = 0; j < nonogram.nonogram[0].size(); ++j) {
        for (int i = 0; i < height; i++) {
            if (nonogram.nonogram[i][j] == 1) {
                y++;
            } else if(y != 0) {
                yLine.push_back(y);
                y = 0;
            }
        }
        if(y!= 0){
            yLine.push_back(y);
        }
        daneY.push_back(yLine);
        yLine.clear();
        y = 0;
    }

    return {daneX, daneY, 4};
}


ostream &operator<<(ostream &o, vector<vector<int>> tab) {
    for (auto row: tab) {
        o << "{ ";
        for (auto item: row) {
            o << item << " ";
        }
        o << " }";
    }
    return o;
};

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
    return {vectorX, vectorY, 4};
}


//
//
//vector<vector<string>> createRandomNonogram(int x, int y, int elements){
//
//
//
//}

int main() {
    vector<int> potentialNeighbour_x;
    vector<int> potentialNeighbour_y;

    Dimensions dimensions = getDataFromFile("../input.txt");
    cout << dimensions.daneX << endl;
    cout << dimensions.daneY << endl;
    dimensions.solve();
    dimensions.showMatrix();

    Nonogram testNonogram = generatePuzzle(4, 4,11);

    Dimensions testDimension = createDimensionsFromNonogram(testNonogram, 4);


    cout << testDimension.daneX << endl;
    cout << testDimension.daneY << endl;
    testDimension.solve();
    testDimension.showMatrix();

    return 0;
}








