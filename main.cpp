#include <iostream>
#include <fstream>
#include<string>
#include <vector>
#include <sstream>
#include <random>
#include <algorithm>
#include <set>

using namespace std;


random_device rd;
mt19937 rng(rd());

bool equals(vector<int> a, vector<int> b) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < (int) a.size(); i++) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

struct Nonogram {
public:
    vector<vector<int>> nonogram;
};

struct Dimensions {

public:
    vector<vector<int>> daneX;
    vector<vector<int>> daneY;


    vector<vector<int>> matrix;
    long dimensionX = daneX.size();
    long dimensionY = daneY.size();

    bool operator==(const Dimensions &rhs) const {
        return daneX == rhs.daneX &&
               daneY == rhs.daneY;
    }

    bool operator!=(const Dimensions &rhs) const {
        return !(rhs == *this);
    }

    void showMatrix() {
        cout << "Resulting grid is: " << endl;
        for (int i = 0; i < dimensionX; i++) {
            for (int j = 0; j < dimensionY; j++) {
                if (j != 0) cout << " ";
                cout << ((matrix[i][j] == 1) ? 'X' : '.');
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

    Nonogram solve() {
        matrix.resize(dimensionY);
        for (int i = 0; i < dimensionY; i++) matrix[i].resize(dimensionX, -1);

        bool finished = false;
        while (not finished) {
            finished = true;
            for (int i = 0; i < dimensionX; i++) {
                if (updateLine(i)) finished = false;
            }
        }
        finished = false;

        while (not finished) {
            finished = true;
            for (int i = 0; i < dimensionY; i++) {
                if (updateColumn(i)) finished = false;
            }
        }
        return {matrix};
    }
};

int randomNumber(int to) {

    std::uniform_int_distribution<int> uni(0, to);

    int random_integer = uni(rng);

    return random_integer;
}


Nonogram generatePuzzle(int height, int width, int allElements) {


    int maxHeight;
    int maxWidth;

    vector<vector<int> > nonogram(
            height,
            vector<int>(width, 0));


    for (auto i = 0; i < allElements; i++) {
        maxHeight = randomNumber(height - 1);
        maxWidth = randomNumber(width - 1);
        if (nonogram[maxHeight][maxWidth] != 1) {
            nonogram[maxHeight][maxWidth] = 1;
        } else {
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

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < nonogram.nonogram[0].size(); ++j) {
            if (nonogram.nonogram[i][j] == 1) {
                x++;
            } else if (x != 0) {
                xLine.push_back(x);
                x = 0;
            }
        }
        if (x != 0) {
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
            } else if (y != 0) {
                yLine.push_back(y);
                y = 0;
            }
        }
        if (y != 0) {
            yLine.push_back(y);
        }
        daneY.push_back(yLine);
        yLine.clear();
        y = 0;
    }

    return {daneX, daneY};
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
    return {vectorX, vectorY};
}


double fitnessFunction(Dimensions original, Dimensions generated) {

    double fit = 0;

    int smallestSize;

    //check correctness of number of groups
    for (int i = 0; i < original.dimensionX; ++i) {
        fit += max(0, (int) (original.daneX[i].size() -
                             abs((int) original.daneX[i].size() - (int) generated.daneX[i].size())));
    }

    //check correctness of sums
    int generatedRowSum = 0;
    int originalRowSum = 0;

    for (int i = 0; i < generated.dimensionX; ++i) {
        for (int j: generated.daneX[i]) {
            generatedRowSum += j;
        }

        for (int j: original.daneX[i]) {
            originalRowSum += j;
        }

        fit += max(0, originalRowSum - abs(originalRowSum - generatedRowSum));

        generatedRowSum = 0;
        originalRowSum = 0;
    }

    //check correctness of each element
    for (int i = 0; i < original.dimensionX; ++i) {
        smallestSize = min(generated.daneX[i].size(), original.daneX[i].size());
        for (int j = 0; j < smallestSize; j++) {
            if (generated.daneX[i][j] == original.daneX[i][j]) {
                fit += 1;
            }
        }
    }

    //check correctness of number of groups
    for (int i = 0; i < original.dimensionY; ++i) {
        fit += max(0, (int) (original.daneY[i].size() -
                             abs((int) original.daneY[i].size() - (int) generated.daneY[i].size())));
    }

    //check correctness of sums
    int generatedColSum = 0;
    int originalColSum = 0;

    for (int i = 0; i < original.dimensionY; ++i) {
        for (int j: generated.daneY[i]) {
            generatedColSum += j;
        }

        for (int j: original.daneY[i]) {
            originalColSum += j;
        }

        fit += max(0, generatedColSum - abs(generatedColSum - originalColSum));

        generatedColSum = 0;
        originalColSum = 0;
    }


    //check correctness of each element
    for (int i = 0; i < original.dimensionY; ++i) {
        smallestSize = min(generated.daneY[i].size(), original.daneY[i].size());
        for (int j = 0; j < smallestSize; j++) {
            if (generated.daneY[i][j] == original.daneY[i][j]) {
                fit += 1;
            }
        }
    }

    return fit;

}


double goal(Dimensions original, Dimensions generated) {

    double ret = 0;
    int biggestSize;

    for (int i = 0; i < generated.dimensionX; i++) {
        biggestSize = max(generated.daneX[i].size(), original.daneX[i].size());
        if (generated.daneX[i].size() < biggestSize) {
            generated.daneX[i].resize(biggestSize, 0);
        }

        for (int j = 0; j < biggestSize; j++) {
            ret += abs(generated.daneX[i][j] - original.daneX[i][j]);
        }

    }

    for (int i = 0; i < generated.dimensionX; i++) {
        biggestSize = max(generated.daneY[i].size(), original.daneY[i].size());
        if (original.daneY[i].size() < biggestSize) {
            original.daneY[i].resize(biggestSize, 0);
        }

        for (int j = 0; j < biggestSize; j++) {
            ret += abs(generated.daneY[i][j] - original.daneY[i][j]);
        }
    }

    return ret;
}

vector<vector<int>> change1DvectorTo2D(vector<int> original, int rowSize) {
    vector<vector<int>> neighbourNonogram;

//    neighbourNonogram.resize(rowSize, vector<int>(columnSize));
    vector<int> neighbourNonogramRow;

    int min = 0;
    for (int j = 0; j < original.size(); j += rowSize) {
        for (int k = j; k < rowSize + j; k++) {
            neighbourNonogramRow.push_back(original[k]);
        }
        min += rowSize;
        neighbourNonogram.push_back(neighbourNonogramRow);
        neighbourNonogramRow.clear();
    }
    return neighbourNonogram;

}

vector<Dimensions> generateNeighbours(Dimensions original) {
    vector<Dimensions> result;
    Nonogram nonogram = original.solve();
    vector<int> bits;

    for (vector<int> el: nonogram.nonogram) {
        for (int j: el) {
            bits.push_back(j);
        }
    }

    result.push_back(original);

    for (int i = 0; i < bits.size(); ++i) {

        for (int j = 0; j < bits.size(); j++) {
            if (j + 1 < bits.size()) {
                bits[j] = bits[j + 1];
                bits[bits.size() - 1] = bits[0];
            }
        }
        nonogram = {change1DvectorTo2D(bits, original.dimensionX)};

        result.push_back(createDimensionsFromNonogram(nonogram, original.dimensionY));
    }
    return result;

}

Dimensions deterministicHillClimbing(Dimensions original) {
    vector<Dimensions> neighbours = generateNeighbours(original);


    Dimensions temporarySolution;
    int temporarySolutionFit;
    Dimensions bestSolution = neighbours[randomNumber(neighbours.size() - 1)];
    int bestSolutionFit;

    do {
        bestSolutionFit = fitnessFunction(original, bestSolution);
        for (auto &neighbour: neighbours) {
            if (fitnessFunction(original, bestSolution) < fitnessFunction(original, neighbour)) {
                temporarySolution = neighbour;
                temporarySolutionFit = fitnessFunction(original, temporarySolution);
            } else {
                temporarySolution = bestSolution;
                temporarySolutionFit = bestSolutionFit;
            }
        }


        if (temporarySolutionFit > bestSolutionFit) {
            bestSolution = temporarySolution;
        } else {
            break;
        }
        neighbours = generateNeighbours(bestSolution);

    } while (true);

    cout << "The best fitness score for deterministic hill climb: " << fitnessFunction(original, bestSolution) << endl;
    return bestSolution;

}

Dimensions randomHillClimbing(Dimensions original, int rounds) {
    vector<Dimensions> neighbours = generateNeighbours(original);


    Dimensions temporarySolution;
    Dimensions bestSolution = neighbours[randomNumber(neighbours.size() - 1)];

    int temporarySolutionFit;
    int bestSolutionFit;

    int randomIndex;
    for (int i = 0; i < rounds; ++i) {
        randomIndex = randomNumber(neighbours.size() - 1);
        bestSolutionFit = fitnessFunction(original, bestSolution);
        if (fitnessFunction(original, bestSolution) < fitnessFunction(original, neighbours[randomIndex])) {
            temporarySolution = neighbours[randomIndex];
            temporarySolutionFit = fitnessFunction(original, temporarySolution);
        } else {
            temporarySolution = bestSolution;
            temporarySolutionFit = bestSolutionFit;
        }

        if (temporarySolutionFit > bestSolutionFit) {
            bestSolution = temporarySolution;
            neighbours = generateNeighbours(bestSolution);
        }
    }
    cout << "The best fitness score for random hill climb: " << bestSolutionFit << endl;
    return bestSolution;

}

Dimensions tabuSearch(Dimensions original, int maxTabuSize, int iteration = -1) {
    vector<Dimensions> neighbours = generateNeighbours(original);

    int i = 0;
    Dimensions solution = neighbours[randomNumber(neighbours.size() - 1)];
    Dimensions bestSolution = solution;
    vector<Dimensions> tabuList = {solution};

    Dimensions temporarySolution;

    int temporarySolutionFit;
    int bestSolutionFit;
    int solutionFit;

    auto whileCondition = [](int counter, int rounds) -> bool {
        if (rounds < 0) {
            return true;
        } else {
            return counter < rounds;
        }
    };


    do {
        solutionFit = fitnessFunction(original, solution);
        bestSolutionFit = fitnessFunction(original, bestSolution);
        for (auto neighbour: neighbours) {
            if (fitnessFunction(original, bestSolution) < fitnessFunction(original, neighbour)) {
                temporarySolution = neighbour;
                temporarySolutionFit = fitnessFunction(original, temporarySolution);
            } else {
                temporarySolution = bestSolution;
                temporarySolutionFit = bestSolutionFit;
            }
        }


        if (find(tabuList.begin(), tabuList.end(), temporarySolution) != tabuList.end() &&
            bestSolutionFit < temporarySolutionFit) {
            bestSolution = temporarySolution;
        }

        if (solutionFit < bestSolutionFit) {
            solution = bestSolution;
        }

        tabuList.push_back(bestSolution);
        if (tabuList.size() > maxTabuSize) {
            tabuList.erase(tabuList.begin());
        }
        neighbours = generateNeighbours(bestSolution);
        i++;
    } while (whileCondition(i, iteration));

    cout << "The best fitness score for tabu search: " << fitnessFunction(original, solution) << endl;

    return solution;
}

Dimensions simulatedAnnealing(Dimensions original, int iterations) {
    vector<Dimensions> neighbours = generateNeighbours(original);
    int temp = 1;


    Dimensions temporarySolution;
    int temporarySolutionFit;
    Dimensions bestSolution = neighbours[randomNumber(neighbours.size() - 1)];
    int bestSolutionFit;

    for(int i = 0; i < iterations; i++) {
       bestSolutionFit = fitnessFunction(original, bestSolution);
       for (auto &neighbour: neighbours) {
           if (fitnessFunction(original, bestSolution) < fitnessFunction(original, neighbour)) {
               temporarySolution = neighbour;
               temporarySolutionFit = fitnessFunction(original, temporarySolution);
           } else {
               temporarySolution = bestSolution;
               temporarySolutionFit = bestSolutionFit;
           }
       }


       if (temporarySolutionFit > bestSolutionFit) {
           bestSolution = temporarySolution;
       } else {
           uniform_real_distribution<double> u(0,1);
           if (u(rng) < exp(-abs(fitnessFunction(original, temporarySolution)-fitnessFunction(original,temporarySolution))/ (1000 / temp))) {
               bestSolution = temporarySolution;
           }
       }

        temp++;
       neighbours = generateNeighbours(bestSolution);
   }

    cout << "The best fitness score for simulated annealing: " << fitnessFunction(original, bestSolution) << endl;

    return bestSolution;
}


int main() {
    vector<int> potentialNeighbour_x;
    vector<int> potentialNeighbour_y;

    Dimensions dimensions = getDataFromFile("../input.txt");
    cout << dimensions.daneX << endl;
    cout << dimensions.daneY << endl;
    dimensions.solve();
    dimensions.showMatrix();


    cout << "Best from random hillclimb: " << endl;
    Dimensions randomHillClmb = randomHillClimbing(dimensions, 10);
    cout << randomHillClmb.daneX << endl;
    cout << randomHillClmb.daneY << endl;
    randomHillClmb.solve();
    randomHillClmb.showMatrix();

    cout << "Best from deterministic hillclimb: " << endl;
    Dimensions deterministicHillClimb = deterministicHillClimbing(dimensions);
    cout << deterministicHillClimb.daneX << endl;
    cout << deterministicHillClimb.daneY << endl;
    deterministicHillClimb.solve();
    deterministicHillClimb.showMatrix();

    cout << "Best from tabu search: " << endl;
    Dimensions tabu = tabuSearch(dimensions, 200, 1000);
    cout << tabu.daneX << endl;
    cout << tabu.daneY << endl;
    tabu.solve();
    tabu.showMatrix();


    cout << "Best from simmulation annealing search: " << endl;
    Dimensions annealing = simulatedAnnealing(dimensions, 200);
    cout << annealing.daneX << endl;
    cout << annealing.daneY << endl;
    annealing.solve();
    annealing.showMatrix();

    return 0;
}








