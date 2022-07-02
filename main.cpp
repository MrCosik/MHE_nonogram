#include <iostream>
#include <fstream>
#include<string>
#include <vector>
#include <sstream>
#include <random>
#include <algorithm>
#include <set>
#include <map>
#include <chrono>
#include <any>

using namespace std;

random_device rd;
mt19937 rng(rd());


auto arg = [](int argc, char **argv, std::string name, auto default_value) -> decltype(default_value) {
    using namespace std;
    string paramname = "";
    any ret = default_value;
    for (auto argument: vector<string>(argv, argv + argc)) {
        if ((argument.size() > 0) && (argument[0] == '-')) {
            if (paramname != "") {
                if (name == argument.substr(1))
                    ret = true;
            }
            paramname = argument.substr(1);
        } else if (name == paramname) {
            if (std::is_same_v<decltype(default_value), int>)
                ret = stoi(argument);
            else if (std::is_same_v<decltype(default_value), double>)
                ret = stod(argument);
            else if (std::is_same_v<decltype(default_value), char>)
                ret = argument.at(0);
            else if (std::is_same_v<decltype(default_value), bool>)
                ret = (argument == "true") || (argument == "1") || (argument == "yes");
            else
                ret = argument;
            paramname = "";
        }
    }
    return std::any_cast<decltype(default_value)>(ret);
};

void show_help() {
#ifndef __OPTIMIZE__
#define OPTIMIZED "optimized "
#else
#define OPTIMIZED ""
#endif
    std::cout << "# " << OPTIMIZED << "binary with date: " << __TIME__ << " " << __DATE__ << std::endl;
    std::cout << "-fname filename" << std::endl;
    std::cout << "-show_progress true/false" << std::endl;
}

bool equals(vector<int> a, vector<int> b) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < (int) a.size(); i++) {
        if (a[i] != b[i]) return false;
    }
    return true;
}


struct Nonogram1D {
public:
    vector<int> bits;
};

struct Nonogram {
public:
    vector<vector<int>> nonogram;

    Nonogram1D convertTo1D() {
        vector<int> result;

        for (vector<int> arr: nonogram) {
            for (int el: arr) {
                result.push_back(el);
            }
        }

        return {result};

    }

};

struct Dimensions {

public:
    vector<vector<int>> daneX;
    vector<vector<int>> daneY;


    vector<vector<int>> matrix;
    long dimensionX = daneX.size();
    long dimensionY = daneY.size();

    bool operator==(Dimensions rhs) const {

        for (int i = 0; i < daneX.size(); ++i) {
            for (int j = 0; j < daneX[i].size(); ++j) {
                if (daneX[i][j] != rhs.daneX[i][j]) {
                    return false;
                }
            }
        }

        for (int i = 0; i < daneY.size(); ++i) {
            for (int j = 0; j < daneY[i].size(); ++j) {
                if (daneY[i][j] != rhs.daneY[i][j]) {
                    return false;
                }
            }
        }

        return true;
    }

    bool operator!=(const Dimensions &rhs) const {
        return !(rhs == *this);
    }

    int getNumberOfElements() {
        int elements = 0;
        for (vector<int> arr: daneX) {
            for (int el: arr) {
                elements += el;
            }
        }
        return elements;
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

Dimensions createDimensionsFromNonogram1D(Nonogram1D nonogram, int height) {
    vector<vector<int>> daneX;
    vector<vector<int>> daneY;

    vector<int> xLine;
    vector<int> yLine;

    int x = 0;
    int y = 0;

    for (int i = 0; i < height; i++) {
        if (nonogram.bits[i] == 1) {
            x++;
        }
        if (x != 0) {
            xLine.push_back(x);
        }
        daneX.push_back(xLine);
        xLine.clear();
        x = 0;
    }


    for (int j = 0; j < nonogram.bits.size(); ++j) {
        for (int i = 0; i < height; i++) {
            if (nonogram.bits[i] == 1) {
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
//    Nonogram nonogram = generatePuzzle(original.dimensionY, original.dimensionX, original.getNumberOfElements());
    vector<int> bits;

    for (vector<int> el: nonogram.nonogram) {
        for (int j: el) {
            bits.push_back(j);
        }
    }


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

vector<Dimensions> generateRandomNeighbours(Dimensions original) {
    vector<Dimensions> result;
//    Nonogram nonogram = original.solve();
    Nonogram nonogram = generatePuzzle(original.dimensionY, original.dimensionX, original.getNumberOfElements());
    vector<int> bits;

    for (vector<int> el: nonogram.nonogram) {
        for (int j: el) {
            bits.push_back(j);
        }
    }


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
    vector<Dimensions> neighbours = generateRandomNeighbours(original);


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

    cout << "The best fitness score for deterministic hill climb: " << bestSolutionFit << endl;
    return bestSolution;

}

Dimensions randomHillClimbing(Dimensions original, int rounds) {
    vector<Dimensions> neighbours = generateRandomNeighbours(original);


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
//                    cout << "Fit function: " << bestSolutionFit << endl;
    }
    cout << "The best fitness score for random hill climb: " << bestSolutionFit << endl;
    return bestSolution;

}

Dimensions tabuSearch(Dimensions original, int maxTabuSize, int iteration = -1) {
    vector<Dimensions> neighbours = generateRandomNeighbours(original);

    int i = 0;
    Dimensions bestSolution = neighbours[randomNumber(neighbours.size() - 1)];
    Dimensions temporarySolution;
    vector<Dimensions> tabuList = {bestSolution};
    vector<Dimensions> tabuSteps = {bestSolution};

    int temporarySolutionFit;
    int bestSolutionFit;

    int bestPossibleFitness = fitnessFunction(original, original);

    bool changed;

    auto whileCondition = [](int counter, int rounds) -> bool {
        if (rounds < 0) {
            return true;
        } else {
            return counter < rounds;
        }
    };

    auto termination = [=](Dimensions tempSol) {
        for (Dimensions el: tabuList) {
            if (el == tempSol) {
                return false;
            }
        }
        return true;
    };

    temporarySolution = bestSolution;

    do {
        changed = false;
        bestSolutionFit = fitnessFunction(original, bestSolution);

        for (auto neighbour: neighbours) {
            if (fitnessFunction(original, temporarySolution) < fitnessFunction(original, neighbour)) {
                temporarySolution = neighbour;
                temporarySolutionFit = fitnessFunction(original, neighbour);
            }
        }


        if (termination(temporarySolution) && temporarySolutionFit > bestSolutionFit) {
            bestSolution = temporarySolution;
            bestSolutionFit = temporarySolutionFit;
            tabuList.push_back(bestSolution);
            tabuSteps.push_back(bestSolution);
            changed = true;
        }

        if (tabuList.size() > maxTabuSize) {
            tabuList.erase(tabuList.begin());
        }

        if (bestSolutionFit == bestPossibleFitness) {
            break;
        }

        if (!changed && tabuSteps.size() > 0) {
            bestSolution = tabuSteps[tabuSteps.size() - 1];
            tabuSteps.pop_back();
        }


        neighbours = generateNeighbours(bestSolution);
        i++;
    } while (whileCondition(i, iteration));

    cout << "The best fitness score for tabu search: " << bestSolutionFit << endl;

    return bestSolution;
}

Dimensions simulatedAnnealing(Dimensions original, int iterations, double temp) {
    vector<Dimensions> neighbours = generateRandomNeighbours(original);

    Dimensions temporarySolution;
    int temporarySolutionFit;
    Dimensions bestSolution = neighbours[randomNumber(neighbours.size() - 1)];
    int bestSolutionFit;
    Dimensions finalSolution = bestSolution;
    int finalSolutionFit = fitnessFunction(original, finalSolution);


    int randomIndex;
    for (int i = 0; i < iterations; i++) {
        neighbours = generateNeighbours(bestSolution);
        randomIndex = randomNumber(neighbours.size() - 1);

        bestSolutionFit = fitnessFunction(original, bestSolution);
        if (fitnessFunction(original, bestSolution) < fitnessFunction(original, neighbours[randomIndex])) {
            temporarySolution = neighbours[randomIndex];
            temporarySolutionFit = fitnessFunction(original, temporarySolution);
        } else {
            uniform_real_distribution<double> distr(0.0, 1.0);
//
//                cout << exp(-abs(
//                        fitnessFunction(original, bestSolution) - fitnessFunction(original, neighbours[randomIndex])) /
//                            (1000 / temp)) << endl;

            if (distr(rng) <
                exp(-abs(fitnessFunction(original, bestSolution) - fitnessFunction(original, neighbours[randomIndex])) /
                    (1000 / temp))) {
                bestSolution = neighbours[randomIndex];
                bestSolutionFit = fitnessFunction(original, neighbours[randomIndex]);
            }
        }

        temp++;

    }

    if (finalSolutionFit < bestSolutionFit) finalSolution = bestSolution;

    cout << "The best fitness score for simulated annealing: " << fitnessFunction(original, bestSolution) << endl;

    return bestSolution;

}


int selectPop(vector<double> arr) {
    std::uniform_int_distribution<int> distr(0, arr.size() - 1);
    int idx1 = distr(rng);
    int idx2 = distr(rng);
    if (arr[idx1] > arr[idx2]) return idx1;
    return idx2;

}

pair<Nonogram1D, Nonogram1D> crossoverFunction(Nonogram1D parent1, Nonogram1D parent2, double probability) {
    uniform_real_distribution<double> distr(0.0, 1.0);
    if (distr(rng) < probability) {
        uniform_int_distribution<int> distr(0, parent1.bits.size() - 1);
        int crossPtk = distr(rng);
        for (int i = 0; i < crossPtk; i++) {
            swap(parent1.bits[i], parent2.bits[i]);
        }
    }

    return {parent1, parent2};
}

pair<Nonogram1D, Nonogram1D> doubleCross(Nonogram1D parent1, Nonogram1D parent2, double probability) {
    uniform_real_distribution<double> distr(0.0, 1.0);
    if (distr(rng) < probability) {
        uniform_int_distribution<int> distr(1, parent1.bits.size() - 2);
        int crossPtk1 = distr(rng);
        int crossPtk2 = distr(rng);
        if (crossPtk1 > crossPtk2) swap(crossPtk1, crossPtk2);

        for (int i = crossPtk1; i <= crossPtk2; i++) {
            swap(parent1.bits[i], parent2.bits[i]);
        }

    }

    return {parent1, parent2};
}


Nonogram1D mutation(Nonogram1D element, double probability) {
    uniform_real_distribution<double> distr(0.0, 1.0);
    if (distr(rng) < probability) {
        uniform_real_distribution<float> distr(0, element.bits.size() - 1);
        int a = distr(rng);
        int b = distr(rng);
        swap(element.bits[a], element.bits[b]);

    }
    return element;
}

Nonogram1D inverseMutation(Nonogram1D element, double probability) {
    uniform_real_distribution<double> distr(0.0, 1.0);
    if (distr(rng) < probability) {
        uniform_real_distribution<float> distr(0, element.bits.size() - 1);
        int rangeStart = distr(rng);
        int rangeEnd = distr(rng);
        if (rangeStart > rangeEnd) {
            swap(rangeStart, rangeEnd);
        }

        int j = rangeEnd;
        for (int i = 0; i <= rangeEnd - rangeStart; i++) {

            iter_swap(element.bits.begin() + rangeStart, element.bits.begin() + rangeEnd);
            rangeStart++;
            rangeEnd--;
        }

    }
    return element;
}

vector<Dimensions> generateInitialPopulation(int amount, int height, int width, int memberSize) {
    vector<Dimensions> result;
    Nonogram populationMember;
    for (int i = 0; i < amount; i++) {
        populationMember = generatePuzzle(height, width, memberSize);
        result.push_back(createDimensionsFromNonogram(populationMember, populationMember.nonogram[0].size()));
    }
    return result;
}

vector<Nonogram1D> convertPopulationTo1D(vector<Dimensions> dimensions) {
    vector<Nonogram1D> result;
    Nonogram temporaryNonogram;

    for (Dimensions el: dimensions) {
        temporaryNonogram = el.solve();
        result.push_back(temporaryNonogram.convertTo1D());
    }

    return result;
}

int getNumberOfElements(Dimensions element) {
    int sum = 0;

    for (vector<int> el: element.daneX) {
        sum += accumulate(el.begin(), el.end(), 0);
    }

    return sum;
}

Dimensions geneticAlgorithm(Dimensions original, int iterations, int populationSize, double crossingP, double mutationP,
                            double fitnessIterations) {
    int i = 0;
    int fitIterations = 0;
    auto whileCondition = [&](int counter) -> bool {
        i++;
        return i < counter;
    };

    auto calculateFitness = [&](auto pop) {
        std::vector<double> ret;
        for (auto e: pop) ret.push_back(fitnessFunction(original, e));
        return ret;
    };

    auto population = generateInitialPopulation(populationSize, original.dimensionX, original.dimensionY,
                                                getNumberOfElements(original));
    vector<double> fitForPop = calculateFitness(population);

    auto populationInBits = convertPopulationTo1D(population);

    while (whileCondition(iterations)) {
        vector<Nonogram1D> parents;
        vector<Nonogram1D> children;

        for (int i = 0; i < population.size(); i++) {
            parents.push_back(populationInBits[selectPop(fitForPop)]);
        }
        for (int i = 0; i < parents.size(); i += 2) {
            auto[a, b] = crossoverFunction(parents[i], parents[i + 1], crossingP);
//            auto[a, b] = doubleCross(parents[i], parents[i + 1], 0.9);
            children.push_back(a);
            children.push_back(b);
        }
        for (int i = 0; i < parents.size(); i++) {
//            children[i] = mutation(children[i], 0.1);
            children[i] = inverseMutation(children[i], mutationP);
        }

        fitForPop.clear();
        for (Dimensions el: population) {
            fitForPop.push_back(fitnessFunction(original, el));
            fitIterations++;
        }
        if (fitIterations == fitnessIterations) {
            break;
        }

    };

    cout << "Fitness function for genetic: " << fitnessFunction(original, population[
            max_element(fitForPop.begin(), fitForPop.end())
            - fitForPop.begin()]) << endl;

    return population[
            max_element(fitForPop.begin(), fitForPop.end())
            - fitForPop.begin()];
}

int main(int argc, char **argv) {
    Dimensions best_solution;

    Dimensions dimensions = getDataFromFile("../input.txt");
    cout << dimensions.daneX << endl;
    cout << dimensions.daneY << endl;
    dimensions.solve();
    dimensions.showMatrix();


    auto fname = arg(argc, argv, "fname", std::string(""));
    auto iterations = arg(argc, argv, "iterations", 1000);
    auto method = arg(argc, argv, "method", std::string("random_hillclimb"));
    auto tabuSize = arg(argc, argv, "tabuSize", 200);
    auto temperature = arg(argc, argv, "temperature", 20);

    auto pop_size = arg(argc, argv, "pop_size", 1000);
    auto genIteration = arg(argc, argv, "genIteration", -1);
    auto crossover_p = arg(argc, argv, "crossover_p", 0.9);
    auto mutation_p = arg(argc, argv, "mutation_p", 0.1);
    auto help = arg(argc, argv, "help", true);

    /// do we need to show some help?
    if (arg(argc, argv, "h", false)) {
        show_help();
        return 0;
    }

    /// load the dimensions at hand
    std::chrono::duration<double> calculation_duration;
//    if (method == "deterministic_hill") {
    cout << "Best from deterministic hillclimb: " << endl;
    Dimensions deterministicHillClimb = deterministicHillClimbing(dimensions);
    cout << deterministicHillClimb.daneX << endl;
    cout << deterministicHillClimb.daneY << endl;
    deterministicHillClimb.solve();
    deterministicHillClimb.showMatrix();
//    } else if (method == "random_hillclimb") {
    cout << "Best from random hillclimb: " << endl;
    Dimensions randomHillClmb = randomHillClimbing(dimensions, iterations);
    cout << randomHillClmb.daneX << endl;
    cout << randomHillClmb.daneY << endl;
    randomHillClmb.solve();
    randomHillClmb.showMatrix();

//    } else if (method == "tabu") {
    cout << "Best from tabu search: " << endl;
    Dimensions tabu = tabuSearch(dimensions, tabuSize, iterations);
    cout << tabu.daneX << endl;
    cout << tabu.daneY << endl;
    tabu.solve();
    tabu.showMatrix();

//    } else if (method == "simulated_annealing") {
    cout << "Best from simmulation annealing search: " << endl;
    Dimensions annealing = simulatedAnnealing(dimensions, iterations, 100);
    cout << annealing.daneX << endl;
    cout << annealing.daneY << endl;
    annealing.solve();
    annealing.showMatrix();
//    } else if (method == "genetic_algorithm") {
    cout << "Best from genetic: " << endl;
    Dimensions genetic = geneticAlgorithm(dimensions, iterations, pop_size, crossover_p, mutation_p, genIteration);
    cout << genetic.daneX << endl;
    cout << genetic.daneY << endl;
    genetic.solve();
    genetic.showMatrix();
//    } else {
//        std::cerr << "unknown method" << std::endl;
//    }
//    cout << "# " << method << ": best_cost: " << fitnessFunction(dimensions, best_solution)
//         << " calculation_time: " << calculation_duration.count() << endl;
//    cout << "# solution: " << best_solution.daneX << endl;
//    cout << "# solution: " << best_solution.daneY << endl;

    return 0;
}









