#include "molarmasscalculator.h"

using namespace std;

MolarMassCalculator::MolarMassCalculator(string elementFile) {
    //setup file input
    ifstream fin(elementFile);
    SkipBOM(fin);

    //load elements & mass into elements vector
    int l = 0;
    while (!fin.eof()) {
        l++;

        //take in input from csv
        string buff;
        getline(fin, buff);
        istringstream ss(buff);
        pair<string, double> temp;

        //element name
        string line;
        getline(ss, line, ',');

        //element mass
        temp.first = line;
        getline(ss, line, ',');
        try {
            temp.second = stof(line); //converts string representation of mass to float
        } catch (const invalid_argument e) {
            cout << "NaN found in file elementMass.csv line " << l << endl;
            e.what();
        }
        elements.push_back(temp);
    }

    sort(elements.begin(), elements.end(), alphaSort);
}

bool MolarMassCalculator::alphaSort(pair<string, double> p1, pair<string, double> p2) {
    for (int i = 0; i < 2; i++) {
        if (p1.first[i] > p2.first[i]) {
            return false;
        } else if (p1.first[i] < p2.first[i]) {
            return true;
        }
    }
    return true;
}

// Skips the Byte Order Mark (BOM) that defines UTF-8 in some text files.
void MolarMassCalculator::SkipBOM(ifstream &in) {
    char test[3] = {0};
    in.read(test, 3);
    if ((unsigned char)test[0] == 0xEF &&
            (unsigned char)test[1] == 0xBB &&
            (unsigned char)test[2] == 0xBF) {
        return;
    }
    in.seekg(0);
}

int MolarMassCalculator::binarySearch(vector<pair<string,double>> elements, string value) {
    int first = 0;
    int last = elements.size() - 1;
    int middle = (first + last) / 2;

    while (first <= last) {
        middle = (first + last) / 2;     // Calculate mid point
        if (elements[middle].first == value) {     // If value is found at mid
            return middle;
        } else if (elements[middle].first > value) // If value is in lower half
            last = middle - 1;
        else
            first = middle + 1;           // If value is in upper half
    }
    return -1;
}

void MolarMassCalculator::solveMassWithPercent(string s) {
    double molarMass = solveMass(s);
    cout << "The molar mass of " << s << " is " << molarMass << " g/mol." << endl;
}

//returns the molar mass of a given molecule
double MolarMassCalculator::solveMass(string inputMolecule) {
    string element = "";
    int elementNum = 0;
    double molarMass = 0;

    //loop through each character in input molecule
    for (int i = 0; i < inputMolecule.size(); i++) {
        if (element.size() == 0) {
            //look at new element
            if (inputMolecule[i] >= 'A' && inputMolecule[i] <= 'Z') {
                element += inputMolecule[i];
            }
            //if bracket, looks for molecule and recurses
            else if (inputMolecule[i] == '(') {
                int startIndex = i;
                int numOpenBrackets = 0;
                //tries to find closing bracket to match open bracket
                while (true) {
                    i++;
                    if (inputMolecule[i] == '(') {
                        numOpenBrackets++;
                    } else if (inputMolecule[i] == ')' && numOpenBrackets != 0) {
                        numOpenBrackets--;
                    } else if (inputMolecule[i] == ')') {
                        break;
                    }
                    if (i > inputMolecule.size() - 1) {
                        throw "Invalid compound!";
                    }
                }
                string subMolecule = inputMolecule.substr(startIndex + 1, i - startIndex - 1);
                int numSubmolecules = 0;

                //finds number after closing bracket
                i++; //skip past closing bracket
                while (inputMolecule[i] >= '0' && inputMolecule[i] <= '9') {
                    numSubmolecules *= 10;
                    numSubmolecules += inputMolecule[i] - '0';
                    i++;
                }
                //solves mass of compound within brackets
                molarMass += solveMass(subMolecule) * max(numSubmolecules, 1);
                i--;
            } else {
                throw "Invalid compound!";
            }
        } else if (element.size() == 1 || element.size() == 2) {
            //second letter of element
            if (inputMolecule[i] >= 'a' && inputMolecule[i] <= 'z' && element.size() == 1 && elementNum == 0) {
                element += inputMolecule[i];
            }
            //numbers after element
            else if (inputMolecule[i] >= '0' && inputMolecule[i] <= '9') {
                elementNum *= 10;
                elementNum += inputMolecule[i] - '0';
            } else if ((inputMolecule[i] >= 'A' && inputMolecule[i] <= 'Z') || inputMolecule[i] == '(') {
                int elementIndex = binarySearch(elements, element);
                if (elementIndex == -1) {
                    throw "Invalid element!";
                } else {
                    molarMass += elements[elementIndex].second * max(elementNum, 1);
                }
                element = "";
                elementNum = 0;
                if (inputMolecule[i] == '(') {
                    i--;
                } else {
                    element += inputMolecule[i];
                }
            } else {
                throw "Invalid compound!";
            }
        } else {
            throw "Invalid compound!";
        }
    }
    if (element != "") {
        int elementIndex = binarySearch(elements, element);
        if (elementIndex == -1) {
            throw "Invalid element!";
        } else {
            molarMass += elements[elementIndex].second * max(elementNum, 1);
        }
    }
    return molarMass;
}

void MolarMassCalculator::solveReaction(string s) {
    vector<pair<string, int>> reactants(1);
    vector<pair<string, int>> products(1);

    int currCompound = 0;
    bool inputtingCompound = false;
    bool isReactant = true;
    string temp = "";
    for (int i = 0; i < s.length(); i++) {
        if ((s[i] <= '9' && s[i] >= '0') && !inputtingCompound) {
            if (isReactant) {
                reactants[currCompound].second *= 10;
                reactants[currCompound].second += s[i] - '0';
            } else {
                products[currCompound].second *= 10;
                products[currCompound].second += s[i] - '0';
            }
        } else if (s[i] == '+') {
            if (isReactant && reactants[currCompound].first != "") {
                currCompound++;
                reactants.resize(reactants.size() + 1);
                inputtingCompound = false;
            } else if (!isReactant && products[currCompound].first != "") {
                currCompound++;
                products.resize(products.size() + 1);
                inputtingCompound = false;
            }
        } else if (i < s.length() - 3 && s[i] == '-' && s[i + 1] == '-' && s[i + 2] == '>') {
            currCompound = 0;
            isReactant = false;
            inputtingCompound = false;
            i += 2;
        } else if (s[i] == '+') {
            throw "Invalid syntax!";
        } else if (s[i] != ' ') {
            inputtingCompound = true;
            if (isReactant) {
                reactants[currCompound].first += s[i];
            } else {
                products[currCompound].first += s[i];
            }
        }
    }
    if (reactants.size() == 0) {
        throw "Reaction must have at least one reactant!";
    }
    if (products.size() == 0) {
        throw "Reaction must have at least one product!";
    }

    //set the number of each compound to at least 1
    for (int i = 0; i < reactants.size(); i++) {
        reactants[i].second = max(reactants[i].second, 1);
    }
    for (int i = 0; i < products.size(); i++) {
        products[i].second = max(products[i].second, 1);
    }

    vector<double> reactantMasses(reactants.size());
    vector<double> productMasses(products.size());

    cout << "Enter the masses of each reactant: " << endl;
    for (int i = 0; i < reactants.size(); i++) {
        string temp;
        cin >> temp;
        reactantMasses[i] = stof(temp);
    }

    int limitingReactant = findLimitingReactant(reactants, reactantMasses);

    cout << "The limiting reactant is " << reactants[limitingReactant].first << endl;

    for (int i = 0; i < reactants.size(); i++) {
        if (i != limitingReactant) {
            reactantMasses[i] = getReactionMass(reactants[i], reactants[limitingReactant], reactantMasses[limitingReactant]);
        }
    }
    for (int i = 0; i < products.size(); i++) {
        productMasses[i] = getReactionMass(products[i], reactants[limitingReactant], reactantMasses[limitingReactant]);
    }

    for (int i = 0; i < reactants.size(); i++) {
        cout << reactants[i].first << ": " << reactantMasses[i] << endl;
    }
    for (int i = 0; i < products.size(); i++) {
        cout << products[i].first << ": " << productMasses[i] << endl;
    }
}

int MolarMassCalculator::findLimitingReactant(vector<pair<string, int>> reactants, vector<double> mass) {
    pair<double, int> limitingReactant;
    limitingReactant.first = mass[0] / (solveMass(reactants[0].first) * reactants[0].second);
    limitingReactant.second = 0;

    for (int i = 1; i < reactants.size(); i++) {
        double relativeMass = mass[i] / (solveMass(reactants[i].first) * reactants[i].second);
        if (relativeMass < limitingReactant.first) {
            limitingReactant.first = relativeMass;
            limitingReactant.second = i;
        }
    }

    return limitingReactant.second;
}

double MolarMassCalculator::getReactionMass(pair<string, int> compound, pair<string, int> limitingCompound, double limitingCompoundMass) {
    return limitingCompoundMass / solveMass(limitingCompound.first) / limitingCompound.second * compound.second * solveMass(compound.first);
}
