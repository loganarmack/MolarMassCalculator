#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include "molarmasscalculator.h"

using namespace std;

int main() {
    MolarMassCalculator m;

    bool ignoreEndline = false;
    while (1) {
        string inputMolecule;
        if (ignoreEndline) cin.ignore(numeric_limits<streamsize>::max(), '\n');
        cout << "Enter reaction or molecule: " << endl;
        getline(cin, inputMolecule);

        try {
            if (inputMolecule.find('+') != string::npos || inputMolecule.find("-->") != string::npos) {
                m.solveReaction(inputMolecule);
                ignoreEndline = true;
            }
            else {
                m.solveMassWithPercent(inputMolecule);
                ignoreEndline = false;
            }
        }
        catch (const char * msg) {
            cout << msg << endl;
        }
        catch (std::invalid_argument &e) {
            cout << "Invalid mass!" << endl;
            ignoreEndline = true;
        }
        cout << endl;
    }

    return 0;
}



