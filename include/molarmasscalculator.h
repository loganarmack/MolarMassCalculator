#ifndef MOLARMASSCALCULATOR_H
#define MOLARMASSCALCULATOR_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cstdlib>

using namespace std;

class MolarMassCalculator
{
    public:
        MolarMassCalculator(string elementFile = "elementMass.csv");

        void solveMassWithPercent(string s);
        void solveReaction(string s);

    private:
        void SkipBOM(ifstream &in);
        static bool alphaSort(pair<string, double> p1, pair<string, double> p2);
        int binarySearch(vector<pair<string,double>> elements, string value);
        int findLimitingReactant(vector<pair<string, int>> reactants, vector<double> mass);
        double getReactionMass(pair<string, int> compound, pair<string, int> limitingCompound, double limitingCompoundMass);
        double solveMass(string inputMolecule);

        vector<pair<string, double>> elements;
        vector<string> molecules;

};

#endif // MOLARMASSCALCULATOR_H
