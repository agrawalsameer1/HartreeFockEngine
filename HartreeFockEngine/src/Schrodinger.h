#include "HartreeFock.h"
#include <iostream>
#include <Eigen/Core>

#pragma once

class Schrodinger {
    public:
        int numBases;
        int numBasesPerOrbital;
        GaussianBasis* bases;
        double* coeffs;

        double* funcVals;

        double rangeToTest = 10;
        const int numDivisions = 100000;

        Schrodinger(HartreeFock hf);
        ~Schrodinger();
        double calcMaxPoint();
};