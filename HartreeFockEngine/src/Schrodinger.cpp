#include "Schrodinger.h"

Schrodinger::Schrodinger(HartreeFock hf) {
    numBases = hf.atom.numBasisFuncs;

    bases = new GaussianBasis[numBases];

    for (int i = 0; i < numBases; i++) {
        bases[i] = GaussianBasis(hf.atom.BasisFuncs[i].a, hf.atom.nucleiCharge);
    }

    coeffs = new double[numBases];

    for (int i = 0; i < numBases; i++) {
        coeffs[i] = hf.C(i, 0);
    }

    funcVals = new double[numDivisions];

    for (int i = 0; i < numDivisions; i++) {
        double x = rangeToTest/numDivisions*i;
        double value = 0;

        for (int j = 0; j < numBases; j++) {
            value += (coeffs[j]*exp(-1*bases[j].a*x*x));
        }
        
        value = value*value;
        value *= (4*M_PI*(x*x));

        funcVals[i] = value;
    }
}

Schrodinger::~Schrodinger() {
    delete[] bases;
    delete[] coeffs;
    delete[] funcVals;
}

double Schrodinger::calcMaxPoint() {
    double* difFuncVals = new double[numDivisions-1];

    for (int i = 1; i < numDivisions; i++) {
        difFuncVals[i-1] = (funcVals[i]-funcVals[i-1])/(rangeToTest/numDivisions);
    }

    for (int i = 0; i < numDivisions - 100; i++) {
        if (abs(difFuncVals[i]) < (rangeToTest*10/numDivisions) && i > (1/(10*rangeToTest)*numDivisions) && abs(difFuncVals[i+1000]) > (rangeToTest*10/numDivisions)) {
            return (rangeToTest/numDivisions*i);
        }
    }
}