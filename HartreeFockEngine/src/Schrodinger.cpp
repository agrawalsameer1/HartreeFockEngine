#include "Schrodinger.h"

Schrodinger::Schrodinger(HartreeFock hf) {
    numBasesPerOrbital = hf.atom.numBasisFuncs/hf.atom.numOrbitals;
    numBases = numBasesPerOrbital;

    bases = new GaussianBasis[numBases];

    for (int i = (hf.atom.numBasisFuncs-numBases); i < hf.atom.numBasisFuncs; i++) {
        bases[i-(hf.atom.numBasisFuncs-numBases)] = GaussianBasis(hf.atom.BasisFuncs[i].a, hf.atom.BasisFuncs[i].coeff, hf.atom.nucleiCharge);
    }

    coeffs = new double[numBases];

    for (int i = (hf.atom.numBasisFuncs-numBases); i < hf.atom.numBasisFuncs; i++) {
        coeffs[i-(hf.atom.numBasisFuncs-numBases)] = hf.C(i, hf.atom.numOrbitals-1);
    }
    

    funcVals = new double[numDivisions];

    for (int i = 0; i < numDivisions; i++) {
        double x = rangeToTest/numDivisions*i;
        double value = 0;
        
        for (int j = 0; j < numBases; j++) {
            value += (coeffs[j]*bases[j].calc(x));
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
    double maxVals[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    int cnt = 0;
    for (int i = 1; i < numDivisions-1; i++) {
        if (funcVals[i-1] < funcVals[i] && funcVals[i] > funcVals[i+1]) {// && funcVals[i] > 0.01) {
            maxVals[cnt] = rangeToTest/numDivisions*i;
            cnt++;
        }
    }
    
    for (int i = 9; i >=0; i--) {
        if (maxVals[i] != 0) {
            return maxVals[i];
        }
    }
}