#include "Coordinate.h"
#include "GaussianBasis.h"

#pragma once

class Atom {
    public:
        int numBasisFuncs;
        GaussianBasis* BasisFuncs;

        Coordinate nucleiCoords;
        int nucleiCharge;
        int numOrbitals;
        int numElectrons;

        Atom();
        Atom(int numBasisFunctions, double* a, double* coeffs, Coordinate nucleiCoordinates, int nucleiQ, int numOrbs, int numElecs);
        Atom(int numBasisFunctions, GaussianBasis* bases, Coordinate nucleiCoordinates, int nucleiQ, int numOrbs, int numElecs);
        double basisFunction(double a);
};