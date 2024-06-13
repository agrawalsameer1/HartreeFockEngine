#include "Atom.h"

Atom::Atom() {
    ;
}

Atom::Atom(int numBasisFunctions, double* a, double* coeffs, Coordinate nucleiCoordinates, int nucleiQ, int numOrbs, int numElecs) {
    numOrbitals = numOrbs;
    numElectrons = numElecs;
    numBasisFuncs = numBasisFunctions;
    BasisFuncs = new GaussianBasis[numBasisFuncs];

    for (int i = 0; i < numBasisFuncs; i++) {
        BasisFuncs[i] = GaussianBasis(a[i], coeffs[i], nucleiQ);
    }

    nucleiCoords = Coordinate(nucleiCoordinates);
    
    nucleiCharge = nucleiQ;
}

Atom::Atom(int numBasisFunctions, GaussianBasis* bases, Coordinate nucleiCoordinates, int nucleiQ, int numOrbs, int numElecs) {
    numOrbitals = numOrbs;
    numElectrons = numElecs;
    numBasisFuncs = numBasisFunctions;
    BasisFuncs = new GaussianBasis[numBasisFuncs];

    for (int i = 0; i < numBasisFuncs; i++) {
        BasisFuncs[i] = GaussianBasis(bases[i].a, bases[i].coeff, nucleiQ);
    }

    nucleiCoords = Coordinate(nucleiCoordinates);
    
    nucleiCharge = nucleiQ;
}