#include <iostream>
#include "HartreeFock.h"
#include "Schrodinger.h"

#define BERYLLIUM // Can be HYDROGEN, HELIUM, or BERYLLIUM

#ifdef HYDROGEN
int main() {
    const int numBasisFuncsPerOrbital = 6;
    const int numElectrons = 1;
    const int numOrbitals = (numElectrons+1)/2;

    const int numBasisFuncs = numBasisFuncsPerOrbital*numElectrons;

    const int nucCharges = 1;

    double aCoefficients[numBasisFuncs] = {35.52322122, 6.513143725, 1.822142904, 0.6259552659, 0.2430767471, 0.1001124280};
    double cCoefficients[numBasisFuncs] = {0.009163596281, 0.04936149294, 0.1685383049, 0.3705627997, 0.4164915298, 0.1303340841};

    Coordinate nucCoords = Coordinate(0, 0, 0);

    GaussianBasis* gauss = new GaussianBasis[numBasisFuncs];

    for (int i = 0; i < numBasisFuncs; i++) {
        gauss[i] = GaussianBasis(aCoefficients[i], cCoefficients[i], nucCharges);
    }

    Atom at = Atom(numBasisFuncs, gauss, nucCoords, nucCharges, numOrbitals, numElectrons);
    HartreeFock hfSolver = HartreeFock(at);

    double finalEnergy = hfSolver.finalSolveEnergy();
    std::cout << "ground state energy: " << std::setprecision(10) << finalEnergy << " hartrees\n";
    std::cout << "ground state energy: " << std::setprecision(10) << finalEnergy*27.2114 << " eV\n";

    Schrodinger sch = Schrodinger(hfSolver);
    std::cout << "atomic radius: " << sch.calcMaxPoint()*52.9177 << " picometers\n";
}
#endif

#ifdef HELIUM
int main() {
    const int numBasisFuncsPerOrbital = 6;
    const int numElectrons = 2;
    const int numOrbitals = (numElectrons+1)/2;

    const int numBasisFuncs = numBasisFuncsPerOrbital*numOrbitals;

    const int nucCharges = 2;

    double aCoefficients[numBasisFuncs] = {65.98456824, 12.09819836, 3.384639924, 1.162715163, 0.4515163224, 0.1859593559};
    double cCoefficients[numBasisFuncs] = {0.009163596281, 0.04936149294, 0.1685383049, 0.3705627997, 0.4164915298, 0.1303340841};

    Coordinate nucCoords = Coordinate(0, 0, 0);

    GaussianBasis* gauss = new GaussianBasis[numBasisFuncs];

    for (int i = 0; i < numBasisFuncs; i++) {
        gauss[i] = GaussianBasis(aCoefficients[i], cCoefficients[i], nucCharges);
    }

    Atom at = Atom(numBasisFuncs, gauss, nucCoords, nucCharges, numOrbitals, numElectrons);
    HartreeFock hfSolver = HartreeFock(at);
    double finalEnergy = hfSolver.finalSolveEnergy();
    std::cout << "ground state energy: " << std::setprecision(10) << finalEnergy << " hartrees\n";
    std::cout << "ground state energy: " << std::setprecision(10) << finalEnergy*27.2114 << " eV\n";

    Schrodinger sch = Schrodinger(hfSolver);
    std::cout << "atomic radius: " << sch.calcMaxPoint()*52.9177 << " picometers\n";
}
#endif

#ifdef BERYLLIUM
int main() {
    const int numBasisFuncsPerOrbital = 6;
    const int numElectrons = 4;

    const int numOrbitals = (numElectrons+1)/2;
    const int numBasisFuncs = numBasisFuncsPerOrbital*numOrbitals;

    const int nucCharges = 4;

    double aCoefficients[numBasisFuncs] = {312.8704937, 57.36446253, 16.04850940, 5.513096119, 2.140896553, 0.8817394283, 13.63324744, 2.698375464, 0.8386530829, 0.3226600698, 0.1401314882, 0.06423251387};
    double cCoefficients[numBasisFuncs] = {0.009163596281, 0.04936149294, 0.1685383049, 0.3705627997, 0.4164915298, 0.1303340841, -0.01325278809, -0.04699171014, -0.03378537151, 0.2502417861, 0.5951172526, 0.2407061763};

    Coordinate nucCoords = Coordinate(0, 0, 0);

    GaussianBasis* gauss = new GaussianBasis[numBasisFuncs];

    for (int i = 0; i < numBasisFuncs; i++) {
        gauss[i] = GaussianBasis(aCoefficients[i], cCoefficients[i], nucCharges);
    }

    Atom at = Atom(numBasisFuncs, gauss, nucCoords, nucCharges, numOrbitals, numElectrons);
    HartreeFock hfSolver = HartreeFock(at);

    double finalEnergy = hfSolver.finalSolveEnergy();
    std::cout << "ground state energy: " << std::setprecision(10) << finalEnergy << " hartrees\n";
    std::cout << "ground state energy: " << std::setprecision(10) << finalEnergy*27.2114 << " eV\n";

    Schrodinger sch = Schrodinger(hfSolver);
    std::cout << "atomic radius: " << sch.calcMaxPoint()*52.9177 << " picometers\n";
}
#endif

