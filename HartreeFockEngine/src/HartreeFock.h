#include "Atom.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <iomanip>

#pragma once

using namespace Eigen;

class HartreeFock {
    public:
        struct sortEigenVectors {
            bool operator()(const std::pair<double,VectorXd> &left, const std::pair<double,VectorXd> &right) {
                return left.first < right.first;
            }
        };

        Atom atom;

        // Overlap matrix
        MatrixXd S;

        // Kinetic Energy Matrix
        MatrixXd T;

        // Nuclear Potential Matrix
        MatrixXd A;

        // Core Matrix = T + A
        MatrixXd H;

        // Density Matrix
        MatrixXd P;
        MatrixXd OldP; // for convergence testing

        // Fock Matrix
        MatrixXd F;

        // C Matrix
        MatrixXd C;

        HartreeFock(Atom at);
        
        MatrixXd makeDiagonal(MatrixXd toTransform);
        double calcDelta(MatrixXd oldMatrix, MatrixXd newMatrix);

        MatrixXd transformMatrix(MatrixXd toTransform, MatrixXd X);
        MatrixXd calculateU();
        MatrixXd calculateV(MatrixXd U);

        void initMatrices();
        void updateDensity();
        void makeFockMatrix();
        double solveEigenvals(MatrixXd toSolve);
        double calcEnergyFromDensity();

        void SCF(bool logging = true);
        double finalSolveEnergy();
};