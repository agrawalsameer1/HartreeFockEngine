#include "HartreeFock.h"

HartreeFock::HartreeFock(Atom at) {
    atom = Atom(at.numBasisFuncs, at.BasisFuncs, at.nucleiCoords, at.nucleiCharge, at.numOrbitals, at.numElectrons);
    S = Eigen::MatrixXd::Zero(atom.numBasisFuncs, atom.numBasisFuncs);
    T = Eigen::MatrixXd::Zero(atom.numBasisFuncs, atom.numBasisFuncs);
    A = Eigen::MatrixXd::Zero(atom.numBasisFuncs, atom.numBasisFuncs);
    H = Eigen::MatrixXd::Zero(atom.numBasisFuncs, atom.numBasisFuncs);
    F = Eigen::MatrixXd::Zero(atom.numBasisFuncs, atom.numBasisFuncs);
    P = Eigen::MatrixXd::Zero(atom.numBasisFuncs, atom.numBasisFuncs);
    OldP = Eigen::MatrixXd::Zero(atom.numBasisFuncs, atom.numBasisFuncs);
}

void HartreeFock::initMatrices() {
    for (int i = 0; i < atom.numBasisFuncs; i++) {
        for (int j = 0; j < atom.numBasisFuncs; j++) {
            S(i, j) = atom.BasisFuncs[i].overlap(atom.BasisFuncs[j]);
            T(i, j) = atom.BasisFuncs[i].KineticEnergy(atom.BasisFuncs[j]);
            A(i, j) = atom.BasisFuncs[i].Potential(atom.BasisFuncs[j]);
        }
    }

    H = T+A;

    // Initializes C matrix as though two-electron integrals did not exist
    if (atom.numElectrons > 1) {
        solveEigenvals(H);
    }
}

void HartreeFock::updateDensity() {
    for (int i = 0; i < atom.numBasisFuncs; i++) {
        for (int j = 0; j < atom.numBasisFuncs; j++) {
            OldP(i,j) = P(i,j);
        }
    }

    for (int i = 0; i < atom.numBasisFuncs; i++) {
        for (int j = 0; j < atom.numBasisFuncs; j++) {
            double sum = 0;

            // sum over all occupied orbitals, NOT BASIS FUNCTIONS in coefficient matrix
            for (int elec = 0; elec < atom.numOrbitals; elec++) {
                sum += 2*C(i, elec)*C(j, elec);
            }

            P(i, j) = sum;
        }
    }
}

void HartreeFock::makeFockMatrix() {
    for (int i = 0; i < atom.numBasisFuncs; i++) {
        for (int j = 0; j < atom.numBasisFuncs; j++) {

            F(i,j) = H(i,j);

            for (int k = 0; k < atom.numBasisFuncs; k++) {
                for (int l = 0; l < atom.numBasisFuncs; l++) {
                    F(i,j) += P(k,l)*(atom.BasisFuncs[i].twoEInt(atom.BasisFuncs[j], atom.BasisFuncs[k], atom.BasisFuncs[l])-0.5*atom.BasisFuncs[i].twoEInt(atom.BasisFuncs[k], atom.BasisFuncs[j], atom.BasisFuncs[l]));
                }
            }
        }
    }
}

MatrixXd HartreeFock::makeDiagonal(MatrixXd toTransform) {
    for (int i = 0; i < toTransform.rows(); i++) {
        for (int j = 0; j < toTransform.cols(); j++) {
            if (i != j) {
                toTransform(i, j) = 0;
            }
        }
    }

    return toTransform;
}

double HartreeFock::calcDelta(MatrixXd oldMatrix, MatrixXd newMatrix) {
    double delta = 0;

    for (int i = 0; i < oldMatrix.rows(); i++) {
        for (int j = 0; j < oldMatrix.cols(); j++) {
            delta += ((newMatrix(i,j)-oldMatrix(i,j))*(newMatrix(i,j)-oldMatrix(i,j)));
        }
    }

    return sqrt(delta);
}

MatrixXd HartreeFock::transformMatrix(MatrixXd toTransform, MatrixXd X) {
    return X.transpose()*toTransform*X;
}

MatrixXd HartreeFock::calculateU() {
    EigenSolver<MatrixXd> es(S);
    MatrixXd U = es.eigenvectors().real();
    return U;
}

MatrixXd HartreeFock::calculateV(MatrixXd U) {
    MatrixXd sqrt_U = U.array().sqrt();
    MatrixXd V = sqrt_U.matrix().inverse();
    return V;
}

double HartreeFock::solveEigenvals(MatrixXd toSolve) {
    MatrixXd U = calculateU();
    MatrixXd s = makeDiagonal(transformMatrix(S, U));
    MatrixXd V = U*calculateV(s);
    MatrixXd Hprime = transformMatrix(toSolve, V);

    EigenSolver<MatrixXd> es2(Hprime);
    VectorXd eigenvalues = es2.eigenvalues().real();
    MatrixXd eigenvectors = es2.eigenvectors().real();

    std::vector<std::pair<double, VectorXd>> eigenPairs;

    for (int i = 0; i < eigenvalues.rows(); i++) {
        VectorXd column = eigenvectors.col(i);
        eigenPairs.push_back(std::make_pair(eigenvalues(i), column));
    }

    std::sort(eigenPairs.begin(), eigenPairs.end(), sortEigenVectors());

    // C matrix should have the eigenvectors corresponding to the lowest energy occupied orbitals, or the total electrons
    // C matrix dimensions should be basisFuncs x electrons
    C = Eigen::MatrixXd::Zero(atom.numBasisFuncs, atom.numOrbitals);

    // Assign C to first N/2 eigenvectors of Hprime sorted by eigenvalue, where N is the number of electrons 
    for (int i = 0; i < atom.numOrbitals; i++) {
        VectorXd colVector = eigenPairs[i].second;
        for (int j = 0; j < atom.numBasisFuncs; j++) {
            // ith column of C (second index in parenthesis indexing) should have ith eigenvector of Hprime
            C(j,i) = colVector(j);
        }
    }

    // Backtransform C
    C = V*C;

    // Normalize C vectors
    for (int i = 0; i < atom.numOrbitals; i++) {
        for (int j = 0; j < atom.numBasisFuncs; j++) {
            C(j,i) /= sqrt((C.col(i).transpose()*S*C.col(i)).value());
        }
    }

    return eigenPairs[0].first;
}

double HartreeFock::calcEnergyFromDensity() {
    double finalEnergy = 0;

    for (int i = 0; i < atom.numBasisFuncs; i++) {
        for (int j = 0; j < atom.numBasisFuncs; j++) {
            finalEnergy += 0.5*P(i,j)*(H(i,j) + F(i,j));
        }
    }

    return finalEnergy;
}

void HartreeFock::SCF(bool logging) {
    double deltaDensity = 1;
    int counter = 0;

    initMatrices();

    updateDensity();

    if (logging) {
        std::cout << "As of iteration # " << counter << ": energy is " << calcEnergyFromDensity() << std::setprecision(10) << " hartrees\n";
    }

    while (deltaDensity > 0.0000001) {
        makeFockMatrix();

        solveEigenvals(F);
        
        updateDensity();

        deltaDensity = calcDelta(P, OldP);

        counter++;

        if (logging) {
            std::cout << "As of iteration # " << counter << ": energy is " << calcEnergyFromDensity() << std::setprecision(10) << " hartrees\n";
        }
    }
    std::cout << "\n";
}

double HartreeFock::finalSolveEnergy(bool logging) {
    if (atom.numElectrons == 1) {
        initMatrices();
        return solveEigenvals(H);
    }

    else {
        SCF(logging);
        return calcEnergyFromDensity();
    }   
}