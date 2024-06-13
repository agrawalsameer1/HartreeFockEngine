#include <math.h> 

#pragma once

enum orbitalType {
    s, p, d, f
};

class GaussianBasis {
    public:
        double a;
        double coeff;
        int nucCharge;

        GaussianBasis();
        GaussianBasis(double aVal, double cVal, int nuclearCharges);
        
        double calc(double r);

        double overlap(GaussianBasis b = GaussianBasis(-1, -1, -1));

        double KineticEnergy(GaussianBasis b = GaussianBasis(-1, -1, -1));

        double Potential(GaussianBasis b = GaussianBasis(-1, -1, -1));

        double twoEInt(GaussianBasis b, GaussianBasis c, GaussianBasis d);
};