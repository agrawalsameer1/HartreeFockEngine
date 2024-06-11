#include "GaussianBasis.h"

GaussianBasis::GaussianBasis() {
    a = 0;
    nucCharge = 0;
}

GaussianBasis::GaussianBasis(double aVal, int nuclearCharges) {
    a = aVal;
    nucCharge = nuclearCharges;
}

double GaussianBasis::calc(double r) {
    return exp(-1*a*r*r);
}

double GaussianBasis::overlap(GaussianBasis b) {
    if (b.a == -1) {
        return 1;
    }
    else {
        return pow((M_PI/(a+b.a)), 1.5);
    }
}

double GaussianBasis::generalOverlap(GaussianBasis b) {
    if (b.a == -1) {
        return 1;
    }
    else {
        return pow((M_PI/(a+b.a)), 1.5);
    }
}

double GaussianBasis::KineticEnergy(GaussianBasis b) {
    if (b.a == -1) {
        return (3*a*a*pow(M_PI, 1.5))/pow((a+a), 2.5);
    }
    else {
        return (3*a*b.a*pow(M_PI, 1.5))/pow((a+b.a), 2.5);
    }
}

double GaussianBasis::Potential(GaussianBasis b) {
    if (b.a == -1) {
        return (-1*nucCharge*2*M_PI)/(a+a);
    }
    else {
        return (-1*nucCharge*2*M_PI)/(a+b.a);
    }
}

double GaussianBasis::twoEInt(GaussianBasis b, GaussianBasis c, GaussianBasis d) {
    return (2*pow(M_PI, 2.5))/((a+b.a)*(c.a+d.a)*sqrt(a+b.a+c.a+d.a));
}