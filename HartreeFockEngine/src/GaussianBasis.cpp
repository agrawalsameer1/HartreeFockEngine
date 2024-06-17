#include "GaussianBasis.h"

GaussianBasis::GaussianBasis() {
    a = 0;
    coeff = 0;
    nucCharge = 0;
}

GaussianBasis::GaussianBasis(double aVal, double cVal, int nuclearCharges) {
    a = aVal;
    coeff = cVal;
    nucCharge = nuclearCharges;
}

double GaussianBasis::calc(double r) {
    return coeff*exp(-1*a*r*r);
}

double GaussianBasis::overlap(GaussianBasis b) {
    return coeff*b.coeff*pow((M_PI/(a+b.a)), 1.5);
}

double GaussianBasis::KineticEnergy(GaussianBasis b) {
    return coeff*b.coeff*(3*a*b.a*pow(M_PI, 1.5))/pow((a+b.a), 2.5);
}

double GaussianBasis::Potential(GaussianBasis b) {
    return coeff*b.coeff*(-1*nucCharge*2*M_PI)/(a+b.a);
}

double GaussianBasis::twoEInt(GaussianBasis b, GaussianBasis c, GaussianBasis d) {
    return coeff*b.coeff*c.coeff*d.coeff*(2*pow(M_PI, 2.5))/((a+b.a)*(c.a+d.a)*sqrt(a+b.a+c.a+d.a));
}