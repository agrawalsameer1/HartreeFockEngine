#include "Coordinate.h"

Coordinate::Coordinate() {
    x = 0;
    y = 0;
    z = 0;
}

Coordinate::Coordinate(double X, double Y, double Z) {
    x = X;
    y = Y;
    z = Z;
}

Coordinate::Coordinate(double* coords) {
    x = coords[0];
    y = coords[1];
    z = coords[2];
}

double* Coordinate::getCoords() {
    double* coords = new double[3];
    coords[0] = x;
    coords[1] = y;
    coords[2] = z;

    return coords;
}

double Coordinate::dist(Coordinate coord) {
    return sqrt((coord.x*coord.x)+(coord.y*coord.y)+(coord.z*coord.z));
}

