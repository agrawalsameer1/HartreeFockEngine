#include <cmath>

#pragma once

class Coordinate {
    public:
        double x, y, z;

        Coordinate();
        Coordinate(double X, double Y, double Z);
        Coordinate(double* coords);

        double* getCoords();
        double dist(Coordinate coord);
};