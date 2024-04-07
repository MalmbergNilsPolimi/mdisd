#include <cmath>
#include "RBFunctions.hpp"


double RBFunctions::multiquadratic(double r, double r0) {
    return sqrt(r*r + r0*r0);
}

double RBFunctions::inverseMultiquadratic(double r, double r0) {
    return 1.0 / sqrt(r*r + r0*r0);
}

double RBFunctions::gaussian(double r, double r0) {
    return exp(-0.5 * r * r / (r0 * r0));
}

double RBFunctions::thinPlateSpline(double r, double r0) {
    double res{};

    if (r != 0) {
        res = r*r * log(r / r0);
    }

    return res;
}