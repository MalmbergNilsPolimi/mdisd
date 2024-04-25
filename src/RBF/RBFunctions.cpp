#include <cmath>
#include <iostream>
#include "RBFunctions.hpp"


double RBFunctions::multiquadratic(double r, double r0) {
    // Multiquadratic radial basis function.
    return sqrt(r*r + r0*r0);
}

double RBFunctions::inverseMultiquadratic(double r, double r0) {
    // Inverse multiquadratic radial basis function.
    return 1. / sqrt(r*r + r0*r0);
}

double RBFunctions::gaussian(double r, double r0) {
    // Gaussian radial basis function.
    return exp(-0.5 * r * r / (r0 * r0));
}

double RBFunctions::thinPlateSpline(double r, double r0) {
    // Thin plate spline radial basis function.
    double res{};

    if (r != 0) {
        res = r*r * log(r / r0);
    } else {
        std::cerr << "Error: thinPlateSpline not defined for negative radial distances." << std::endl;
    }

    return res;
}