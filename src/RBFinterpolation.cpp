#include <cmath>
#include "RBFinterpolation.hpp"

RBFInterpolation::RBFInterpolation(std::function<double(double, double)> rbfunction) : rbfunction(rbfunction) {

}

double RBFInterpolation::rbMultiquadratic(double r, double r0) {
    return sqrt(r*r + r0*r0);
}

double RBFInterpolation::rbInverseMultiquadratic(double r, double r0) {
    return 1.0 / sqrt(r*r + r0*r0);
}

double RBFInterpolation::rbThinPlateSpline(double r, double r0) {
    double res{};

    if (r != 0) {
        res = r*r * log(r / r0);
    }

    return res;
}

double RBFInterpolation::rbGaussian(double r, double r0) {
    return exp(-0.5 * r * r / (r0 * r0));
}