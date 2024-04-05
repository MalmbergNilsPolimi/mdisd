#include "RBFinterpolation.hpp"

RBFInterpolation::RBFInterpolation(std::function<double(double, double)> rbfunction) : rbfunction(rbfunction) {

}

double RBFInterpolation::rbMultiquadratic(double r, double r0) {
    
}

double RBFInterpolation::rbInverseMultiquadratic(double r, double r0) {

}

double RBFInterpolation::rbThinPlateSpline(double r, double r0) {

}

double RBFInterpolation::rbGaussian(double r, double r0) {

}