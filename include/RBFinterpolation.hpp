#ifndef RBFINTERPOLATION_HPP
#define RBFINTERPOLATION_HPP

#include "interpolationMethod.hpp"
#include <functional>

class RBFInterpolation : public InterpolationMethod {

private:
    std::function<double(double, double)> rbfunction; // Radial basis function

public:
    RBFInterpolation(std::function<double(double, double)> rbfunction);
    
    double interpolate(const std::vector<double>& parametersFORinterp, 
                       const std::vector<std::vector<double>>& parameters,
                       const std::vector<double>& measurements)
                       const override;

    static double rbMultiquadratic(double r, double r0);
    static double rbInverseMultiquadratic(double r, double r0);
    static double rbThinPlateSpline(double r, double r0);
    static double rbGaussian(double r, double r0);

};

#endif // RBFINTERPOLATION_HPP
