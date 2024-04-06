#include <cmath>
#include "RBFinterpolation.hpp"
#include "solveSVD.hpp"

RBFInterpolation::RBFInterpolation(std::function<double(double, double)> rbfunction) : rbfunction(rbfunction) {
}

double RBFInterpolation::interpolate(const Eigen::VectorXd& parametersFORinterp, 
                                     const Eigen::MatrixXd& parameters,
                                     const Eigen::VectorXd& measurements) const {

    double res{};
    double r0{0.5};
    size_t num_measures{static_cast<size_t>(measurements.size())};

    // Computation of the weights
    Eigen::VectorXd weights(num_measures);
    Eigen::MatrixXd coeff(num_measures,num_measures);
    
    for (size_t i = 0; i < num_measures; i++)
    {
        for (size_t j = 0; i < num_measures; j++)
        {
            coeff(i,j)=rbfunction((parameters.row(i)-parameters.row(j)).norm(), r0);
        }
    }
    
    weights = solveSVD(coeff, measurements);

    // Computation of the interpolated value
    for (size_t i = 0; i < num_measures; ++i) {
        res += weights(i) * rbfunction((parametersFORinterp-parameters.row(i)).norm(), r0);
    }

    return res;
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