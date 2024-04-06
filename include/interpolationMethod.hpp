#ifndef INTERPOLATIONMETHOD_HPP
#define INTERPOLATIONMETHOD_HPP

#include <Eigen/Core>

class InterpolationMethod {
public:
    virtual double interpolate(const Eigen::VectorXd& parametersFORinterp, const Eigen::MatrixXd& parameters, const Eigen::VectorXd& measurements) const = 0;
};

#endif // INTERPOLATIONMETHOD_HPP