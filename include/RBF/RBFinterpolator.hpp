#pragma once

#include <Eigen/Dense>
#include "interpolator.hpp"

/**
 * @brief Interpolator using Radial Basis Function (RBF) interpolation method.
 */
class RBFInterpolator : public Interpolator {
private:
    std::function<double(double, double)> rbfunction; /**< RBF function. */
    double r0; /**< Scale factor. */

public:
    /**
     * @brief Constructs an RBFInterpolator object.
     * 
     * @param rbfunction The RBF function.
     * @param r0 The scale factor.
     */
    RBFInterpolator(std::function<double(double, double)> rbfunction, double r0)
        : rbfunction(rbfunction), r0(r0) {}

    /**
     * @brief Interpolates a value based on given parameters and measurements using RBF method.
     * 
     * @param parametersFORinterp The points we want to interpolate.
     * @param parameters The matrix of parameters.
     * @param measurements The measurements corresponding to the parameters.
     * @return The interpolated value.
     */
    Eigen::VectorXd interpolate(const Eigen::MatrixXd& parametersFORinterp,
                                const Eigen::MatrixXd& parameters,
                                const Eigen::VectorXd& measurements) const override;
};