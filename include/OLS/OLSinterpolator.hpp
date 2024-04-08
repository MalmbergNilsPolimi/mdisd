#pragma once

#include <Eigen/Dense>
#include "interpolator.hpp"

/**
 * @brief Interpolator using Radial Basis Function (RBF) interpolation method.
 */
class OLSInterpolator : public Interpolator {

public:
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