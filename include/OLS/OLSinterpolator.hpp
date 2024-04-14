#pragma once

#include <Eigen/Dense>
#include "interpolator.hpp"

/**
 * @brief Interpolator using Ordinary Least Squares (OLS) approximation method.
 */
class OLSInterpolator : public Interpolator {

public:
    /**
     * @brief Interpolates a value based on given parameters and measurements using OLS method.
     * 
     * @param parametersFORinterp The points we want to interpolate.
     * @param parameters The matrix of parameters.
     * @param measurements The measurements corresponding to the parameters.
     * @return The interpolated value.
     */
    Eigen::VectorXd interpolate(const Eigen::MatrixXd& parametersFORinterp,
                                const Eigen::MatrixXd& parameters,
                                const Eigen::VectorXd& measurements) const override;

    Eigen::VectorXd interpolate(const Eigen::MatrixXd& parametersFORinterp,
                            const Eigen::MatrixXd& parameters,
                            const Eigen::VectorXd& measurements,
                            Eigen::VectorXd* regression) const override;
};