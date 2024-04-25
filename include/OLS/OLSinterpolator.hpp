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
     * @param parametersFORinterp Th matrix containing the points to interpolate.
     * @param parameters The matrix containing the known parameters.
     * @param measurements The measurements corresponding to the known parameters.
     * @param regression Optional vector where the weights will be stored.
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