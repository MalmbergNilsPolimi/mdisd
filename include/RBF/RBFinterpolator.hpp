#pragma once

#include <Eigen/Dense>
#include "interpolator.hpp"

/**
 * @brief Interpolator using Radial Basis Function (RBF) interpolation method.
 */
class RBFInterpolator : public Interpolator {
private:
    std::function<double(double, double)> rbfunction; /**< Radial basis function. */
    double r0; /**< Scale factor. */
    bool normalizeRBF; /**< Normalization of RBF if true. */
    bool polynomialRBF; /**< Add a linear polynomial term if true. */

public:
    /**
     * @brief Constructs an RBFInterpolator object.
     * 
     * @param rbfunction The radial basis function.
     * @param r0 The scale factor.
     * @param normalizeRBF Boolean to activate normalization of the RBF.
     * @param polynomialRBF Boolean to activate the additional polynomial term in RBF.
     */
    RBFInterpolator(std::function<double(double, double)> rbfunction, double r0,
                    bool normalizeRBF=false, bool polynomialRBF=false)
        : rbfunction(rbfunction), r0(r0), normalizeRBF(normalizeRBF), polynomialRBF(polynomialRBF) {}

    /**
     * @brief Interpolates a value based on given parameters and measurements using RBF method.
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