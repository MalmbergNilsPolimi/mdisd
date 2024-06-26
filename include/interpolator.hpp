#pragma once

#include <Eigen/Dense>
#include <functional>

/**
 * @brief Base class for interpolation methods.
 */
class Interpolator {
public:
    /**
     * @brief Interpolates a value based on given parameters and measurements.
     * 
     * @param parametersFORinterp Th matrix containing the points to interpolate.
     * @param parameters The matrix containing the known parameters.
     * @param measurements The measurements corresponding to the known parameters.
     * @param regression Optional vector where the weights will be stored.
     * @return The interpolated value.
     */
    virtual Eigen::VectorXd interpolate(const Eigen::MatrixXd& parametersFORinterp,
                                        const Eigen::MatrixXd& parameters,
                                        const Eigen::VectorXd& measurements) const = 0;

    // Override because virtual functions can't have by default argurments as Eigen::VectorXd* regression = nullptr
    virtual Eigen::VectorXd interpolate(const Eigen::MatrixXd& parametersFORinterp,
                                        const Eigen::MatrixXd& parameters,
                                        const Eigen::VectorXd& measurements,
                                        Eigen::VectorXd* regression) const = 0;
    
    virtual ~Interpolator() = default;
};