#pragma once

/**
 * @brief Collection of Radial Basis Function (RBF) functions.
 */
class RBFunctions {
public:
    /**
     * @brief Multiquadratic radial basis function.
     * 
     * @param r The radial distance.
     * @param r0 The scale factor.
     * @return The computed value.
     */
    static double multiquadratic(double r, double r0);

    /**
     * @brief Inverse Multiquadratic radial basis function.
     * 
     * @param r The radial distance.
     * @param r0 The scale factor.
     * @return The computed value.
     */
    static double inverseMultiquadratic(double r, double r0);

    /**
     * @brief Gaussian radial basis function.
     * 
     * @param r The radial distance.
     * @param r0 The scale factor.
     * @return The computed value.
     */
    static double gaussian(double r, double r0);

    /**
     * @brief Thin Plate Spline radial basis function.
     * 
     * @param r The radial distance.
     * @param r0 The scale factor.
     * @return The computed value.
     */
    static double thinPlateSpline(double r, double r0);
};