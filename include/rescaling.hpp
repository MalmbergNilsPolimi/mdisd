#pragma once

#include <Eigen/Dense>

/**
 * @brief Collection of rescaling methods for data preprocessing.
 */
class Rescaling {
public:
    /**
     * @brief Mean normalization of data.
     * 
     * @param data1 Matrix to be normalized, coefficients will be compute from this one.
     * @param data2 Matrix to be normalized, coefficients will be compute from data1.
     * @return Normalized data1 and data2.
     */
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> meanNormalization(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2 = nullptr);
    
    /**
     * @brief Min-Max normalization of data.
     * 
     * @param data1 Matrix to be normalized, coefficients will be compute from this one.
     * @param data2 Matrix to be normalized, coefficients will be compute from data1.
     * @return Normalized data1 and data2.
     */
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> minMaxNormalization(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2 = nullptr);
    
    /**
     * @brief Z-score normalization/standardization of data.
     * 
     * @param data1 Matrix to be normalized, coefficients will be compute from this one.
     * @param data2 Matrix to be normalized, coefficients will be compute from data1.
     * @return Normalized data1 and data2.
     */
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> zScoreNormalization(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2 = nullptr);

private:
    /**
     * @brief Compute the mean of each column of an Eigen::MatrixXd.
     * 
     * @param data The matrix used for means computations.
     * @return A vector containing the means of each column of data.
     */
    Eigen::VectorXd computeColumnMeans(const Eigen::MatrixXd& data);

    /**
     * @brief Compute the standard deviation of each column of an Eigen::MatrixXd.
     * 
     * @param data The matrix used for standard deviations computations.
     * @return A vector containing the standard deviations of each column of data.
     */
    Eigen::VectorXd computeColumnStdDevs(const Eigen::MatrixXd& data, const Eigen::VectorXd& means);
    
    /**
     * @brief Find the minimum element of each column of an Eigen::MatrixXd.
     * 
     * @param data The matrix used for minimums computations.
     * @return A vector containing the minimums of each column of data.
     */
    Eigen::VectorXd computeColumnMin(const Eigen::MatrixXd& data);

    /**
     * @brief Find the maximum element of each column of an Eigen::MatrixXd.
     * 
     * @param data The matrix used for maximums computations.
     * @return A vector containing the maximums of each column of data.
     */
    Eigen::VectorXd computeColumnMax(const Eigen::MatrixXd& data);
};