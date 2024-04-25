#include <algorithm>
#include "rescaling.hpp"


std::pair<Eigen::MatrixXd, Eigen::MatrixXd> Rescaling::meanNormalization(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2) {
    // Mean normalization of data1 and data2 using data1 based coefficients.
    Eigen::VectorXd means = computeColumnMeans(data1);
    Eigen::VectorXd min_values = computeColumnMin(data1);
    Eigen::VectorXd max_values = computeColumnMax(data1);
    
    // Mean normalization.
    Eigen::MatrixXd normalized_data1 = (data1.rowwise() - means.transpose()).array().rowwise() / (max_values - min_values).transpose().array();

    Eigen::MatrixXd normalized_data2;
    if (data2) {
        // Mean normalization.
        normalized_data2 = (data2->rowwise() - means.transpose()).array().rowwise() / (max_values - min_values).transpose().array();
    }
    return {normalized_data1, normalized_data2};
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> Rescaling::minMaxNormalization(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2) {
    // Min-Max normalization of data1 and data2 using data1 based coefficients.
    Eigen::VectorXd min_values = computeColumnMin(data1);
    Eigen::VectorXd max_values = computeColumnMax(data1);

    // Min-Max normalization.
    Eigen::MatrixXd normalized_data1 = (data1.rowwise() - min_values.transpose()).array().rowwise() / (max_values - min_values).transpose().array();
    
    Eigen::MatrixXd normalized_data2;
    if (data2) {
        // Min-Max normalization.
        normalized_data2 = (data2->rowwise() - min_values.transpose()).array().rowwise() / (max_values - min_values).transpose().array();
    }
    return {normalized_data1, normalized_data2};
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> Rescaling::zScoreNormalization(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2) {
    // Z-score normalization/standardization of data1 and data2 using data1 based coefficients.
    Eigen::VectorXd means = computeColumnMeans(data1);
    Eigen::VectorXd std_devs = computeColumnStdDevs(data1, means);

    // Z-score normalization.
    Eigen::MatrixXd normalized_data1 = ((data1.rowwise() - means.transpose()).array().rowwise() / std_devs.transpose().array()).eval();
    
    Eigen::MatrixXd normalized_data2;
    if (data2) {
        // Z-score normalization.
        normalized_data2 = ((data2->rowwise() - means.transpose()).array().rowwise() / std_devs.transpose().array()).eval();
    }
    return {normalized_data1, normalized_data2};
}



Eigen::VectorXd Rescaling::computeColumnMeans(const Eigen::MatrixXd& data) {
    // Compute the mean of each column of data.
    return data.colwise().mean();
}

Eigen::VectorXd Rescaling::computeColumnStdDevs(const Eigen::MatrixXd& data, const Eigen::VectorXd& means) {
    // Compute the standard deviation of each column of data.
    Eigen::VectorXd squared_diffs = (data.rowwise() - means.transpose()).array().square().colwise().sum();
    return (squared_diffs / static_cast<double>(data.rows())).array().sqrt();
}

Eigen::VectorXd Rescaling::computeColumnMin(const Eigen::MatrixXd& data) {
    // Compute the minimum value of each column of data.
    return data.colwise().minCoeff();
}

Eigen::VectorXd Rescaling::computeColumnMax(const Eigen::MatrixXd& data) {
    // Compute the maximum value of each column of data.
    return data.colwise().maxCoeff();
}