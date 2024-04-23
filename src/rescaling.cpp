#include "rescaling.hpp"
#include <algorithm>

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> Rescaling::meanNormalization(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2) {
    Eigen::VectorXd means = computeColumnMeans(data1);
    Eigen::VectorXd min_values = computeColumnMin(data1);
    Eigen::VectorXd max_values = computeColumnMax(data1);
    
    Eigen::MatrixXd normalized_data1 = (data1.rowwise() - means.transpose()).array().rowwise() / (max_values - min_values).transpose().array();

    Eigen::MatrixXd normalized_data2;
    if (data2) {
        normalized_data2 = (data2->rowwise() - means.transpose()).array().rowwise() / (max_values - min_values).transpose().array();
    }

    return {normalized_data1, normalized_data2};
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> Rescaling::minMaxNormalization(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2) {
    Eigen::VectorXd min_values = computeColumnMin(data1);
    Eigen::VectorXd max_values = computeColumnMax(data1);

    Eigen::MatrixXd normalized_data1 = (data1.rowwise() - min_values.transpose()).array().rowwise() / (max_values - min_values).transpose().array();
    Eigen::MatrixXd normalized_data2;
    if (data2) {
        normalized_data2 = (data2->rowwise() - min_values.transpose()).array().rowwise() / (max_values - min_values).transpose().array();
    }

    return {normalized_data1, normalized_data2};
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> Rescaling::zScoreNormalization(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2) {
    Eigen::VectorXd means = computeColumnMeans(data1);
    Eigen::VectorXd std_devs = computeColumnStdDevs(data1, means);

    Eigen::MatrixXd normalized_data1 = ((data1.rowwise() - means.transpose()).array().rowwise() / std_devs.transpose().array()).eval();
    Eigen::MatrixXd normalized_data2;
    if (data2) {
        normalized_data2 = ((data2->rowwise() - means.transpose()).array().rowwise() / std_devs.transpose().array()).eval();
    }

    return {normalized_data1, normalized_data2};
}





Eigen::VectorXd Rescaling::computeColumnMeans(const Eigen::MatrixXd& data) {
    return data.colwise().mean();
}

Eigen::VectorXd Rescaling::computeColumnStdDevs(const Eigen::MatrixXd& data, const Eigen::VectorXd& means) {
    Eigen::VectorXd squared_diffs = (data.rowwise() - means.transpose()).array().square().colwise().sum();
    return (squared_diffs / static_cast<double>(data.rows())).array().sqrt();
}

Eigen::VectorXd Rescaling::computeColumnMin(const Eigen::MatrixXd& data) {
    return data.colwise().minCoeff();
}

Eigen::VectorXd Rescaling::computeColumnMax(const Eigen::MatrixXd& data) {
    return data.colwise().maxCoeff();
}
