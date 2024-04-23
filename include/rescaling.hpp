#pragma once

#include <Eigen/Dense>

class Rescaling {
public:

    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> meanNormalization(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2 = nullptr);
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> minMaxNormalization(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2 = nullptr);
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> zScoreNormalization(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2 = nullptr);

private:
    Eigen::VectorXd computeColumnMeans(const Eigen::MatrixXd& data);
    Eigen::VectorXd computeColumnStdDevs(const Eigen::MatrixXd& data, const Eigen::VectorXd& means);
    Eigen::VectorXd computeColumnMin(const Eigen::MatrixXd& data);
    Eigen::VectorXd computeColumnMax(const Eigen::MatrixXd& data);
};
