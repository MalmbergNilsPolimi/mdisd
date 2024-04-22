#pragma once

#include <Eigen/Dense>

class Rescaling {
public:

    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> meanNormalization(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2 = nullptr);
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> minMaxNormalization(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2 = nullptr);
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> zScoreNormalization(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2 = nullptr);

private:
    Eigen::MatrixXd combineData(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2 = nullptr);
    Eigen::VectorXd computeColumnMeans(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2 = nullptr);
    Eigen::VectorXd computeColumnStdDevs(const Eigen::MatrixXd& data1, const Eigen::VectorXd& means, const Eigen::MatrixXd* data2 = nullptr);
    Eigen::VectorXd computeColumnMin(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2 = nullptr);
    Eigen::VectorXd computeColumnMax(const Eigen::MatrixXd& data1, const Eigen::MatrixXd* data2 = nullptr);
};
