#include "solveSVD.hpp"

Eigen::VectorXd solveSVD(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd x = svd.solve(b);
    return x;
}