#ifndef SOLVESVD_HPP
#define SOLVESVD_HPP

#include <Eigen/Dense>

Eigen::VectorXd solveSVD(const Eigen::MatrixXd& A, const Eigen::VectorXd& b);

#endif // SOLVESVD_HPP