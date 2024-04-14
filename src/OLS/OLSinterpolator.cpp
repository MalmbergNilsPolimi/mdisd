#include "OLSinterpolator.hpp"

Eigen::VectorXd OLSInterpolator::interpolate(const Eigen::MatrixXd& parametersFORinterp,
                                             const Eigen::MatrixXd& parameters,
                                             const Eigen::VectorXd& measurements,
                                             Eigen::VectorXd* regression) const {

    size_t num_params{static_cast<size_t>(parameters.cols())};
    size_t num_measures{static_cast<size_t>(measurements.size())};
    size_t num_points{static_cast<size_t>(parametersFORinterp.rows())};
    
    Eigen::VectorXd results = Eigen::VectorXd::Zero(num_points);

    Eigen::MatrixXd transpose_parameters{parameters.transpose()};
    Eigen::VectorXd beta_coeff(num_params);

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(transpose_parameters*parameters, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd inverse_matrix = svd.matrixV() * svd.singularValues().asDiagonal().inverse() * svd.matrixU().transpose();

    beta_coeff = inverse_matrix * transpose_parameters * measurements;

    // Storage of the weights if the user define the pointer to the VectorXd
    if (regression) {
        *regression = beta_coeff;
    }

    // Computation of the interpolated value
    for (size_t i = 0; i < num_points; ++i)
    {
        for (size_t j = 0; j < num_params; ++j) {
            results(i) += beta_coeff(j) * parametersFORinterp(i,j);
        }
    }
    return results;
}

Eigen::VectorXd OLSInterpolator::interpolate(const Eigen::MatrixXd& parametersFORinterp,
                                              const Eigen::MatrixXd& parameters,
                                              const Eigen::VectorXd& measurements) const {
    return interpolate(parametersFORinterp, parameters, measurements, nullptr);
}