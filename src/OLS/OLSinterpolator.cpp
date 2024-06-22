#include "OLSinterpolator.hpp"

Eigen::VectorXd OLSInterpolator::interpolate(const Eigen::MatrixXd& parametersFORinterp,
                                             const Eigen::MatrixXd& parameters,
                                             const Eigen::VectorXd& measurements,
                                             Eigen::VectorXd* regression) const {
    // Interpolation using OLS method.

    size_t num_params{static_cast<size_t>(parameters.cols())}; // number of variables of the function to interpolate.
    size_t num_points{static_cast<size_t>(parametersFORinterp.rows())}; // number of points to interpolate.
    
    // Add a column with ones in the first position.
    Eigen::MatrixXd newparameters(parameters.rows(), parameters.cols() + 1);
    newparameters << Eigen::MatrixXd::Ones(parameters.rows(), 1), parameters;

    Eigen::MatrixXd transpose_parameters{newparameters.transpose()};
    
    // Compute each matrix of the SVD decomposition.
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(transpose_parameters*newparameters, Eigen::ComputeThinU | Eigen::ComputeThinV);
    
    // Compute the weights of the OLS interpolated function.
    Eigen::VectorXd beta_coeff(num_params);
    Eigen::MatrixXd inverse_matrix = svd.matrixV() * svd.singularValues().asDiagonal().inverse() * svd.matrixU().transpose();
    beta_coeff = inverse_matrix * transpose_parameters * measurements;

    // Storage of the weights if the user define the pointer to the VectorXd.
    if (regression) {
        *regression = beta_coeff;
    }

    // Computation of the interpolated value.
    Eigen::VectorXd results{Eigen::VectorXd::Constant(num_points, beta_coeff(0))} ;

    for (size_t i = 0; i < num_points; ++i)
    {
        for (size_t j = 0; j < num_params; ++j)
        {
            results(i) += beta_coeff(j+1) * parametersFORinterp(i,j);
        }
    }
    return results;
}

Eigen::VectorXd OLSInterpolator::interpolate(const Eigen::MatrixXd& parametersFORinterp,
                                              const Eigen::MatrixXd& parameters,
                                              const Eigen::VectorXd& measurements) const {
    // Override of the function in case of no definition of a vector to store the weights.
    return interpolate(parametersFORinterp, parameters, measurements, nullptr);
}