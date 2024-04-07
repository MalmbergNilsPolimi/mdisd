#include "RBFinterpolator.hpp"
#include "RBFunctions.hpp"

Eigen::VectorXd RBFInterpolator::interpolate(const Eigen::MatrixXd& parametersFORinterp,
                                             const Eigen::MatrixXd& parameters,
                                             const Eigen::VectorXd& measurements) const {

    size_t num_measures{static_cast<size_t>(measurements.size())};
    size_t num_points{static_cast<size_t>(parametersFORinterp.rows())};

    Eigen::VectorXd results = Eigen::VectorXd::Zero(num_points);

    // Computation of the weights
    Eigen::MatrixXd coeff(num_measures, num_measures);

    for (size_t i = 0; i < num_measures; ++i)
    {
        for (size_t j = 0; j < num_measures; ++j)
        {
            coeff(i, j)=rbfunction((parameters.row(i)-parameters.row(j)).norm(), r0);
        }
    }
    
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(coeff, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd weights = svd.solve(measurements);
    
    // Computation of the interpolated value
    for (size_t k = 0; k < num_points; ++k)
    {
        for (size_t l = 0; l < num_measures; ++l) {
            results(k) += weights(l) * rbfunction((parametersFORinterp.row(k)-parameters.row(l)).norm(), r0);
        }
    }
    
    return results;
}