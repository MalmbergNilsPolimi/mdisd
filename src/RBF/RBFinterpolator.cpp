#include "RBFinterpolator.hpp"
#include "RBFunctions.hpp"

Eigen::VectorXd RBFInterpolator::interpolate(const Eigen::MatrixXd& parametersFORinterp,
                                             const Eigen::MatrixXd& parameters,
                                             const Eigen::VectorXd& measurements,
                                             Eigen::VectorXd* regression) const {

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
    Eigen::VectorXd weights(parameters.cols());
    
    if (normalizeRBF)
    {
        Eigen::VectorXd NEWmeasurements = Eigen::VectorXd::Zero(num_measures);
        for (size_t i = 0; i < num_measures; ++i)
        {
            double sum1{};
            for (size_t j = 0; j < num_measures; ++j)
            {
                sum1 += rbfunction((parameters.row(i)-parameters.row(j)).norm(), r0);
            }  
            NEWmeasurements(i) = measurements(i) * sum1;
        }

        
        weights = svd.solve(NEWmeasurements);
        
    } else {
        weights = svd.solve(measurements);
    }

    // Storage of the weights if the user define the pointer to the VectorXd
    if (regression) {
        *regression = weights;
    }
    
    // Computation of the interpolated value
    Eigen::VectorXd normalize_part = Eigen::VectorXd::Ones(num_points);
    


    for (size_t k = 0; k < num_points; ++k)
    {
        if (normalizeRBF)
        {
            normalize_part(k) = 0;
            for (size_t l = 0; l < num_measures; ++l)
            {
                normalize_part(k) += rbfunction((parametersFORinterp.row(k)-parameters.row(l)).norm(), r0);
            }
            
        }
        
        for (size_t l = 0; l < num_measures; ++l) {           
            results(k) += weights(l) * rbfunction((parametersFORinterp.row(k)-parameters.row(l)).norm(), r0) / normalize_part(k);
        }
    }
    
    return results;
}

Eigen::VectorXd RBFInterpolator::interpolate(const Eigen::MatrixXd& parametersFORinterp,
                                              const Eigen::MatrixXd& parameters,
                                              const Eigen::VectorXd& measurements) const {
    return interpolate(parametersFORinterp, parameters, measurements, nullptr);
}