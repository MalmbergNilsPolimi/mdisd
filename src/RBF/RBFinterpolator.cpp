#include <iostream>
#include "RBFinterpolator.hpp"
#include "RBFunctions.hpp"

Eigen::VectorXd RBFInterpolator::interpolate(const Eigen::MatrixXd& parametersFORinterp,
                                             const Eigen::MatrixXd& parameters,
                                             const Eigen::VectorXd& measurements,
                                             Eigen::VectorXd* regression) const {
    
    size_t num_params{static_cast<size_t>(parameters.cols())};
    size_t num_measures{static_cast<size_t>(measurements.size())};
    size_t num_points{static_cast<size_t>(parametersFORinterp.rows())};

    Eigen::VectorXd results = Eigen::VectorXd::Zero(num_points);

    // Computation of the weights

    Eigen::MatrixXd coeff(num_measures + (polynomialRBF ? num_params + 1 : 0), num_measures + (polynomialRBF ? num_params + 1 : 0));
    coeff.setZero();

    for (size_t i = 0; i < num_measures; ++i)
    {
        for (size_t j = 0; j < num_measures; ++j)
        {
            coeff(i, j)=rbfunction((parameters.row(i)-parameters.row(j)).norm(), r0);            
        }

        if (polynomialRBF)
        {
            for (size_t k = num_measures; k < num_measures+num_params; ++k)
            {
                coeff(i,k)=parameters(i,k-num_measures);
                coeff(k,i)=(parameters.transpose())(k-num_measures,i);
            }
            
            coeff(i , num_measures+num_params) = 1;
            coeff(num_measures+num_params , i) = 1;            
        }


    }    

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(coeff, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::VectorXd weights(num_params  + (polynomialRBF ? num_params + 1 : 0));
    
    if (normalizeRBF)
    {
        if (polynomialRBF)
        {
            std::cerr << "Error: can't use normalized RBF method with polynomial addition simultaneously, as normalization only affects radial basis functions and not polynomials." << std::endl;
        }
        
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
        
    } else if (polynomialRBF)
    {
        Eigen::VectorXd NEWmeasurements = Eigen::VectorXd::Zero(num_measures+num_params+1);
        for (size_t i = 0; i < num_measures; ++i)
        {
            NEWmeasurements(i)=measurements(i);
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

        if (polynomialRBF)
        {
            for (size_t i = num_measures; i < num_measures+num_params; ++i)
            {
                results(k) += weights(i)*parametersFORinterp(k,i-num_measures);
            }
            results(k) += weights(num_measures+num_params);
        }
              
    }
    
    return results;
}

Eigen::VectorXd RBFInterpolator::interpolate(const Eigen::MatrixXd& parametersFORinterp,
                                              const Eigen::MatrixXd& parameters,
                                              const Eigen::VectorXd& measurements) const {
    return interpolate(parametersFORinterp, parameters, measurements, nullptr);
}