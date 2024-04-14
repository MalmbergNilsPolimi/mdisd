#include <iostream>
#include <cmath>
#include <random>

#include "RBFunctions.hpp"
#include "RBFinterpolator.hpp"
#include "OLSinterpolator.hpp"

template<typename T>
T generateRandomNumber(T min, T max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dis(min, max);
    return dis(gen);
}

int main() {

    ////////////////////////////////////////////////
    //////// TEST : Using random sampling //////////
    ////////////////////////////////////////////////

    // Function that return many sets of parameters.
    auto fillParameters = [](int num_sets, int num_params, double inf, double sup) {
        Eigen::MatrixXd parameters(num_sets, num_params);
        for (size_t i = 0; i < num_sets; ++i) {
            for (size_t j = 0; j < num_params; ++j) {
                
                double a{generateRandomNumber(inf , sup)};
                double b{generateRandomNumber(inf , sup)};

                if (a > b)
                {
                    double tmp = a;
                    a = b;
                    b = tmp;
                }
                
                parameters(i, j) = generateRandomNumber(a, b);
            }
        }
        return parameters;
    };

    // Definition of the function to interpolate (THE USER CAN CHANGE HERE THE FUNCTION)
    auto funToInterpolate = [](const Eigen::VectorXd& params) {
        double res = 0;
        for (int i = 0; i < params.size(); ++i) {
            res += std::pow(params(i), i);
        }
        return res;
    };


    // Dimensions of the problem
    size_t num_params{4};       // number of parameters taken by the function
    size_t num_measures{100};   // number of known points
    size_t num_points{3};       // number of points to interpolate
    
    // [inf, sup] is the interval of definition of the parameter
    double inf{0.};
    double sup{20.};

    // Declaration of the known parameters, known measurements and set of parameters used for the interpolation
    Eigen::MatrixXd parameters(num_measures, num_params);
    Eigen::VectorXd measurements(num_measures);
    Eigen::MatrixXd parametersFORinterp(num_points, num_params);

    // Creation of a random matrix containing the parameters for each measurement    
    parameters = fillParameters(num_measures, num_params, inf, sup);

    // Creation of the measurements 
    for (size_t i = 0; i < num_measures; i++)
    {
        measurements(i) = funToInterpolate(parameters.row(i));
    }

    // Creation of the parameters for which the user want the estimated output
    parametersFORinterp = fillParameters(num_points, num_params, inf, sup);


    // Use of RBF interpolation method
    double scale_factor{0.5};   // (THE USER CAN CHANGER HERE THE VALUE OF THE SCALE FACTOR)
    RBFInterpolator interpolatorRBF(&RBFunctions::multiquadratic, scale_factor); // (THE USER CAN CHANGE HERE THE USED RBFUNCTION)
    Eigen::VectorXd RBF_points_interpolated = interpolatorRBF.interpolate(parametersFORinterp, parameters, measurements);

    // Use of OLS approximation method
    OLSInterpolator interpolatorOLS;
    Eigen::VectorXd OLS_points_interpolated = interpolatorOLS.interpolate(parametersFORinterp, parameters, measurements);
    Eigen::VectorXd points_real(num_points);

    // Computing of the real values
    for (size_t i = 0; i < num_points; ++i)
    {
        points_real(i) = funToInterpolate(parametersFORinterp.row(i));
    }


    //////////////////////////////////////////////////
    ////////// PRINT OF RBF AND OLS RESULTS //////////
    //////////////////////////////////////////////////

    std::cout << "____________________USING THE REAL FUNCTION___________________" << std::endl;
    std::cout << "Real value:                 " << points_real.transpose() << std::endl;
    
    std::cout << std::endl;
    
    std::cout << "__________________USING RADIAL BASIS FUNCTIONS__________________" << std::endl;
    std::cout << "Interpolated value:         " << RBF_points_interpolated.transpose() << std::endl;
    std::cout << "Algebraic relative error:   " << ((RBF_points_interpolated - points_real).array() / points_real.array().abs()).transpose() << std::endl;

    std::cout << std::endl;

    std::cout << "__________________USING ORDINARY LEAST SQUARES__________________" << std::endl;
    std::cout << "Interpolated value:         " << OLS_points_interpolated.transpose() << std::endl;
    std::cout << "Algebraic relative error:   " << ((OLS_points_interpolated - points_real).array() / points_real.array().abs()).transpose() << std::endl;

    return 0;
}