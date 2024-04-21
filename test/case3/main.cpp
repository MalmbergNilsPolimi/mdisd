#include <iostream>
#include <cmath>
#include <iomanip>
#include <random>
#include <fstream>
#include <filesystem>

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

void printTable(const Eigen::VectorXd& dimensions,
                const Eigen::VectorXd& real_values,
                const Eigen::VectorXd& RBF_interpolations,
                const Eigen::VectorXd& RBF_errors,
                const Eigen::VectorXd& RBF_norm_interpolations,
                const Eigen::VectorXd& RBF_norm_errors,
                const Eigen::VectorXd& RBF_poly_interpolations,
                const Eigen::VectorXd& RBF_poly_errors,
                const Eigen::VectorXd& OLS_interpolations,
                const Eigen::VectorXd& OLS_errors) {
    // Affichage du tableau
    std::cout << "____________________________________________________________________________________________________________________________________________" << std::endl;
    std::cout << std::setw(4) << std::left << "|| "
              << std::setw(14) << std::left << "|| Real value"
              << std::setw(16) << std::left << "|| RBF interp."
              << std::setw(14) << std::left << "| RBF error"
              << std::setw(16) << std::left << "|| NRBF interp."
              << std::setw(14) << std::left << "| NRBF error"
              << std::setw(16) << std::left << "|| RBFP interp."
              << std::setw(14) << std::left << "| RBFP error"
              << std::setw(16) << std::left << "|| OLS interp."
              << std::setw(14) << std::left << "| OLS error"
              << "||" << std::endl;
    std::cout << "____________________________________________________________________________________________________________________________________________" << std::endl;
    for (int i = 0; i < dimensions.size(); ++i) {
        std::cout << "||" << std::setw(2) << std::left << dimensions(i)
                  << "|| " << std::setw(11) << std::left << real_values(i)
                  << "|| " << std::setw(13) << std::left << RBF_interpolations(i)
                  << "| " << std::setw(12) << std::left << RBF_errors(i)
                  << "|| " << std::setw(13) << std::left << RBF_norm_interpolations(i)
                  << "| " << std::setw(12) << std::left << RBF_norm_errors(i)
                  << "|| " << std::setw(13) << std::left << RBF_poly_interpolations(i)
                  << "| " << std::setw(12) << std::left << RBF_poly_errors(i)
                  << "|| " << std::setw(13) << std::left << OLS_interpolations(i)
                  << "| " << std::setw(12) << std::left << OLS_errors(i)
                  << "|| " << std::endl;
    }
    std::cout << "____________________________________________________________________________________________________________________________________________" << std::endl;
}


void plotData(const Eigen::VectorXd& error_RBF, const Eigen::VectorXd& error_NRBF,
              const Eigen::VectorXd& error_RBFP, const Eigen::VectorXd& error_OLS, 
              const bool EXPORT) {

    std::filesystem::create_directories("./plot/");
    std::filesystem::create_directories("./plot/files/");
    std::filesystem::create_directories("./plot/figures/");

    // RBF error
    std::ofstream dataFileRBF("./plot/files/error_RBF.dat");
    if (!dataFileRBF.is_open()) {
        std::cerr << "Error: Unable to open data file." << std::endl;
        return;
    }

    for (int i = 0; i < error_RBF.size(); ++i) {
        dataFileRBF << i+1 << " " << error_RBF(i) << std::endl;
    }
    dataFileRBF.close();

    // NRBF error
    std::ofstream dataFileRBFnorm("./plot/files/error_NRBF.dat");
    if (!dataFileRBFnorm.is_open()) {
        std::cerr << "Error: Unable to open data file." << std::endl;
        return;
    }

    for (int i = 0; i < error_NRBF.size(); ++i) {
        dataFileRBFnorm << i+1 << " " << error_NRBF(i) << std::endl;
    }
    dataFileRBFnorm.close();

    // RBFP error
    std::ofstream dataFileRBFpoly("./plot/files/error_RBFP.dat");
    if (!dataFileRBFpoly.is_open()) {
        std::cerr << "Error: Unable to open data file." << std::endl;
        return;
    }

    for (int i = 0; i < error_RBFP.size(); ++i) {
        dataFileRBFpoly << i+1 << " " << error_RBFP(i) << std::endl;
    }
    dataFileRBFpoly.close();

    // OLS error
    std::ofstream dataFileOLS("./plot/files/error_OLS.dat");
    if (!dataFileOLS.is_open()) {
        std::cerr << "Error: Unable to open data file." << std::endl;
        return;
    }

    for (int i = 0; i < error_OLS.size(); ++i) {
        dataFileOLS << i+1 << " " << error_OLS(i) << std::endl;
    }
    dataFileOLS.close();

    std::string title{"Relative error on a 4D interpolation"};

    // GNUplot commands
    std::ofstream gnuplotScript("./plot/files/plot_script.gnu");
    
    if (EXPORT)
    {
        // To export the plot in svg format
        gnuplotScript << "set terminal svg" << std::endl;
        gnuplotScript << "set output \"" << "./plot/figures/case3_plot.svg" << "\"" << std::endl;
    }

    gnuplotScript << "set title \"" << title << "\"" << std::endl;
    gnuplotScript << "set xlabel \"interpolated point number\"" << std::endl;
    gnuplotScript << "set ylabel \"relative error\"" << std::endl;
    gnuplotScript << "set key box" << std::endl;
    gnuplotScript << "set style line 1 lw 2 lc rgb '#990042'" << std::endl;
    gnuplotScript << "set style line 2 lw 2 lc rgb '#31f120'" << std::endl;
    gnuplotScript << "set style line 3 lw 2 lc rgb '#0044a5'" << std::endl;
    gnuplotScript << "set style line 4 lw 2 lc rgb '#888888'" << std::endl;
    gnuplotScript << "plot \"" << "./plot/files/error_RBF.dat" << "\" with lines linestyle 1 title \"RBF error\","
                << " \"" << "./plot/files/error_NRBF.dat" << "\" with lines linestyle 2 title \"NRBF error\","
                << " \"" << "./plot/files/error_RBFP.dat" << "\" with lines linestyle 3 title \"RBFP error\","
                << " \"" << "./plot/files/error_OLS.dat" << "\" with lines linestyle 4 title \"OLS error\""
                << std::endl;
    gnuplotScript << "set key top right" << std::endl;
        
    
    if (EXPORT)
    {   
        // The export "cancel" the plot
        // Replot to see it in a window
        gnuplotScript << "set terminal wxt" << std::endl;
        gnuplotScript << "replot" << std::endl;
    }

    gnuplotScript.close();

    // Execute GNUplot
    system("gnuplot -persist ./plot/files/plot_script.gnu");
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
    size_t num_points{10};       // number of points to interpolate
    
    // [inf, sup] is the interval of definition of the parameter
    double inf{0.};
    double sup{1.};

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

    // Use of normalized RBF interpolation method
    bool normalize{true};
    RBFInterpolator interpolatorRBFnorm(&RBFunctions::multiquadratic, scale_factor, normalize);
    Eigen::VectorXd RBF_norm_points_interpolated = interpolatorRBFnorm.interpolate(parametersFORinterp, parameters, measurements);

    // Use of RBF augmented with polynomial
    normalize=false;
    bool polynomial{true};
    RBFInterpolator interpolatorRBFpoly(&RBFunctions::multiquadratic, scale_factor, normalize, polynomial);
    Eigen::VectorXd RBF_poly_points_interpolated = interpolatorRBFpoly.interpolate(parametersFORinterp, parameters, measurements);

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

    Eigen::VectorXd RBF_relative_errors = ((RBF_points_interpolated - points_real).array() / points_real.array().abs());
    Eigen::VectorXd RBF_norm_relative_errors = ((RBF_norm_points_interpolated - points_real).array() / points_real.array().abs());
    Eigen::VectorXd RBF_poly_relative_errors = ((RBF_poly_points_interpolated - points_real).array() / points_real.array().abs());
    Eigen::VectorXd OLS_relative_errors = ((OLS_points_interpolated - points_real).array() / points_real.array().abs());

    Eigen::VectorXd points = Eigen::VectorXd::LinSpaced(num_points, 1, num_points);
    printTable(points, points_real, RBF_points_interpolated, RBF_relative_errors, RBF_norm_points_interpolated, RBF_norm_relative_errors, RBF_poly_points_interpolated, RBF_poly_relative_errors, OLS_points_interpolated, OLS_relative_errors);
    
    bool EXPORT{true};
    plotData(RBF_relative_errors, RBF_norm_relative_errors, RBF_poly_relative_errors, OLS_relative_errors, EXPORT);
    return 0;
}