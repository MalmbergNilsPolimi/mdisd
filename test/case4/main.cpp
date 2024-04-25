#include <iostream>
#include <cmath>
#include <iomanip>
#include <random>
#include <fstream>
#include <filesystem>

#include "RBFunctions.hpp"
#include "RBFinterpolator.hpp"
#include "OLSinterpolator.hpp"
#include "rescaling.hpp"

template<typename T>
T generateRandomNumber(T min, T max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dis(min, max);
    return dis(gen);
}

void printTable(const Eigen::VectorXd& dimensions,
                const Eigen::VectorXd& real_values,
                const Eigen::VectorXd& RBF_interpolations_simple,
                const Eigen::VectorXd& RBF_errors_simple,
                const Eigen::VectorXd& RBF_interpolations_minmax,
                const Eigen::VectorXd& RBF_errors_minmax,
                const Eigen::VectorXd& RBF_interpolations_mean,
                const Eigen::VectorXd& RBF_errors_mean,
                const Eigen::VectorXd& RBF_interpolations_zscore,
                const Eigen::VectorXd& RBF_errors_zscore) {
    // Affichage du tableau
    std::cout << "____________________________________________________________________________________________________________________________________________________________________________________________" << std::endl;
    std::cout << std::setw(6) << std::left << "|| "
              << std::setw(20) << std::left << "|| Real value"
              << std::setw(20) << std::left << "|| no norm. interp."
              << std::setw(20) << std::left << "| no norm. error"
              << std::setw(20) << std::left << "|| min-max interp."
              << std::setw(20) << std::left << "| min-max error"
              << std::setw(20) << std::left << "|| mean interp."
              << std::setw(20) << std::left << "| mean error"
              << std::setw(20) << std::left << "|| z-score interp."
              << std::setw(20) << std::left << "| z-score error"
              << "||" << std::endl;
    std::cout << "____________________________________________________________________________________________________________________________________________________________________________________________" << std::endl;
    for (int i = 0; i < dimensions.size(); ++i) {
        std::cout << "||" << std::setw(4) << std::left << dimensions(i)
                  << "|| " << std::setw(17) << std::left << real_values(i)
                  << "|| " << std::setw(17) << std::left << RBF_interpolations_simple(i)
                  << "| " << std::setw(18) << std::left << RBF_errors_simple(i)
                  << "|| " << std::setw(17) << std::left << RBF_interpolations_minmax(i)
                  << "| " << std::setw(18) << std::left << RBF_errors_minmax(i)
                  << "|| " << std::setw(17) << std::left << RBF_interpolations_mean(i)
                  << "| " << std::setw(18) << std::left << RBF_errors_mean(i)
                  << "|| " << std::setw(17) << std::left << RBF_interpolations_zscore(i)
                  << "| " << std::setw(18) << std::left << RBF_errors_zscore(i)
                  << "|| " << std::endl;
    }
    std::cout << "____________________________________________________________________________________________________________________________________________________________________________________________" << std::endl;
}


void plotData(const Eigen::VectorXd& error_RBF_BR, const Eigen::VectorXd& error_RBF_AR1,
              const Eigen::VectorXd& error_RBF_AR2, const Eigen::VectorXd& error_RBF_AR3,
              const std::string& title, const std::string& RBF, 
              const bool EXPORT) {

    std::filesystem::create_directories("./plot/");
    std::filesystem::create_directories("./plot/files/");
    std::filesystem::create_directories("./plot/figures/");

    // RBF error BR.
    std::ofstream dataFileRBF_BR("./plot/files/error_RBF_"+RBF+"_BR.dat");
    if (!dataFileRBF_BR.is_open()) {
        std::cerr << "Error: Unable to open data file." << std::endl;
        return;
    }

    for (int i = 0; i < error_RBF_BR.size(); ++i) {
        dataFileRBF_BR << i+1 << " " << error_RBF_BR(i) << std::endl;
    }
    dataFileRBF_BR.close();

    // RBF error AR1.
    std::ofstream dataFileRBF_AR1("./plot/files/error_RBF_"+RBF+"_AR1.dat");
    if (!dataFileRBF_AR1.is_open()) {
        std::cerr << "Error: Unable to open data file." << std::endl;
        return;
    }

    for (int i = 0; i < error_RBF_AR1.size(); ++i) {
        dataFileRBF_AR1 << i+1 << " " << error_RBF_AR1(i) << std::endl;
    }
    dataFileRBF_AR1.close();

    // RBF error AR2.
    std::ofstream dataFileRBF_AR2("./plot/files/error_RBF_"+RBF+"_AR2.dat");
    if (!dataFileRBF_AR2.is_open()) {
        std::cerr << "Error: Unable to open data file." << std::endl;
        return;
    }

    for (int i = 0; i < error_RBF_AR2.size(); ++i) {
        dataFileRBF_AR2 << i+1 << " " << error_RBF_AR2(i) << std::endl;
    }
    dataFileRBF_AR2.close();

    // RBF error AR3.
    std::ofstream dataFileRBF_AR3("./plot/files/error_RBF_"+RBF+"_AR3.dat");
    if (!dataFileRBF_AR3.is_open()) {
        std::cerr << "Error: Unable to open data file." << std::endl;
        return;
    }

    for (int i = 0; i < error_RBF_AR3.size(); ++i) {
        dataFileRBF_AR3 << i+1 << " " << error_RBF_AR3(i) << std::endl;
    }
    dataFileRBF_AR3.close();

    //std::string title{"Relative error on a 4D interpolation"};

    // GNUplot commands
    std::ofstream gnuplotScript("./plot/files/plot_script_"+RBF+".gnu");
    
    if (EXPORT)
    {
        // To export the plot in svg format
        gnuplotScript << "set terminal svg" << std::endl;
        gnuplotScript << "set output \"" << "./plot/figures/case4_plot_"+RBF+".svg" << "\"" << std::endl;
    }

    gnuplotScript << "set title \"" << title << "\"" << std::endl;
    gnuplotScript << "set xlabel \"interpolated point number\"" << std::endl;
    gnuplotScript << "set ylabel \"relative error\"" << std::endl;
    gnuplotScript << "set key box" << std::endl;
    gnuplotScript << "set style line 1 lw 2 lc rgb '#990042'" << std::endl;
    gnuplotScript << "set style line 2 lw 2 lc rgb '#31f120'" << std::endl;
    gnuplotScript << "set style line 3 lw 2 lc rgb '#0044a5'" << std::endl;
    gnuplotScript << "set style line 4 lw 2 lc rgb '#888888'" << std::endl;
    gnuplotScript << "plot \"" << "./plot/files/error_RBF_"+RBF+"_BR.dat" << "\" with lines linestyle 1 title \"without rescaling\","
                << " \"" << "./plot/files/error_RBF_"+RBF+"_AR1.dat" << "\" with points pointtype 7 title \" Min-Max rescaling\","
                << " \"" << "./plot/files/error_RBF_"+RBF+"_AR2.dat" << "\" with lines linestyle 3 title \"Mean rescaling\","
                << " \"" << "./plot/files/error_RBF_"+RBF+"_AR3.dat" << "\" with lines linestyle 4 title \"Z-score rescaling\""
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
    std::string path{"gnuplot -persist ./plot/files/plot_script_"+RBF+".gnu"};
    system(path.c_str());
}




int main() {

    ////////////////////////////////////////////////
    /////////// SAMPLING OF THE POINTS /////////////
    ////////////////////////////////////////////////

    // Function that return many sets of parameters.
    auto fillParameters = [](int num_sets, int num_params, Eigen::VectorXd tab_inf, Eigen::VectorXd tab_sup) {
        Eigen::MatrixXd parameters(num_sets, num_params);
        for (size_t i = 0; i < num_sets; ++i) {
            for (size_t j = 0; j < num_params; ++j) {
                parameters(i, j) = generateRandomNumber(tab_inf(j), tab_sup(j));
            }
        }
        return parameters;
    };

    // Definition of the function to interpolate (THE USER CAN CHANGE HERE THE FUNCTION)
    auto funToInterpolate = [](const Eigen::VectorXd& params) {
        double res{};

        if (params.size()==4)
        {
            res = (3./2)*params(0)*params(0)*cos(params(1)*M_PI)*sin(params(3)-params(2)) - params(3)*exp(-(params(0)+params(2)/2.)) + log(5*abs(params(2))) - 18.12*std::pow(abs(params(1)), params(3));
        } else {
            for (int i = 0; i < params.size(); ++i) {
                res += params(i)*exp(params(i)/2.);
            }
        }
        return res;
    };


    // Dimensions of the problem
    size_t num_params{4};       // number of parameters taken by the function
    size_t num_measures{150};   // number of known points
    size_t num_points{10};      // number of points to interpolate
    
    // [inf, sup] is the interval of definition of the parameter
    Eigen::VectorXd tab_inf(num_params);
    tab_inf << -1. , 10. , 1. , 0.; 
    Eigen::VectorXd tab_sup(num_params);
    tab_inf << 1. , 20. , 5. , 0.5; 

    // Declaration of the known parameters, known measurements and set of parameters used for the interpolation
    Eigen::MatrixXd parameters(num_measures, num_params);
    Eigen::VectorXd measurements(num_measures);
    Eigen::MatrixXd parametersFORinterp(num_points, num_params);

    // Creation of a random matrix containing the parameters for each measurement    
    parameters = fillParameters(num_measures, num_params, tab_inf, tab_sup);

    // Creation of the measurements 
    for (size_t i = 0; i < num_measures; i++)
    {
        measurements(i) = funToInterpolate(parameters.row(i));
    }

    // Creation of the parameters for which the user want the estimated output
    parametersFORinterp = fillParameters(num_points, num_params, tab_inf, tab_sup);

    Eigen::VectorXd points_real(num_points);
    // Computing of the real values
    for (size_t i = 0; i < num_points; ++i)
    {
        points_real(i) = funToInterpolate(parametersFORinterp.row(i));
    }

    ////////////////////////////////////////////////
    //////////// BR : BEFORE RESCALING /////////////
    ////////////////////////////////////////////////

    // Use of RBF interpolation method
    double scale_factor{0.5};
    RBFInterpolator interpolatorRBF(&RBFunctions::gaussian, scale_factor);
    Eigen::VectorXd RBF_points_interpolated_BR = interpolatorRBF.interpolate(parametersFORinterp, parameters, measurements);

    // Use of normalized RBF interpolation method
    bool normalize{true};
    RBFInterpolator interpolatorRBFnorm(&RBFunctions::gaussian, scale_factor, normalize);
    Eigen::VectorXd RBF_norm_points_interpolated_BR = interpolatorRBFnorm.interpolate(parametersFORinterp, parameters, measurements);

    // Use of RBF augmented with polynomial
    normalize=false;
    bool polynomial{true};
    RBFInterpolator interpolatorRBFpoly(&RBFunctions::gaussian, scale_factor, normalize, polynomial);
    Eigen::VectorXd RBF_poly_points_interpolated_BR = interpolatorRBFpoly.interpolate(parametersFORinterp, parameters, measurements);

    // Use of OLS approximation method
    OLSInterpolator interpolatorOLS;
    Eigen::VectorXd OLS_points_interpolated_BR = interpolatorOLS.interpolate(parametersFORinterp, parameters, measurements);


    Eigen::VectorXd RBF_relative_errors_BR = ((RBF_points_interpolated_BR - points_real).array() / points_real.array().abs());
    Eigen::VectorXd RBF_norm_relative_errors_BR = ((RBF_norm_points_interpolated_BR - points_real).array() / points_real.array().abs());
    Eigen::VectorXd RBF_poly_relative_errors_BR = ((RBF_poly_points_interpolated_BR - points_real).array() / points_real.array().abs());
    Eigen::VectorXd OLS_relative_errors_BR = ((OLS_points_interpolated_BR - points_real).array() / points_real.array().abs());



    ////////////////////////////////////////////////
    //////// AR1 : AFTER MIN-MAX RESCALING /////////
    ////////////////////////////////////////////////

    Rescaling rescaling_minmax;
    auto normalizedData_minmax = rescaling_minmax.minMaxNormalization(parameters, &parametersFORinterp);

    Eigen::VectorXd RBF_points_interpolated_AR1 = interpolatorRBF.interpolate(normalizedData_minmax.second, normalizedData_minmax.first, measurements);
    Eigen::VectorXd RBF_norm_points_interpolated_AR1 = interpolatorRBFnorm.interpolate(normalizedData_minmax.second, normalizedData_minmax.first, measurements);
    Eigen::VectorXd RBF_poly_points_interpolated_AR1 = interpolatorRBFpoly.interpolate(normalizedData_minmax.second, normalizedData_minmax.first, measurements);
    Eigen::VectorXd OLS_points_interpolated_AR1 = interpolatorOLS.interpolate(normalizedData_minmax.second, normalizedData_minmax.first, measurements);


    Eigen::VectorXd RBF_relative_errors_AR1 = ((RBF_points_interpolated_AR1 - points_real).array() / points_real.array().abs());
    Eigen::VectorXd RBF_norm_relative_errors_AR1 = ((RBF_norm_points_interpolated_AR1 - points_real).array() / points_real.array().abs());
    Eigen::VectorXd RBF_poly_relative_errors_AR1 = ((RBF_poly_points_interpolated_AR1 - points_real).array() / points_real.array().abs());
    Eigen::VectorXd OLS_relative_errors_AR1 = ((OLS_points_interpolated_AR1 - points_real).array() / points_real.array().abs());

    ////////////////////////////////////////////////
    ////////// AR2 : AFTER MEAN RESCALING //////////
    ////////////////////////////////////////////////

    //  Use of RBF interpolation method with mean normalization
    Rescaling rescaling_mean;
    auto normalizedData_mean = rescaling_mean.meanNormalization(parameters, &parametersFORinterp);

    Eigen::VectorXd RBF_points_interpolated_AR2 = interpolatorRBF.interpolate(normalizedData_mean.second, normalizedData_mean.first, measurements);
    Eigen::VectorXd RBF_norm_points_interpolated_AR2 = interpolatorRBFnorm.interpolate(normalizedData_mean.second, normalizedData_mean.first, measurements);
    Eigen::VectorXd RBF_poly_points_interpolated_AR2 = interpolatorRBFpoly.interpolate(normalizedData_mean.second, normalizedData_mean.first, measurements);
    Eigen::VectorXd OLS_points_interpolated_AR2 = interpolatorOLS.interpolate(normalizedData_mean.second, normalizedData_mean.first, measurements);


    Eigen::VectorXd RBF_relative_errors_AR2 = ((RBF_points_interpolated_AR2 - points_real).array() / points_real.array().abs());
    Eigen::VectorXd RBF_norm_relative_errors_AR2 = ((RBF_norm_points_interpolated_AR2 - points_real).array() / points_real.array().abs());
    Eigen::VectorXd RBF_poly_relative_errors_AR2 = ((RBF_poly_points_interpolated_AR2 - points_real).array() / points_real.array().abs());
    Eigen::VectorXd OLS_relative_errors_AR2 = ((OLS_points_interpolated_AR2 - points_real).array() / points_real.array().abs());

    ////////////////////////////////////////////////
    ///////// AR3 : AFTER Z-SCORE RESCALING ////////
    ////////////////////////////////////////////////

    //  Use of RBF interpolation method with mean normalization
    Rescaling rescaling_zscore;
    auto normalizedData_zscore = rescaling_zscore.zScoreNormalization(parameters, &parametersFORinterp);

    Eigen::VectorXd RBF_points_interpolated_AR3 = interpolatorRBF.interpolate(normalizedData_zscore.second, normalizedData_zscore.first, measurements);
    Eigen::VectorXd RBF_norm_points_interpolated_AR3 = interpolatorRBFnorm.interpolate(normalizedData_zscore.second, normalizedData_zscore.first, measurements);
    Eigen::VectorXd RBF_poly_points_interpolated_AR3 = interpolatorRBFpoly.interpolate(normalizedData_zscore.second, normalizedData_zscore.first, measurements);
    Eigen::VectorXd OLS_points_interpolated_AR3 = interpolatorOLS.interpolate(normalizedData_zscore.second, normalizedData_zscore.first, measurements);


    Eigen::VectorXd RBF_relative_errors_AR3 = ((RBF_points_interpolated_AR3 - points_real).array() / points_real.array().abs());
    Eigen::VectorXd RBF_norm_relative_errors_AR3 = ((RBF_norm_points_interpolated_AR3 - points_real).array() / points_real.array().abs());
    Eigen::VectorXd RBF_poly_relative_errors_AR3 = ((RBF_poly_points_interpolated_AR3 - points_real).array() / points_real.array().abs());
    Eigen::VectorXd OLS_relative_errors_AR3 = ((OLS_points_interpolated_AR3 - points_real).array() / points_real.array().abs());




    ////////////////////////////////////////////////
    ///////// PRINT AND PLOT OF THE RESULTS ////////
    ////////////////////////////////////////////////

    Eigen::VectorXd points = Eigen::VectorXd::LinSpaced(num_points, 1, num_points);
    
    bool PLOT{true};
    bool EXPORT{true};

    std::cout << "\n-----------------------------------------------------------------------------------------SIMPLE RBF-----------------------------------------------------------------------------------------" << std::endl;

    printTable(points, points_real, RBF_points_interpolated_BR, RBF_relative_errors_BR, RBF_points_interpolated_AR1, RBF_relative_errors_AR1, RBF_points_interpolated_AR2, RBF_relative_errors_AR2, RBF_points_interpolated_AR3, RBF_relative_errors_AR3);
    std::string title{"Relative error on a 4D interpolation with the RBF method"};
    std::string RBF{"simple"};
    
    if (PLOT)
    {
        plotData(RBF_relative_errors_BR, RBF_relative_errors_AR1, RBF_relative_errors_AR2, RBF_relative_errors_AR3, title, RBF, EXPORT);
    }

    std::cout << "\n---------------------------------------------------------------------------------------NORMALIZED RBF---------------------------------------------------------------------------------------" << std::endl;

    printTable(points, points_real, RBF_norm_points_interpolated_BR, RBF_norm_relative_errors_BR, RBF_norm_points_interpolated_AR1, RBF_norm_relative_errors_AR1, RBF_norm_points_interpolated_AR2, RBF_norm_relative_errors_AR2, RBF_norm_points_interpolated_AR3, RBF_norm_relative_errors_AR3);
    title = "Relative error on a 4D interpolation with the NRBF method";
    RBF = "norm";
    
    if (PLOT)
    {
        plotData(RBF_norm_relative_errors_BR, RBF_norm_relative_errors_AR1, RBF_norm_relative_errors_AR2, RBF_norm_relative_errors_AR3, title, RBF, EXPORT);
    }



    std::cout << "\n-------------------------------------------------------------------------------------RBF WITH POLYNOMIAL------------------------------------------------------------------------------------" << std::endl;

    printTable(points, points_real, RBF_poly_points_interpolated_BR, RBF_poly_relative_errors_BR, RBF_poly_points_interpolated_AR1, RBF_poly_relative_errors_AR1, RBF_poly_points_interpolated_AR2, RBF_poly_relative_errors_AR2, RBF_poly_points_interpolated_AR3, RBF_poly_relative_errors_AR3);
    title = "Relative error on a 4D interpolation with the RBFP method";
    RBF = "poly";

    if (PLOT)
    {
        plotData(RBF_poly_relative_errors_BR, RBF_poly_relative_errors_AR1, RBF_poly_relative_errors_AR2, RBF_poly_relative_errors_AR3, title, RBF, EXPORT);
    }

    return 0;
}