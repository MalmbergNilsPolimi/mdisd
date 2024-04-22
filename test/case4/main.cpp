#include <iostream>
#include <cmath>
#include <iomanip>
#include <random>
#include <fstream>
#include <filesystem>

#include "RBFunctions.hpp"
#include "RBFinterpolator.hpp"
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
                const Eigen::VectorXd& RBF_interpolations,
                const Eigen::VectorXd& RBF_errors,
                const Eigen::VectorXd& RBF_mean_interpolations,
                const Eigen::VectorXd& RBF_mean_errors,
                const Eigen::VectorXd& RBF_minmax_interpolations,
                const Eigen::VectorXd& RBF_minmax_errors,
                const Eigen::VectorXd& RBF_zscore_interpolations,
                const Eigen::VectorXd& RBF_zscore_errors) {
    // Affichage du tableau
    std::cout << "____________________________________________________________________________________________________________________________________________" << std::endl;
    std::cout << std::setw(4) << std::left << "|| "
              << std::setw(14) << std::left << "|| real value"
              << std::setw(16) << std::left << "|| simple"
              << std::setw(14) << std::left << "| error"
              << std::setw(16) << std::left << "|| mean norm."
              << std::setw(14) << std::left << "| error"
              << std::setw(16) << std::left << "|| min-max norm."
              << std::setw(14) << std::left << "| error"
              << std::setw(16) << std::left << "|| z-score norm."
              << std::setw(14) << std::left << "| error"
              << "||" << std::endl;
    std::cout << "____________________________________________________________________________________________________________________________________________" << std::endl;
    for (int i = 0; i < dimensions.size(); ++i) {
        std::cout << "||" << std::setw(2) << std::left << dimensions(i)
                  << "|| " << std::setw(11) << std::left << real_values(i)
                  << "|| " << std::setw(13) << std::left << RBF_interpolations(i)
                  << "| " << std::setw(12) << std::left << RBF_errors(i)
                  << "|| " << std::setw(13) << std::left << RBF_mean_interpolations(i)
                  << "| " << std::setw(12) << std::left << RBF_mean_errors(i)
                  << "|| " << std::setw(13) << std::left << RBF_minmax_interpolations(i)
                  << "| " << std::setw(12) << std::left << RBF_minmax_errors(i)
                  << "|| " << std::setw(13) << std::left << RBF_zscore_interpolations(i)
                  << "| " << std::setw(12) << std::left << RBF_zscore_errors(i)
                  << "|| " << std::endl;
    }
    std::cout << "____________________________________________________________________________________________________________________________________________" << std::endl;
}


void plotData(const Eigen::VectorXd& error_RBF, const Eigen::VectorXd& error_RBF_mean,
              const Eigen::VectorXd& error_RBF_minmax, const Eigen::VectorXd& error_RBF_zscore, 
              const bool EXPORT) {

    std::filesystem::create_directories("./plot/");
    std::filesystem::create_directories("./plot/files/");
    std::filesystem::create_directories("./plot/figures/");

    // simple RBF error
    std::ofstream dataFileRBF("./plot/files/error_RBF.dat");
    if (!dataFileRBF.is_open()) {
        std::cerr << "Error: Unable to open data file." << std::endl;
        return;
    }

    for (int i = 0; i < error_RBF.size(); ++i) {
        dataFileRBF << i+1 << " " << error_RBF(i) << std::endl;
    }
    dataFileRBF.close();

    // mean normalized RBF error
    std::ofstream dataFileRBFmean("./plot/files/error_RBF_mean.dat");
    if (!dataFileRBFmean.is_open()) {
        std::cerr << "Error: Unable to open data file." << std::endl;
        return;
    }

    for (int i = 0; i < error_RBF_mean.size(); ++i) {
        dataFileRBFmean << i+1 << " " << error_RBF_mean(i) << std::endl;
    }
    dataFileRBFmean.close();

    // min-max normalized RBF error
    std::ofstream dataFileRBFminmax("./plot/files/error_RBF_minmax.dat");
    if (!dataFileRBFminmax.is_open()) {
        std::cerr << "Error: Unable to open data file." << std::endl;
        return;
    }

    for (int i = 0; i < error_RBF_minmax.size(); ++i) {
        dataFileRBFminmax << i+1 << " " << error_RBF_minmax(i) << std::endl;
    }
    dataFileRBFminmax.close();

    // z-score normalized RBF error
    std::ofstream dataFileRBFzscore("./plot/files/error_RBF_zscore.dat");
    if (!dataFileRBFzscore.is_open()) {
        std::cerr << "Error: Unable to open data file." << std::endl;
        return;
    }

    for (int i = 0; i < error_RBF_zscore.size(); ++i) {
        dataFileRBFzscore << i+1 << " " << error_RBF_zscore(i) << std::endl;
    }
    dataFileRBFzscore.close();

    std::string title{"Relative error on a 4D interpolation"};

    // GNUplot commands
    std::ofstream gnuplotScript("./plot/files/plot_script.gnu");
    
    if (EXPORT)
    {
        // To export the plot in svg format
        gnuplotScript << "set terminal svg" << std::endl;
        gnuplotScript << "set output \"" << "./plot/figures/case4_plot.svg" << "\"" << std::endl;
    }

    gnuplotScript << "set title \"" << title << "\"" << std::endl;
    gnuplotScript << "set xlabel \"interpolated point number\"" << std::endl;
    gnuplotScript << "set ylabel \"relative error\"" << std::endl;
    gnuplotScript << "set key box bottom right" << std::endl;
    gnuplotScript << "set style line 1 lw 2 lc rgb '#990042'" << std::endl;
    gnuplotScript << "set style line 2 lw 2 lc rgb '#31f120'" << std::endl;
    gnuplotScript << "set style line 3 lw 2 lc rgb '#0044a5'" << std::endl;
    gnuplotScript << "set style line 4 lw 2 lc rgb '#888888'" << std::endl;
    gnuplotScript << "plot \"" << "./plot/files/error_RBF.dat" << "\" with lines linestyle 1 title \"simple RBF error\","
                << " \"" << "./plot/files/error_RBF_mean.dat" << "\" with points pointtype 7 title \"mean norm. RBF error\","
                << " \"" << "./plot/files/error_RBF_minmax.dat" << "\" with lines linestyle 3 title \" min-max norm. RBF error\","
                << " \"" << "./plot/files/error_RBF_zscore.dat" << "\" with lines linestyle 4 title \"z-score norm. RBF error\""
                << std::endl;
        
    
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
            res += params(i)*exp(params(i)/2.);
        }
        return res;
    };


    // Dimensions of the problem
    size_t num_params{10};       // number of parameters taken by the function
    size_t num_measures{150};   // number of known points
    size_t num_points{30};       // number of points to interpolate
    
    // [inf, sup] is the interval of definition of the parameter
    double inf{-5.};
    double sup{5.};

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


    // Use of RBF interpolation method without rescaling
    double scale_factor{3};   // (THE USER CAN CHANGER HERE THE VALUE OF THE SCALE FACTOR)
    RBFInterpolator interpolatorRBF(&RBFunctions::gaussian, scale_factor); // (THE USER CAN CHANGE HERE THE USED RBFUNCTION)
    Eigen::VectorXd RBF_points_interpolated = interpolatorRBF.interpolate(parametersFORinterp, parameters, measurements);

    //  Use of RBF interpolation method with mean normalization
    RBFInterpolator interpolatorRBFmean(&RBFunctions::gaussian, scale_factor);
    Rescaling rescaling_mean;
    auto normalizedData_mean = rescaling_mean.meanNormalization(parameters, &parametersFORinterp);
    Eigen::VectorXd RBF_points_interpolated_mean = interpolatorRBFmean.interpolate(normalizedData_mean.second, normalizedData_mean.first, measurements);

    //  Use of RBF interpolation method with min-max normalization
    RBFInterpolator interpolatorRBFminmax(&RBFunctions::gaussian, scale_factor);
    Rescaling rescaling_minmax;
    auto normalizedData_minmax = rescaling_minmax.minMaxNormalization(parameters, &parametersFORinterp);
    Eigen::VectorXd RBF_points_interpolated_minmax = interpolatorRBFminmax.interpolate(normalizedData_minmax.second, normalizedData_minmax.first, measurements);

    //  Use of RBF interpolation method with z-score normalization
    RBFInterpolator interpolatorRBFzscore(&RBFunctions::gaussian, scale_factor);
    Rescaling rescaling_zscore;
    auto normalizedData_zscore = rescaling_zscore.zScoreNormalization(parameters, &parametersFORinterp);
    Eigen::VectorXd RBF_points_interpolated_zscore = interpolatorRBFzscore.interpolate(normalizedData_zscore.second, normalizedData_zscore.first, measurements);



    // Computing of the real values
    Eigen::VectorXd points_real(num_points);
    for (size_t i = 0; i < num_points; ++i)
    {
        points_real(i) = funToInterpolate(parametersFORinterp.row(i));
    }


    //////////////////////////////////////////////////
    ////////// PRINT OF RBF AND OLS RESULTS //////////
    //////////////////////////////////////////////////

    Eigen::VectorXd RBF_relative_errors = ((RBF_points_interpolated - points_real).array() / points_real.array().abs());
    Eigen::VectorXd RBF_relative_errors_mean = ((RBF_points_interpolated_mean - points_real).array() / points_real.array().abs());
    Eigen::VectorXd RBF_relative_errors_minmax = ((RBF_points_interpolated_minmax - points_real).array() / points_real.array().abs());
    Eigen::VectorXd RBF_relative_errors_zscore = ((RBF_points_interpolated_zscore - points_real).array() / points_real.array().abs());

    Eigen::VectorXd points = Eigen::VectorXd::LinSpaced(num_points, 1, num_points);
    printTable(points, points_real, RBF_points_interpolated, RBF_relative_errors, RBF_points_interpolated_mean, RBF_relative_errors_mean, RBF_points_interpolated_minmax, RBF_relative_errors_minmax, RBF_points_interpolated_zscore, RBF_relative_errors_zscore);
    
    bool EXPORT{true};
    plotData(RBF_relative_errors, RBF_relative_errors_mean, RBF_relative_errors_minmax, RBF_relative_errors_zscore, EXPORT);
    return 0;
}