#include <iostream>
#include <cmath>
#include <fstream>
#include <filesystem>

#include "RBFunctions.hpp"
#include "RBFinterpolator.hpp"
#include "OLSinterpolator.hpp"


void plotData(const Eigen::MatrixXd& x, const Eigen::VectorXd& y_interpolatedRBF, const Eigen::VectorXd& y_interpolatedOLS, const Eigen::MatrixXd& x_data, const Eigen::VectorXd& y_data, const bool EXPORT) {
    std::filesystem::create_directories("./plot/");
    std::filesystem::create_directories("./plot/files/");
    std::filesystem::create_directories("./plot/figures/");

    // RBF
    std::ofstream dataFileRBF("./plot/files/interpolated_points_RBF.dat");
    if (!dataFileRBF.is_open()) {
        std::cerr << "Error: Unable to open data file." << std::endl;
        return;
    }

    for (int i = 0; i < x.rows(); ++i) {
        dataFileRBF << x(i, 0) << " " << y_interpolatedRBF(i) << std::endl;
    }
    dataFileRBF.close();

    // OLS
    std::ofstream dataFileOLS("./plot/files/interpolated_points_OLS.dat");
    if (!dataFileOLS.is_open()) {
        std::cerr << "Error: Unable to open data file." << std::endl;
        return;
    }

    for (int i = 0; i < x.rows(); ++i) {
        dataFileOLS << x(i, 0) << " " << y_interpolatedOLS(i) << std::endl;
    }
    dataFileOLS.close();

    // DATA FROM [du Toit]
    std::ofstream dataFileData("./plot/files/data.dat");
    if (!dataFileData.is_open()) {
        std::cerr << "Error: Unable to open data file for data." << std::endl;
        return;
    }

    for (int i = 0; i < x_data.rows(); ++i) {
        dataFileData << x_data(i, 0) << " " << y_data(i) << std::endl;
    }
    dataFileData.close();

    std::string title{"1D interpolation : test case of a linear function y = a*x + b"};

    // GNUplot commands
    std::ofstream gnuplotScript("./plot/files/plot_script.gnu");
    
    if (EXPORT)
    {
        // To export the plot in svg format
        gnuplotScript << "set terminal svg" << std::endl;
        gnuplotScript << "set output \"" << "./plot/figures/case2_plot.svg" << "\"" << std::endl;
    }

    gnuplotScript << "set title \"" << title << "\"" << std::endl;
    gnuplotScript << "set xlabel \"x\"" << std::endl;
    gnuplotScript << "set ylabel \"interpolated f(x)\"" << std::endl;
    gnuplotScript << "plot \"" << "./plot/files/interpolated_points_RBF.dat" << "\" with lines title \"RBF interpolation\","
                << " \"" << "./plot/files/interpolated_points_OLS.dat" << "\" with lines title \"OLS interpolation\","
                << " \"" << "./plot/files/data.dat" << "\" with points pointtype 7 title \"Regressors\"" 
                << std::endl; // Plot command ends here, before the legend key
    gnuplotScript << "set key bottom right" << std::endl; // Setting the position of the legend
        
    
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

    //////////////////////////////////////////////////
    ////// TEST : 1D OLS with a linear function //////
    //////////////////////////////////////////////////

    size_t num_params{1};      // 1D scattered data : y = f(x)
    size_t num_measures{5};    // 3 known points
    size_t num_points{10};     // number of points to interpolate
    
    auto funToInterpolate = [](double& a, double& b, double& x) {
        return a*x + b;
    };

    // y = 0.5 * x - 4.3
    double a{0.5};
    double b{-4.3};

    // Create the matrix with the known parameters
    Eigen::MatrixXd parameters(num_measures, num_params);
    parameters <<  -2,
                  3.7,
                  0.1,
                   -6,
                 18.2;
    
    // Create the vector with the known results corresponding to parameters
    Eigen::VectorXd measurements(num_measures);
    
    for (size_t i = 0; i < num_measures; ++i)
    {
        measurements(i) = funToInterpolate(a, b, parameters(i));
    }    

    // Create the values of the parameters for which we want the interpolation
    Eigen::MatrixXd parametersFORinterp(num_points, num_params);

    double inf{-10};
    double sup{20};

    int j{};
    for (size_t i = 0; i < num_points; ++i)
    {   
        parametersFORinterp(i,j) = inf + i * (sup - inf) / (num_points-1);
    }

    // Define the vectors to store the coefficients of the interpolations
    Eigen::VectorXd regressionRBF;
    Eigen::VectorXd regressionOLS;

    // Use of RBF interpolation method
    double scale_factor{0};
    RBFInterpolator interpolatorRBF(&RBFunctions::multiquadratic, scale_factor);
    Eigen::VectorXd RBF_points_interpolated = interpolatorRBF.interpolate(parametersFORinterp, parameters, measurements, &regressionRBF);

    // Use of OLS approximation method
    OLSInterpolator interpolatorOLS;
    Eigen::VectorXd OLS_points_interpolated = interpolatorOLS.interpolate(parametersFORinterp, parameters, measurements, &regressionOLS);


    //////////////////////////////////////////////////
    ////////// PRINT OF RBF AND OLS RESULTS //////////
    //////////////////////////////////////////////////

    bool PRINT{true};

    if (PRINT)
    {
        std::cout << "_USING RADIAL BASIS FUNCTIONS_" << std::endl;

        if (regressionRBF.size() != 0)

        {
            std::cout << "______________________________" << std::endl;

            std::cout << "||" << std::setw(5) << ""
                      << std::setw(21) << std::left << "|| Weights"
                      << "||" << std::endl;

            std::cout << "______________________________" << std::endl;
            for (int i = 0; i < regressionRBF.size(); ++i) {
                std::cout << "|| " << std::setw(4) << std::left << i
                          << "|| " << std::setw(18) << std::left << regressionRBF(i)
                          << "|| " << std::endl;
            }
            std::cout << "______________________________" << std::endl;

        }
        
        //std::cout << "Interpolated value: " << RBF_points_interpolated.transpose() << std::endl;

        std::cout << std::endl;

        std::cout << "_USING ORDINARY LEAST SQUARES_" << std::endl;

        if (regressionOLS.size() != 0)
        {
                        std::cout << "______________________________" << std::endl;

            std::cout << "||" << std::setw(5) << ""
                      << std::setw(21) << std::left << "|| Coefficients"
                      << "||" << std::endl;

            std::cout << "______________________________" << std::endl;
            for (int i = 0; i < regressionOLS.size(); ++i) {
                std::cout << "|| " << std::setw(4) << std::left << i
                          << "|| " << std::setw(18) << std::left << regressionOLS(i)
                          << "|| " << std::endl;
            }
            std::cout << "______________________________" << std::endl;
        }

        //std::cout << "Interpolated value: " << OLS_points_interpolated.transpose() << std::endl;
    }


    //////////////////////////////////////////////////
    /////////// PLOT OF RBF AND OLS RESULTS //////////
    //////////////////////////////////////////////////

    bool EXPORT{true};
    plotData(parametersFORinterp, RBF_points_interpolated, OLS_points_interpolated, parameters, measurements, EXPORT);

    return 0;
}