#include <iostream>
#include <cmath>
#include <fstream>
#include <filesystem>

#include "RBFunctions.hpp"
#include "RBFinterpolator.hpp"
#include "OLSinterpolator.hpp"


void plotData(const Eigen::MatrixXd& x, const Eigen::VectorXd& y_interpolatedRBF,
              const Eigen::VectorXd& y_interpolatedRBFpoly, const Eigen::MatrixXd& x_data,
              const Eigen::VectorXd& y_data, const bool EXPORT) {

    std::filesystem::create_directories("./plot/");
    std::filesystem::create_directories("./plot/files/");
    std::filesystem::create_directories("./plot/figures/");

    // RBF
    std::ofstream dataFileRBF("./plot/files/interpolated_points_RBF.dat");
    if (!dataFileRBF.is_open()) {
        std::cerr << "Error: Unable to open data file for RBF." << std::endl;
        return;
    }

    for (int i = 0; i < x.rows(); ++i) {
        dataFileRBF << x(i, 0) << " " << y_interpolatedRBF(i) << std::endl;
    }
    dataFileRBF.close();

    // RBF with polynomial
    std::ofstream dataFileRBFpoly("./plot/files/interpolated_points_RBF_poly.dat");
    if (!dataFileRBFpoly.is_open()) {
        std::cerr << "Error: Unable to open data file for RBF poly." << std::endl;
        return;
    }

    for (int i = 0; i < x.rows(); ++i) {
        dataFileRBFpoly << x(i, 0) << " " << y_interpolatedRBFpoly(i) << std::endl;
    }
    dataFileRBFpoly.close();

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


    std::string title{"1D interpolation : test case from [du Toit]"};


    // GNUplot commands
    std::ofstream gnuplotScript("./plot/files/plot_script.gnu");
    
    if (EXPORT)
    {
        // To export the plot in svg format
        gnuplotScript << "set terminal svg" << std::endl;
        gnuplotScript << "set output \"" << "./plot/figures/case4_plot.svg" << "\"" << std::endl;
    }

    gnuplotScript << "set title \"" << title << "\"" << std::endl;
    gnuplotScript << "set xlabel \"x\"" << std::endl;
    gnuplotScript << "set ylabel \"interpolated f(x)\"" << std::endl;
    gnuplotScript << "set key box" << std::endl;
    gnuplotScript << "set style line 1 lc rgb 'blue' lt 1 lw 2" << std::endl;
    gnuplotScript << "set style line 2 lc rgb 'orange' lt 1 lw 2" << std::endl;
    
    gnuplotScript << "plot \""
                  << "./plot/files/interpolated_points_RBF.dat" << "\" with lines linestyle 1 title \"RBF interpolation\","
                  << " \"./plot/files/interpolated_points_RBF_poly.dat" << "\" with lines linestyle 2 title \"RBFP interpolation\","
                  << " \"./plot/files/data.dat\" with points pointtype 7 title \"Regressors\""
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

    //////////////////////////////////////////////////
    /////// TEST : 1D from Wilna du Toit p.4-5 ///////
    //////////////////////////////////////////////////

    size_t num_params{1};      // 1D scattered data : y = f(x)
    size_t num_measures{3};    // 3 known points
    size_t num_points{100};    // number of points to interpolate
    
    // [inf, sup] is the interval of definition of the parameter
    double inf{-2}; 
    double sup{7};
    
    // Create the matrix with the known parameters
    Eigen::MatrixXd parameters(num_measures, num_params);
    parameters <<   1,
                    3,
                  3.5;
    
    // Create the vector with the known results corresponding to parameters
    Eigen::VectorXd measurements(num_measures);
    measurements <<   1,
                    0.2,
                    0.1;

    // Create the values of the parameters for which we want the interpolation
    Eigen::MatrixXd parametersFORinterp(num_points, num_params);

    for (size_t i = 0; i < num_points; ++i)
    {   
        parametersFORinterp(i,0) = inf + i * (sup - inf) / (num_points-1);
    }

    // Define the vectors to store the coefficients of the interpolations
    Eigen::VectorXd regressionRBF;
    Eigen::VectorXd regressionRBFpoly;

    // Use of RBF interpolation method
    double scale_factor{sqrt(0.5)};
    RBFInterpolator interpolatorRBF(&RBFunctions::gaussian, scale_factor);
    Eigen::VectorXd RBF_points_interpolated = interpolatorRBF.interpolate(parametersFORinterp, parameters, measurements, &regressionRBF);

    // Use of RBF interpolation method augmented with polynomial
    bool normalizeRBF{false};
    bool polynomialRBF{true};
    RBFInterpolator interpolatorRBFpoly(&RBFunctions::gaussian, scale_factor, normalizeRBF, polynomialRBF);
    Eigen::VectorXd RBF_poly_points_interpolated = interpolatorRBFpoly.interpolate(parametersFORinterp, parameters, measurements, &regressionRBFpoly);


    //////////////////////////////////////////////////
    ////////////// PRINT OF RBF RESULTS //////////////
    //////////////////////////////////////////////////

    bool PRINT{true};

    if (PRINT)
    {
        

        if (regressionRBF.size() != 0)
        {
            std::cout << "_USING RADIAL BASIS FUNCTIONS_" << std::endl;
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
        
        std::cout << "\n" << std::endl;
        //std::cout << "RBF interpolated value: " << RBF_points_interpolated.transpose() << std::endl;

        if (regressionRBFpoly.size() != 0)
        {
            std::cout << "_USING POLYNOMIAL RADIAL BASIS FUNCTIONS_" << std::endl;
            std::cout << "______________________________" << std::endl;

            std::cout << "||" << std::setw(5) << ""
                      << std::setw(21) << std::left << "|| Weights"
                      << "||" << std::endl;

            std::cout << "______________________________" << std::endl;
            for (int i = 0; i < regressionRBFpoly.size(); ++i) {
                std::cout << "|| " << std::setw(4) << std::left << i
                          << "|| " << std::setw(18) << std::left << regressionRBFpoly(i)
                          << "|| " << std::endl;
            }
            std::cout << "______________________________" << std::endl;

        }
    }


    //////////////////////////////////////////////////
    /////////////// PLOT OF RBF RESULTS //////////////
    //////////////////////////////////////////////////

    bool EXPORT{true};
    plotData(parametersFORinterp, RBF_points_interpolated, RBF_poly_points_interpolated, parameters, measurements, EXPORT);

    return 0;
}