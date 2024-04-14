#include <iostream>
#include <cmath>
#include <fstream>
#include <filesystem>

#include "RBFunctions.hpp"
#include "RBFinterpolator.hpp"


void plotData(const Eigen::MatrixXd& x, const Eigen::VectorXd& y_interpolatedRBF,
              const Eigen::MatrixXd& x_data, const Eigen::VectorXd& y_data, const bool EXPORT) {

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
        gnuplotScript << "set output \"" << "./plot/figures/case1_plot.svg" << "\"" << std::endl;
    }

    gnuplotScript << "set title \"" << title << "\"" << std::endl;
    gnuplotScript << "set xlabel \"x\"" << std::endl;
    gnuplotScript << "set ylabel \"interpolated f(x)\"" << std::endl;
    gnuplotScript << "plot \"" << "./plot/files/interpolated_points_RBF.dat" << "\" with lines title \"RBF interpolation\","
                  << " \"./plot/files/data.dat\" with points pointtype 7 title \"Regressors\"" << std::endl;
    
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
                    2,
                  3.5;
    
    // Create the vector with the known results corresponding to parameters
    Eigen::VectorXd measurements(num_measures);
    measurements <<   1,
                    0.2,
                    0.1;

    // Create the values of the parameters for which we want the interpolation
    Eigen::MatrixXd parametersFORinterp(num_points, num_params);

    int j{};
    for (size_t i = 0; i < num_points; ++i)
    {   
        parametersFORinterp(i,j) = inf + i * (sup - inf) / (num_points-1);
    }

    // Define the vectors to store the coefficients of the interpolations
    Eigen::VectorXd regressionRBF;

    // Use of RBF interpolation method
    double scale_factor{sqrt(0.5)};
    RBFInterpolator interpolatorRBF(&RBFunctions::gaussian, scale_factor);
    Eigen::VectorXd RBF_points_interpolated = interpolatorRBF.interpolate(parametersFORinterp, parameters, measurements, &regressionRBF);


    //////////////////////////////////////////////////
    ////////////// PRINT OF RBF RESULTS //////////////
    //////////////////////////////////////////////////

    bool PRINT{true};

    if (PRINT)
    {
        std::cout << "__________________USING RADIAL BASIS FUNCTIONS__________________" << std::endl;

        if (regressionRBF.size() != 0)
        {
            std::cout << "Coefficients: " << regressionRBF.transpose() << std::endl;
        }
        
        //std::cout << "Interpolated value: " << RBF_points_interpolated.transpose() << std::endl;
    }


    //////////////////////////////////////////////////
    /////////////// PLOT OF RBF RESULTS //////////////
    //////////////////////////////////////////////////

    bool EXPORT{true};
    plotData(parametersFORinterp, RBF_points_interpolated, parameters, measurements, EXPORT);

    return 0;
}