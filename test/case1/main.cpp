#include <iostream>
#include <cmath>
#include <fstream>
#include <filesystem>

#include "RBFunctions.hpp"
#include "RBFinterpolator.hpp"
#include "OLSinterpolator.hpp"


void plotData(const Eigen::MatrixXd& x, const Eigen::VectorXd& y_interpolatedRBF,
              const Eigen::VectorXd& y_interpolatedRBFnorm, const Eigen::VectorXd& y_interpolatedRBFpoly,
              const Eigen::VectorXd& y_interpolatedOLS, const Eigen::MatrixXd& x_data, const Eigen::VectorXd& y_data,
              const Eigen::VectorXd& point1, const Eigen::VectorXd& point2, const Eigen::VectorXd& point3, const bool EXPORT) {

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

    // RBF normalized
    std::ofstream dataFileRBFnorm("./plot/files/interpolated_points_RBF_norm.dat");
    if (!dataFileRBFnorm.is_open()) {
        std::cerr << "Error: Unable to open data file for RBF norm." << std::endl;
        return;
    }

    for (int i = 0; i < x.rows(); ++i) {
        dataFileRBFnorm << x(i, 0) << " " << y_interpolatedRBFnorm(i) << std::endl;
    }
    dataFileRBFnorm.close();

    // RBF polynomial
    std::ofstream dataFileRBFpoly("./plot/files/interpolated_points_RBF_poly.dat");
    if (!dataFileRBFpoly.is_open()) {
        std::cerr << "Error: Unable to open data file for RBF poly." << std::endl;
        return;
    }

    for (int i = 0; i < x.rows(); ++i) {
        dataFileRBFpoly << x(i, 0) << " " << y_interpolatedRBFpoly(i) << std::endl;
    }
    dataFileRBFpoly.close();

    // OLS
    std::ofstream dataFileOLS("./plot/files/interpolated_points_OLS.dat");
    if (!dataFileOLS.is_open()) {
        std::cerr << "Error: Unable to open data file for OLS." << std::endl;
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

    // Point 1 contribution
    std::ofstream point1Data("./plot/files/point1.dat");
    if (!point1Data.is_open()) {
        std::cerr << "Error: Unable to open data file for point1." << std::endl;
        return;
    }

    for (int i = 0; i < x.rows(); ++i) {
        point1Data << x(i, 0) << " " << point1(i) << std::endl;
    }
    point1Data.close();

    // Point 2 contribution
    std::ofstream point2Data("./plot/files/point2.dat");
    if (!point2Data.is_open()) {
        std::cerr << "Error: Unable to open data file for point2." << std::endl;
        return;
    }

    for (int i = 0; i < x.rows(); ++i) {
        point2Data << x(i, 0) << " " << point2(i) << std::endl;
    }
    point2Data.close();

    // Point 3 contribution
    std::ofstream point3Data("./plot/files/point3.dat");
    if (!point3Data.is_open()) {
        std::cerr << "Error: Unable to open data file for point3." << std::endl;
        return;
    }

    for (int i = 0; i < x.rows(); ++i) {
        point3Data << x(i, 0) << " " << point3(i) << std::endl;
    }
    point3Data.close();



    std::string title{"1D interpolation : test case from [Wilna du Toit]"};


    // GNUplot commands
    std::ofstream gnuplotScript1("./plot/files/plot_script1.gnu");
    
    if (EXPORT)
    {
        // To export the plot in svg format
        gnuplotScript1 << "set terminal svg" << std::endl;
        gnuplotScript1 << "set output \"" << "./plot/figures/case1_plot1.svg" << "\"" << std::endl;
    }

    gnuplotScript1 << "set title \"" << title << "\"" << std::endl;
    gnuplotScript1 << "set xlabel \"x\"" << std::endl;
    gnuplotScript1 << "set ylabel \"interpolated f(x)\"" << std::endl;
    gnuplotScript1 << "set key box" << std::endl;
    gnuplotScript1 << "set style line 1 lc rgb 'blue' lt 1 lw 2" << std::endl;
    gnuplotScript1 << "set style line 2 lc rgb 'orange' dt 2 lw 2" << std::endl;
    gnuplotScript1 << "set style line 3 lc rgb 'purple' dt 2 lw 2" << std::endl;
    gnuplotScript1 << "set style line 4 lc rgb 'red' dt 2 lw 2" << std::endl;
    
    gnuplotScript1 << "plot \""
                  << "./plot/files/data.dat" << "\" with points pointtype 7 title \"Regressors\"," 
                  << " \"./plot/files/interpolated_points_RBF.dat\" with lines linestyle 1 title \"RBF interpolation\","   
                  << " \"./plot/files/point1.dat\" with lines linestyle 2 title \"1st point RBF contribution\","
                  << " \"./plot/files/point2.dat\" with lines linestyle 3 title \"2nd point RBF contribution\","
                  << " \"./plot/files/point3.dat\" with lines linestyle 4 title \"3rd point RBF contribution\""
                  << std::endl;


    
    if (EXPORT)
    {   
        // The export "cancel" the plot
        // Replot to see it in a window
        gnuplotScript1 << "set terminal wxt" << std::endl;
        gnuplotScript1 << "replot" << std::endl;
    }

    gnuplotScript1.close();




    // GNUplot commands
    std::ofstream gnuplotScript2("./plot/files/plot_script2.gnu");
    
    if (EXPORT)
    {
        // To export the plot in svg format
        gnuplotScript2 << "set terminal svg" << std::endl;
        gnuplotScript2 << "set output \"" << "./plot/figures/case1_plot2.svg" << "\"" << std::endl;
    }

    gnuplotScript2 << "set title \"" << title << "\"" << std::endl;
    gnuplotScript2 << "set xlabel \"x\"" << std::endl;
    gnuplotScript2 << "set ylabel \"interpolated f(x)\"" << std::endl;
    gnuplotScript2 << "set key box" << std::endl;
    gnuplotScript2 << "set style line 1 lc rgb 'blue' lt 1 lw 2" << std::endl;
    gnuplotScript2 << "set style line 2 lc rgb 'orange' lt 1 lw 2" << std::endl;
    gnuplotScript2 << "set style line 3 lc rgb 'red' lt 1 lw 2" << std::endl;
    gnuplotScript2 << "set style line 4 lc rgb 'black' dt 2 lw 2" << std::endl;
    
    gnuplotScript2 << "plot \""
                  << "./plot/files/data.dat" << "\" with points pointtype 7 title \"Regressors\","
                  << " \"./plot/files/interpolated_points_RBF.dat" << "\" with lines linestyle 1 title \"RBF interpolation\","
                  << " \"./plot/files/interpolated_points_RBF_norm.dat" << "\" with lines linestyle 2 title \"NRBF interpolation\","
                  << " \"./plot/files/interpolated_points_RBF_poly.dat" << "\" with lines linestyle 3 title \"RBFP interpolation\","
                  << " \"./plot/files/interpolated_points_OLS.dat" << "\" with lines linestyle 4 title \"OLS interpolation\""
                  << std::endl;


    
    if (EXPORT)
    {   
        // The export "cancel" the plot
        // Replot to see it in a window
        gnuplotScript2 << "set terminal wxt" << std::endl;
        gnuplotScript2 << "replot" << std::endl;
    }

    gnuplotScript2.close();


    // Execute GNUplot
    system("gnuplot -persist ./plot/files/plot_script1.gnu");
    system("gnuplot -persist ./plot/files/plot_script2.gnu");
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
    Eigen::VectorXd regressionRBFnorm;
    Eigen::VectorXd regressionRBFpoly;
    Eigen::VectorXd regressionOLS;

    // Use of RBF interpolation method
    double scale_factor{sqrt(0.5)};
    RBFInterpolator interpolatorRBF(&RBFunctions::gaussian, scale_factor);
    Eigen::VectorXd RBF_points_interpolated = interpolatorRBF.interpolate(parametersFORinterp, parameters, measurements, &regressionRBF);

    // Use of normalized RBF interpolation method
    bool normalizeRBF{true};
    RBFInterpolator interpolatorRBFnorm(&RBFunctions::gaussian, scale_factor, normalizeRBF);
    Eigen::VectorXd RBF_norm_points_interpolated = interpolatorRBFnorm.interpolate(parametersFORinterp, parameters, measurements, &regressionRBFnorm);

    // Use of RBF augmented with polynomials
    normalizeRBF=false;
    bool polynomialRBF{true};
    RBFInterpolator interpolatorRBFpoly(&RBFunctions::gaussian, scale_factor, normalizeRBF, polynomialRBF);
    Eigen::VectorXd RBF_poly_points_interpolated = interpolatorRBFpoly.interpolate(parametersFORinterp, parameters, measurements, &regressionRBFpoly);


    // Use of OLS interpolation method
    OLSInterpolator interpolatorOLS;
    Eigen::VectorXd OLS_points_interpolated = interpolatorOLS.interpolate(parametersFORinterp, parameters, measurements, &regressionOLS);


    // Datapoint's contribution to the interpolant

    Eigen::VectorXd point1(num_points);
    Eigen::VectorXd point2(num_points);
    Eigen::VectorXd point3(num_points);


    for (size_t i = 0; i < num_points; ++i)
    {
        point1(i) = regressionRBF(0) * RBFunctions::gaussian((parametersFORinterp(i,0) - parameters(0,0)) , scale_factor);
        point2(i) = regressionRBF(1) * RBFunctions::gaussian((parametersFORinterp(i,0) - parameters(1,0)) , scale_factor);
        point3(i) = regressionRBF(2) * RBFunctions::gaussian((parametersFORinterp(i,0) - parameters(2,0)) , scale_factor);
    }


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

        if (regressionRBFnorm.size() != 0)
        {
            std::cout << "_USING NORMALIZED RADIAL BASIS FUNCTIONS_" << std::endl;
            std::cout << "______________________________" << std::endl;

            std::cout << "||" << std::setw(5) << ""
                      << std::setw(21) << std::left << "|| Weights"
                      << "||" << std::endl;

            std::cout << "______________________________" << std::endl;
            for (int i = 0; i < regressionRBFnorm.size(); ++i) {
                std::cout << "|| " << std::setw(4) << std::left << i
                          << "|| " << std::setw(18) << std::left << regressionRBFnorm(i)
                          << "|| " << std::endl;
            }
            std::cout << "______________________________" << std::endl;

        }

        std::cout << "\n" << std::endl;

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
        
        std::cout << "\n" << std::endl;

        if (regressionOLS.size() != 0)
        {
            std::cout << "_USING ORDINARY LEAST SQUARES_" << std::endl;
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
    }


    //////////////////////////////////////////////////
    /////////////// PLOT OF RBF RESULTS //////////////
    //////////////////////////////////////////////////

    bool EXPORT{true};
    plotData(parametersFORinterp, RBF_points_interpolated, RBF_norm_points_interpolated, RBF_poly_points_interpolated, OLS_points_interpolated, parameters, measurements, point1, point2, point3, EXPORT);

    return 0;
}