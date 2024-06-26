#include <Eigen/Dense>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <iomanip>
#include "RBFunctions.hpp"

//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// TEST CASE 0 ///////////////////////////////////
/////////// PLOT OF RADIAL BASIS FUNCTIONS FOR DIFFERENT SCALE FACTORS ///////////
//////////////////////////////////////////////////////////////////////////////////

int main() {

    // Definition of the scale factors.
    Eigen::VectorXd scale_factor(9);
    scale_factor << 0.5, 0.75, 1, 1.25, 1.50, 2, 3, 4, 5;

    // Definition of the x axis.
    Eigen::VectorXd tab_x = Eigen::VectorXd::LinSpaced(100, -10, 10);

    // Definition of y axis.
    Eigen::MatrixXd gaussian_mat(scale_factor.size(), tab_x.size());
    Eigen::MatrixXd multiquadratic_mat(scale_factor.size(), tab_x.size());
    Eigen::MatrixXd inverseMultiquadratic_mat(scale_factor.size(), tab_x.size());
    Eigen::MatrixXd thinPlateSpline_mat(scale_factor.size(), tab_x.size());

    // Computation of the radial basis functions.
    for (size_t i = 0; i < scale_factor.size(); ++i)
    {
        for (size_t j = 0; j < tab_x.size(); ++j)
        {
            gaussian_mat(i,j) = RBFunctions::gaussian(tab_x(j), scale_factor(i));
            multiquadratic_mat(i,j) = RBFunctions::multiquadratic(tab_x(j), scale_factor(i));
            inverseMultiquadratic_mat(i,j) = RBFunctions::inverseMultiquadratic(tab_x(j), scale_factor(i));
            thinPlateSpline_mat(i,j) = RBFunctions::thinPlateSpline(tab_x(j), scale_factor(i));
        }
    }

    // Boolean used to export in .svg the plots.
    bool EXPORT{true};

    std::filesystem::create_directories("./plot/");
    std::filesystem::create_directories("./plot/files/");
    std::filesystem::create_directories("./plot/figures/");


    //////////////////////////////
    ////////// GAUSSIAN //////////
    //////////////////////////////

    std::ofstream dataFileGauss("./plot/files/data_gaussian.dat");
    if (!dataFileGauss.is_open()) {
        std::cerr << "Error: Unable to open data_gaussian.dat file" << std::endl;
        return 1;
    }

    for (size_t i = 0; i < scale_factor.size(); ++i)
    {
        for (size_t j = 0; j < tab_x.size(); ++j)
        {
            dataFileGauss << tab_x(j) << " " << gaussian_mat(i,j) << " " << scale_factor(i) << " " << std::endl;
        }
        dataFileGauss << "\n\n";
    }
    dataFileGauss.close();

    std::ofstream gnuplotScriptGauss("./plot/files/plot_gaussian.gnu");
    if (!gnuplotScriptGauss.is_open()) {
        std::cerr << "Error: Unable to open GNUplot script for gaussian file." << std::endl;
        return 1;
    }

    if (EXPORT)
    {
        gnuplotScriptGauss << "set terminal svg" << std::endl;
        gnuplotScriptGauss << "set output \"" << "./plot/figures/gaussian_plot.svg" << "\"" << std::endl;
    }

    gnuplotScriptGauss << "set title \"Gaussian basis function\"" << std::endl;
    gnuplotScriptGauss << "set key box" << std::endl;
    gnuplotScriptGauss << "plot ";
    for (size_t i = 0; i < scale_factor.size(); ++i)
    {
        if (i > 0) {
            gnuplotScriptGauss << ", ";
        }
        gnuplotScriptGauss << "'./plot/files/data_gaussian.dat' index " << i << " using 1:2:3 with lines title sprintf('r_0 = %.2f', " << scale_factor(i) << ")";
    }
    gnuplotScriptGauss << std::endl;

    if (EXPORT)
    {   
        // Need to replot the plot to see it on the screen.
        gnuplotScriptGauss << "set terminal wxt" << std::endl;
        gnuplotScriptGauss << "replot" << std::endl;
    }

    gnuplotScriptGauss.close();


    //////////////////////////////
    /////// MULTIQUADRATIC ///////
    //////////////////////////////

    std::ofstream dataFilemultiquad("./plot/files/data_multiquad.dat");
    if (!dataFilemultiquad.is_open()) {
        std::cerr << "Error: Unable to open data_multiquad.dat file." << std::endl;
        return 1;
    }

    for (size_t i = 0; i < scale_factor.size(); ++i)
    {
        for (size_t j = 0; j < tab_x.size(); ++j)
        {
            dataFilemultiquad << tab_x(j) << " " << multiquadratic_mat(i,j) << " " << scale_factor(i) << " "  << std::endl;
        }
        dataFilemultiquad << "\n\n";
    }
    dataFilemultiquad.close();

    std::ofstream gnuplotScriptmultiquad("./plot/files/plot_multiquad.gnu");
    if (!gnuplotScriptmultiquad.is_open()) {
        std::cerr << "Error: Unable to open GNUplot script for multiquad file." << std::endl;
        return 1;
    }

    if (EXPORT)
    {
        gnuplotScriptmultiquad << "set terminal svg" << std::endl;
        gnuplotScriptmultiquad << "set output \"" << "./plot/figures/multiquad_plot.svg" << "\"" << std::endl;
    }

    gnuplotScriptmultiquad << "set title \"Multiquadratic basis function\"" << std::endl;
    gnuplotScriptmultiquad << "set key box bottom right" << std::endl;
    gnuplotScriptmultiquad << "plot ";
    for (size_t i = 0; i < scale_factor.size(); ++i)
    {
        if (i > 0) {
            gnuplotScriptmultiquad << ", ";
        }
        gnuplotScriptmultiquad << "'./plot/files/data_multiquad.dat' index " << i << " using 1:2:3 with lines title sprintf('r_0 = %.2f', " << scale_factor(i) << ")";
    }
    gnuplotScriptmultiquad << std::endl;

    if (EXPORT)
    {   
        gnuplotScriptmultiquad << "set terminal wxt" << std::endl;
        gnuplotScriptmultiquad << "replot" << std::endl;
    }
    
    gnuplotScriptmultiquad.close();


    //////////////////////////////
    ///// INVMULTIQUADRATIC //////
    //////////////////////////////

    std::ofstream dataFileinvmultiquad("./plot/files/data_invmultiquad.dat");
    if (!dataFileinvmultiquad.is_open()) {
        std::cerr << "Error: Unable to open data_invmultiquad.dat file." << std::endl;
        return 1;
    }

    for (size_t i = 0; i < scale_factor.size(); ++i)
    {
        for (size_t j = 0; j < tab_x.size(); ++j)
        {
            dataFileinvmultiquad << tab_x(j) << " " << inverseMultiquadratic_mat(i,j) << " " << scale_factor(i) << " "  << std::endl;
        }
        dataFileinvmultiquad << "\n\n";
    }
    dataFileinvmultiquad.close();

    std::ofstream gnuplotScriptinvmultiquad("./plot/files/plot_invmultiquad.gnu");
    if (!gnuplotScriptinvmultiquad.is_open()) {
        std::cerr << "Error: Unable to open GNUplot script for invmultiquad file." << std::endl;
        return 1;
    }

    if (EXPORT)
    {
        gnuplotScriptinvmultiquad << "set terminal svg" << std::endl;
        gnuplotScriptinvmultiquad << "set output \"" << "./plot/figures/invmultiquad_plot.svg" << "\"" << std::endl;
    }

    gnuplotScriptinvmultiquad << "set title \"Inverse multiquadratic basis function\"" << std::endl;
    gnuplotScriptinvmultiquad << "set key box top right" << std::endl;
    gnuplotScriptinvmultiquad << "plot ";
    for (size_t i = 0; i < scale_factor.size(); ++i)
    {
        if (i > 0) {
            gnuplotScriptinvmultiquad << ", ";
        }
        gnuplotScriptinvmultiquad << "'./plot/files/data_invmultiquad.dat' index " << i << " using 1:2:3 with lines title sprintf('r_0 = %.2f', " << scale_factor(i) << ")";
    }
    gnuplotScriptinvmultiquad << std::endl;

    if (EXPORT)
    {   
        gnuplotScriptinvmultiquad << "set terminal wxt" << std::endl;
        gnuplotScriptinvmultiquad << "replot" << std::endl;
    }
    
    gnuplotScriptinvmultiquad.close();


    //////////////////////////////
    ///// Thin Plate Spline //////
    //////////////////////////////

    std::ofstream dataFilethinPlateSpline("./plot/files/data_thinPlateSpline.dat");
    if (!dataFilethinPlateSpline.is_open()) {
        std::cerr << "Error: Unable to open data_thinPlateSpline.dat file." << std::endl;
        return 1;
    }

    for (size_t i = 0; i < scale_factor.size(); ++i)
    {
        for (size_t j = 0; j < tab_x.size(); ++j)
        {
            dataFilethinPlateSpline << tab_x(j) << " " << thinPlateSpline_mat(i,j) << " " << scale_factor(i) << " "  << std::endl;
        }
        dataFilethinPlateSpline << "\n\n";
    }
    dataFilethinPlateSpline.close();

    std::ofstream gnuplotScriptthinPlateSpline("./plot/files/plot_thinPlateSpline.gnu");
    if (!gnuplotScriptthinPlateSpline.is_open()) {
        std::cerr << "Error: Unable to open GNUplot script for thinPlateSpline file." << std::endl;
        return 1;
    }

    if (EXPORT)
    {
        gnuplotScriptthinPlateSpline << "set terminal svg" << std::endl;
        gnuplotScriptthinPlateSpline << "set output \"" << "./plot/figures/thinPlateSpline_plot.svg" << "\"" << std::endl;
    }

    gnuplotScriptthinPlateSpline << "set title \"Thin plate spline basis function\"" << std::endl;
    gnuplotScriptthinPlateSpline << "set key box top left" << std::endl;
    gnuplotScriptthinPlateSpline << "plot [0:10]";
    for (size_t i = 0; i < scale_factor.size(); ++i)
    {
        if (i > 0) {
            gnuplotScriptthinPlateSpline << ", ";
        }
        gnuplotScriptthinPlateSpline << "'./plot/files/data_thinPlateSpline.dat' index " << i << " using 1:2:3 with lines title sprintf('r_0 = %.2f', " << scale_factor(i) << ")";
    }
    gnuplotScriptthinPlateSpline << std::endl;

    if (EXPORT)
    {   
        gnuplotScriptthinPlateSpline << "set terminal wxt" << std::endl;
        gnuplotScriptthinPlateSpline << "replot" << std::endl;
    }

    gnuplotScriptthinPlateSpline.close();




    system("gnuplot -persist ./plot/files/plot_gaussian.gnu");
    system("gnuplot -persist ./plot/files/plot_multiquad.gnu");
    system("gnuplot -persist ./plot/files/plot_invmultiquad.gnu");
    system("gnuplot -persist ./plot/files/plot_thinPlateSpline.gnu");

    return 0;
}