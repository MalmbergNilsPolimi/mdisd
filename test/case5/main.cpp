#include <iostream>
#include <cmath>
#include <iomanip>
#include <random>
#include <fstream>
#include <filesystem>
#include <unistd.h>
#include <chrono>

#include "RBFunctions.hpp"
#include "RBFinterpolator.hpp"


//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// TEST CASE 5 ///////////////////////////////////
///////////////////////// ABOUT CONVERGENCE OF RBF ERROS /////////////////////////
//////////////////////////////////////////////////////////////////////////////////


template<typename T>
T generateRandomNumber(T inf, T sup) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dis(inf, sup);
    return dis(gen);
}


auto rbfunForCompact(double r, double r0) {
    (void)r0;
    double diff{};
    if ((1-r) >= 0)
    {
        diff = 1-r;
    } else {
        diff = 0;
    }
    return std::pow(diff, 6)*(35*r*r + 18*r + 3);
}


void plotData(const Eigen::VectorXi& KNOWN_POINTS, const Eigen::VectorXd& MSE, const std::string dim, const std::string error, const bool EXPORT) {

    std::filesystem::create_directories("./plot/");
    std::filesystem::create_directories("./plot/files/");
    std::filesystem::create_directories("./plot/figures/");

    std::ofstream dataFile("./plot/files/"+error+"_RBF_"+dim+".dat");
    if (!dataFile.is_open()) {
        std::cerr << "Error: Unable to open " << error << " data file." << std::endl;
        return;
    }

    for (int j = 0; j < KNOWN_POINTS.size(); ++j) {
        dataFile << KNOWN_POINTS(j) << " " << MSE(j) << std::endl;
    }
    dataFile.close();

    std::string title{"Convergence of the RBF "+error+" error"};
    std::string label_y{error+" error on interpolated points"};

    std::ofstream gnuplotScript("./plot/files/plot_script_"+error+"_"+dim+".gnu");

    if (EXPORT)
    {
        // To export the plot in svg format
        gnuplotScript << "set terminal svg" << std::endl;
        gnuplotScript << "set output \"" << "./plot/figures/plot_"+error+"error_"+dim+".svg" << "\"" << std::endl;
    }

    gnuplotScript << "set title \"" << title << "\"" << std::endl;
    gnuplotScript << "set xlabel \"number of known points\"" << std::endl;
    gnuplotScript << "set ylabel \"" << label_y << "\"" << std::endl;
    gnuplotScript << "set logscale y\n";
    gnuplotScript << "set format y '%.1e'\n";
    gnuplotScript << "set key box" << std::endl;
    
    gnuplotScript << "plot './plot/files/"+error+"_RBF_"+dim+".dat' with lines lw 2 title \"";
        
    gnuplotScript << std::endl;
    
    if (EXPORT)
    {   
        gnuplotScript << "set terminal wxt" << std::endl;
        gnuplotScript << "replot" << std::endl;
    }

    gnuplotScript.close();

    std::string path{"gnuplot -persist ./plot/files/plot_script_"+error+"_"+dim+".gnu"};
    system(path.c_str());
}

void plotRelativeError(const Eigen::MatrixXd& RE, const std::string known_points, const std::string dimensions){
    // Create directories
    std::filesystem::create_directories("./plot/");
    std::filesystem::create_directories("./plot/files/");
    std::filesystem::create_directories("./plot/figures/");

    // Open a file to write the data
    std::ofstream data_file("./plot/files/relative_error_"+dimensions+"_"+known_points+"knownpoints.dat");
    if (!data_file) {
        std::cerr << "Error opening file to write data.\n";
    }

    // Generate the data
    for (int i = 0; i < RE.rows(); ++i) {
        data_file << RE(i,0) << " " << RE(i,1) << " " << RE(i,2) << "\n";
        data_file << "\n";
    }

    data_file.close();

    // Create a gnuplot script file
    std::ofstream gnuplot_script("./plot/files/relative_error_"+dimensions+"_"+known_points+"knownpoints.gnuplot");
    if (!gnuplot_script) {
        std::cerr << "Error opening file to write gnuplot script.\n";
    }

    gnuplot_script << "set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n";
    gnuplot_script << "set output './plot/figures/relative_error_"+dimensions+"_"+known_points+"knownpoints.png'\n";
    gnuplot_script << "set title 'Relative Error in function of the position of the point'\n";
    gnuplot_script << "set xlabel 'x'\n";
    gnuplot_script << "set ylabel 'y'\n";
    gnuplot_script << "set zlabel 'relative error'\n";
    gnuplot_script << "set pm3d at s\n";
    gnuplot_script << "set palette defined (0 'blue', 1 'green', 2 'yellow', 3 'red')\n";
    gnuplot_script << "set xrange [0:1]\n";
    gnuplot_script << "set yrange [1:0]\n";
    // gnuplot_script << "set zrange [0:1.6]\n";
    gnuplot_script << "set view 60, 30, 1, 1\n";
    gnuplot_script << "set xyplane at 0\n";
    gnuplot_script << "splot './plot/files/relative_error_"+dimensions+"_"+known_points+"knownpoints.dat' using 1:2:3 with points pointtype 7 pointsize 1 palette title 'Points'\n";

    gnuplot_script.close();

    // Run gnuplot with the script
    std::string path{"gnuplot ./plot/files/relative_error_"+dimensions+"_"+known_points+"knownpoints.gnuplot"};
    system(path.c_str());
}

int main() {

    ////////////////////////////////////////////////
    // TEST : Using results from Lazzaro's paper ///
    ////////////////////////////////////////////////

    ///// PLOT OF FRANKE'S FUNCTION (function to interpolate)

    // Definition of the function to interpolate.
    auto Franke_function = [](double x, double y) {
        double res{};
        res = 0.75 * exp(-(9 * x - 2) * (9 * x - 2) / 4 - (9 * y - 2) * (9 * y - 2) / 4) + 
            0.75 * exp(-(9 * x + 1) * (9 * x + 1) / 49 - (9 * y + 1) / 10) + 
            0.5 * exp(-(9 * x - 7) * (9 * x - 7) / 4 - (9 * y - 3) * (9 * y - 3) / 4) - 
            0.2 * exp(-(9 * x - 4) * (9 * x - 4) - (9 * y - 7) * (9 * y - 7));
        return res;
    };

    // Grid parameters
    int grid_size = 100;
    double x_min = 0.0, x_max = 1.0;
    double y_min = 0.0, y_max = 1.0;
    double dx = (x_max - x_min) / (grid_size - 1);
    double dy = (y_max - y_min) / (grid_size - 1);

    // Create directories
    std::filesystem::create_directories("./plot/");
    std::filesystem::create_directories("./plot/files/");
    std::filesystem::create_directories("./plot/figures/");

    // Open a file to write the data
    std::ofstream data_file("./plot/files/franke_function_data.dat");
    if (!data_file) {
        std::cerr << "Error opening file to write data.\n";
        return 1;
    }

    // Generate the data
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            double x = x_min + i * dx;
            double y = y_min + j * dy;
            double z = Franke_function(x, y);
            data_file << x << " " << y << " " << z << "\n";
        }
        data_file << "\n";  // New line to separate data blocks in gnuplot
    }

    data_file.close();

    // Create a gnuplot script file
    std::ofstream gnuplot_script("./plot/files/plot_franke_function.gnuplot");
    if (!gnuplot_script) {
        std::cerr << "Error opening file to write gnuplot script.\n";
        return 1;
    }

    gnuplot_script << "set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n";
    gnuplot_script << "set output './plot/figures/franke_function_plot.png'\n";
    gnuplot_script << "set title 'Franke Function'\n";
    gnuplot_script << "set xlabel 'x'\n";
    gnuplot_script << "set ylabel 'y'\n";
    gnuplot_script << "set zlabel 'F(x,y)'\n";
    gnuplot_script << "set pm3d at s\n";
    gnuplot_script << "set palette defined (0 'blue', 1 'green', 2 'yellow', 3 'red')\n";
    gnuplot_script << "set xrange [0:1]\n";
    gnuplot_script << "set yrange [1:0]\n";
    gnuplot_script << "set zrange [0:1.6]\n";
    gnuplot_script << "set view 60, 30, 1, 1\n";
    gnuplot_script << "set xyplane at 0\n";
    gnuplot_script << "splot './plot/files/franke_function_data.dat' with pm3d\n";

    gnuplot_script.close();

    // Run gnuplot with the script
    system("gnuplot ./plot/files/plot_franke_function.gnuplot");


    ///// INTERPOLATION

    // Function that return many sets of parameters.
    auto fillParameters = [](int num_sets, int num_params, double inf, double sup) {
        Eigen::MatrixXd parameters(num_sets, num_params);
        for (size_t i = 0; i < num_sets; ++i) {
            for (size_t j = 0; j < num_params; ++j) {
                parameters(i, j) = generateRandomNumber(inf, sup);
            }
        }
        return parameters;
    };    

    auto funToInterpolate = [](const Eigen::VectorXd& params) {
        double res{};
        int num_dimensions{static_cast<int>(params.size())};
        if (num_dimensions == 2)
        {
            double x{params(0)};
            double y{params(1)};
            res = 0.75 * exp(-(9 * x - 2) * (9 * x - 2) / 4 - (9 * y - 2) * (9 * y - 2) / 4) + 
            0.75 * exp(-(9 * x + 1) * (9 * x + 1) / 49 - (9 * y + 1) / 10) + 
            0.5 * exp(-(9 * x - 7) * (9 * x - 7) / 4 - (9 * y - 3) * (9 * y - 3) / 4) - 
            0.2 * exp(-(9 * x - 4) * (9 * x - 4) - (9 * y - 7) * (9 * y - 7));
        } else if (num_dimensions == 3)
        {
            double x{params(0)};
            double y{params(1)};
            double z{params(2)};
            res = 0.75 * exp(-(9 * x - 2) * (9 * x - 2) / 4 + (9 * y - 2) * (9 * y - 2) / 4 + (9 * z - 2) * (9 * z - 2) / 4) + 
            0.75 * exp(-(9 * x + 1) * (9 * x + 1) / 49 - (9 * y + 1) / 10 - (9 * z + 1) / 10) + 
            0.5 * exp(-(9 * x - 7) * (9 * x - 7) / 4 - (9 * y - 3) * (9 * y - 3) / 4 + (9 * z - 5) * (9 * z - 5) / 4) - 
            0.2 * exp(-(9 * x - 4) * (9 * x - 4) - (9 * y - 7) * (9 * y - 7) - (9 * z - 5) * (9 * z - 5));
        } else {
            std::cerr << "Error: Franke function not defined for this number of dimensions" << std::endl;
        }
        return res;
    };

    auto start = std::chrono::high_resolution_clock::now();

    size_t num_params{2};
    std::string dim{"2D"};

    Eigen::VectorXi KNOWN_POINTS(7);
    KNOWN_POINTS << 100, 500, 1000, 2000, 3000, 4000, 5000;

    int steps{static_cast<int>(KNOWN_POINTS.size())};
   
    size_t num_points{2500}; // number of points to interpolate

    // Mean Squared Error
    Eigen::VectorXd MSE(KNOWN_POINTS.size());
    MSE.setZero();

    // Max Error
    Eigen::VectorXd MaxE(KNOWN_POINTS.size());
    MaxE.setZero();

    // Relative Error 
    // Column : X Y ERROR
    Eigen::MatrixXd RE(num_points, num_params+1);
    RE.setZero();

    Eigen::MatrixXd parameters(KNOWN_POINTS(KNOWN_POINTS.size()-1), num_params);
    Eigen::MatrixXd parametersFORinterp(num_points, num_params);

    double inf{0};
    double sup{1};
    parametersFORinterp = fillParameters(num_points, num_params, inf, sup);
    parameters = fillParameters(KNOWN_POINTS(KNOWN_POINTS.size()-1), num_params, inf, sup);

    RE.block(0,0,num_points, num_params) = parametersFORinterp.block(0,0,num_points, num_params);

    for (size_t j = 0; j < KNOWN_POINTS.size(); ++j)
    {
        std::cout << "\rstep " << (j+1) << "/" << steps;
        std::cout.flush();
        usleep(100000);

        size_t num_measures{static_cast<size_t>(KNOWN_POINTS(j))};

        Eigen::VectorXd measurements(num_measures);
        for (size_t m = 0; m < num_measures; ++m)
        {
            measurements(m) = funToInterpolate((parameters.block(0, 0, num_measures, num_params)).row(m));
        }

        double scale_factor{0};
        RBFInterpolator interpolatorRBF(rbfunForCompact, scale_factor);
        Eigen::VectorXd RBF_points_interpolated = interpolatorRBF.interpolate(parametersFORinterp.block(0, 0, num_points, num_params), parameters.block(0, 0, num_measures, num_params), measurements);

        Eigen::VectorXd points_real(num_points);
        for (size_t k = 0; k < num_points; ++k)
        {
            points_real(k) = funToInterpolate((parametersFORinterp.block(0, 0, num_points, num_params)).row(k));
        }
        
        for (size_t l = 0; l < points_real.size(); ++l)
        {
            MSE(j) += (points_real(l) - RBF_points_interpolated(l))*(points_real(l) - RBF_points_interpolated(l));
            RE(l,num_params) = (RBF_points_interpolated(l) - points_real(l)) / abs(points_real(l));

            double err{abs(RBF_points_interpolated(l) - points_real(l))};
            if (err > MaxE(j))
            {
                MaxE(j) = err;
            }            
        }
        MSE(j) /= num_points;

        plotRelativeError(RE, std::to_string(num_measures), dim);
    }

    std::cout << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
    duration -= hours;
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
    duration -= minutes;
    auto seconds = duration;

    std::cout << "Computation time (only interpolation): " << std::setfill('0') << std::setw(2) << hours.count() << ":"
              << std::setfill('0') << std::setw(2) << minutes.count() << ":"
              << std::setfill('0') << std::setw(2) << seconds.count() << std::endl;

    bool EXPORT{true};
    std::string error{"MSE"};
    plotData(KNOWN_POINTS, MSE, dim, error, EXPORT);

    std::string error2{"Max"};
    plotData(KNOWN_POINTS, MaxE, dim, error2, EXPORT);

    return 0;
}