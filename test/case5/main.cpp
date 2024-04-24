#include <iostream>
#include <cmath>
#include <iomanip>
#include <random>
#include <fstream>
#include <filesystem>
#include <chrono>

#include "RBFunctions.hpp"
#include "RBFinterpolator.hpp"

template<typename T>
T generateRandomNumber(T min, T max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dis(min, max);
    return dis(gen);
}


void plotData(const Eigen::VectorXi& DIMENSIONS, const Eigen::VectorXi& KNOWN_POINTS, const Eigen::MatrixXd& MSE, const std::string error, const bool EXPORT) {

    std::filesystem::create_directories("./plot/");
    std::filesystem::create_directories("./plot/files/");
    std::filesystem::create_directories("./plot/figures/");


    for (size_t i = 0; i < DIMENSIONS.size(); ++i)
    {
        std::ofstream dataFile("./plot/files/MSE_RBF_"+std::to_string(DIMENSIONS(i))+"D.dat");
        if (!dataFile.is_open()) {
            std::cerr << "Error: Unable to open data file." << std::endl;
            return;
        }

        for (int j = 0; j < KNOWN_POINTS.size(); ++j) {
            dataFile << KNOWN_POINTS(j) << " " << MSE(i,j) << std::endl;
        }
        dataFile.close();
    }

    std::string title{"Convergence of the RBF "+error+" for different dimensions"};

    std::ofstream gnuplotScript("./plot/files/plot_script.gnu");

    if (EXPORT)
    {
        // To export the plot in svg format
        gnuplotScript << "set terminal svg" << std::endl;
        gnuplotScript << "set output \"" << "./plot/figures/case5_plot.svg" << "\"" << std::endl;
    }

    gnuplotScript << "set title \"" << title << "\"" << std::endl;
    gnuplotScript << "set xlabel \"number of known points\"" << std::endl;
    gnuplotScript << "set ylabel \"Mean Squared Error on the interpolated points\"" << std::endl;
    gnuplotScript << "set key box" << std::endl;
    
    gnuplotScript << "plot ";
    for (size_t i = 0; i < DIMENSIONS.size(); ++i)
    {
        if (i > 0)
        {
            gnuplotScript << ", ";
        }
        gnuplotScript << " './plot/files/MSE_RBF_" << std::to_string(DIMENSIONS(i)) << "D.dat' with lines lw 2 title \"" << DIMENSIONS(i) << "D\"";
        
    }
    gnuplotScript << std::endl;
    
    if (EXPORT)
    {   
        gnuplotScript << "set terminal wxt" << std::endl;
        gnuplotScript << "replot" << std::endl;
    }

    gnuplotScript.close();
    system("gnuplot -persist ./plot/files/plot_script.gnu");

}



int main() {

    auto start = std::chrono::high_resolution_clock::now();

    ////////////////////////////////////////////////
    //////// TEST : Using random sampling //////////
    ////////////////////////////////////////////////

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

    // Definition of the function to interpolate (THE USER CAN CHANGE HERE THE FUNCTION)
    auto funToInterpolate = [](const Eigen::VectorXd& params, const Eigen::MatrixXd& coeff) {
        double res{1.};
        double inf{-10.};
        double sup{10.};

        for (size_t i = 0; i < params.size()-1; ++i)
        {
            //res += params(i) * params(i) * params(i) * exp(-params(i+1)/2);
            res += coeff(0,i)*params(i)*cos(coeff(1,i)*M_PI*params(i+1));
        }
        return res;
    };


    // Eigen::VectorXi DIMENSIONS(3);
    // DIMENSIONS << 2, 4, 6;

    Eigen::VectorXi DIMENSIONS(6);
    DIMENSIONS << 1, 2, 4, 6, 8, 10;

    Eigen::VectorXi KNOWN_POINTS(9);
    KNOWN_POINTS << 10, 50, 100, 150, 200, 250, 300, 350, 400;

    // Eigen::VectorXi KNOWN_POINTS(13);
    // KNOWN_POINTS << 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000;

    size_t num_points{50};       // number of points to interpolate

    Eigen::MatrixXd MSE(DIMENSIONS.size(), KNOWN_POINTS.size());
    MSE.setZero();
    Eigen::MatrixXd MSE2(DIMENSIONS.size(), KNOWN_POINTS.size());
    MSE2.setZero();
    Eigen::MatrixXd MSE3(DIMENSIONS.size(), KNOWN_POINTS.size());
    MSE3.setZero();
    Eigen::MatrixXd MSE4(DIMENSIONS.size(), KNOWN_POINTS.size());
    MSE4.setZero();

    double scale_factor{0.5};
    RBFInterpolator interpolatorRBF(&RBFunctions::gaussian, scale_factor);


    Eigen::MatrixXd parameters(KNOWN_POINTS(KNOWN_POINTS.size()-1), DIMENSIONS(DIMENSIONS.size()-1));
    Eigen::MatrixXd parametersFORinterp(num_points, DIMENSIONS(DIMENSIONS.size()-1));

    double inf{0.};
    double sup{1.};
    // Creation of the parameters for which the user want the estimated output
    parametersFORinterp = fillParameters(num_points, DIMENSIONS(DIMENSIONS.size()-1), inf, sup);

    Eigen::MatrixXd coeff_fun = Eigen::MatrixXd::Random(2, DIMENSIONS(DIMENSIONS.size()-1) - 1);


    for (size_t i = 0; i < DIMENSIONS.size(); ++i)
    {
        size_t num_params{static_cast<size_t>(DIMENSIONS(i))};

        for (size_t j = 0; j < KNOWN_POINTS.size(); ++j)
        {
            size_t num_measures{static_cast<size_t>(KNOWN_POINTS(j))};

            Eigen::VectorXd measurements(num_measures);
            for (size_t i = 0; i < num_measures; ++i)
            {
                measurements(i) = funToInterpolate((parameters.block(0, 0, num_measures, num_params)).row(i), coeff_fun);
            }

            Eigen::VectorXd RBF_points_interpolated = interpolatorRBF.interpolate(parametersFORinterp.block(0, 0, num_points, num_params), parameters.block(0, 0, num_measures, num_params), measurements);

            Eigen::VectorXd points_real(num_points);
            for (size_t k = 0; k < num_points; ++k)
            {
                points_real(k) = funToInterpolate((parametersFORinterp.block(0, 0, num_points, num_params)).row(k), coeff_fun);
            }

            for (size_t l = 0; l < points_real.size(); ++l)
            {
                MSE(i,j) += (points_real(l) - RBF_points_interpolated(l))*(points_real(l) - RBF_points_interpolated(l)) / num_points;
                MSE2(i,j) += abs(points_real(l) - RBF_points_interpolated(l)) / num_points;
                MSE4(i,j) += (points_real(l) - RBF_points_interpolated(l)) / num_points;
            }
            MSE3(i,j) = sqrt(MSE(i,j));
            
            //MSE(i,j) = (points_real - RBF_points_interpolated).array().pow(2).sum() / num_points;
            //MSE(i,j) = ((RBF_points_interpolated - points_real).array() / points_real.array().abs()).mean();
        }
        
    }



    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
    duration -= hours;
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
    duration -= minutes;
    auto seconds = duration;

    std::cout << "Computation time: " << std::setfill('0') << std::setw(2) << hours.count() << ":"
              << std::setfill('0') << std::setw(2) << minutes.count() << ":"
              << std::setfill('0') << std::setw(2) << seconds.count() << std::endl;


    bool EXPORT{true};
    std::string error{"MSE"};
    plotData(DIMENSIONS, KNOWN_POINTS, MSE, error, EXPORT);

    error = "MAD";
    plotData(DIMENSIONS, KNOWN_POINTS, MSE2, error, EXPORT);

    error = "RMSE";
    plotData(DIMENSIONS, KNOWN_POINTS, MSE3, error, EXPORT);

    error = "MEAN";
    plotData(DIMENSIONS, KNOWN_POINTS, MSE4, error, EXPORT);
    
    return 0;
}