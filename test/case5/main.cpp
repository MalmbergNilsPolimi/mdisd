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
T generateRandomNumber(T mean, T std) {
    std::random_device rd;
    std::mt19937 gen(rd());;
    std::normal_distribution<T> dis(mean, std);
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

    std::string title{"Convergence of the RBF "+error+" error for different dimensions"};
    std::string label_y{error+" error on interpolated points"};

    std::ofstream gnuplotScript("./plot/files/plot_script.gnu");

    if (EXPORT)
    {
        // To export the plot in svg format
        gnuplotScript << "set terminal svg" << std::endl;
        gnuplotScript << "set output \"" << "./plot/figures/case5_plot.svg" << "\"" << std::endl;
    }

    gnuplotScript << "set title \"" << title << "\"" << std::endl;
    gnuplotScript << "set xlabel \"number of known points\"" << std::endl;
    gnuplotScript << "set ylabel \"" << label_y << "\"" << std::endl;
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
    auto fillParameters = [](int num_sets, int num_params, double mean, double std) {
        Eigen::MatrixXd parameters(num_sets, num_params);
        for (size_t i = 0; i < num_sets; ++i) {
            for (size_t j = 0; j < num_params; ++j) {
                parameters(i, j) = generateRandomNumber(mean, std);
            }
        }
        return parameters;
    };

    // Definition of the function to interpolate (THE USER CAN CHANGE HERE THE FUNCTION)
    auto funToInterpolate = [](const Eigen::VectorXd& params, const Eigen::MatrixXd& coeff) {
        double res{};

        for (size_t i = 0; i < params.size(); ++i)
        {
            res += params(i);
        }
        return res;
    };

    Eigen::VectorXi DIMENSIONS(6);
    DIMENSIONS << 1, 2, 4, 6, 8, 10;

    Eigen::VectorXi KNOWN_POINTS(14);
    KNOWN_POINTS << 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000, 2000;

    size_t num_points{50};       // number of points to interpolate

    Eigen::MatrixXd MSE(DIMENSIONS.size(), KNOWN_POINTS.size());
    MSE.setZero();
    Eigen::MatrixXd MAD(DIMENSIONS.size(), KNOWN_POINTS.size());
    MAD.setZero();
    Eigen::MatrixXd RMSE(DIMENSIONS.size(), KNOWN_POINTS.size());
    RMSE.setZero();
    Eigen::MatrixXd ME(DIMENSIONS.size(), KNOWN_POINTS.size());
    ME.setZero();

    double meanGauss{};
    double stdGauss{1};

    Eigen::MatrixXd parameters(KNOWN_POINTS(KNOWN_POINTS.size()-1), DIMENSIONS(DIMENSIONS.size()-1));
    Eigen::MatrixXd parametersFORinterp(num_points, DIMENSIONS(DIMENSIONS.size()-1));

    // Creation of the parameters for which the user want the estimated output
    parametersFORinterp = fillParameters(num_points, DIMENSIONS(DIMENSIONS.size()-1), meanGauss, stdGauss);
    parameters = fillParameters(KNOWN_POINTS(KNOWN_POINTS.size()-1), DIMENSIONS(DIMENSIONS.size()-1), meanGauss, stdGauss);

    Eigen::MatrixXd coeff_fun = Eigen::MatrixXd::Random(2, DIMENSIONS(DIMENSIONS.size()-1));
    coeff_fun = (coeff_fun.array() + 1.0) * 5.0;


    for (size_t i = 0; i < DIMENSIONS.size(); ++i)
    {
        size_t num_params{static_cast<size_t>(DIMENSIONS(i))};

        for (size_t j = 0; j < KNOWN_POINTS.size(); ++j)
        {
            size_t num_measures{static_cast<size_t>(KNOWN_POINTS(j))};

            Eigen::VectorXd measurements(num_measures);
            for (size_t m = 0; m < num_measures; ++m)
            {
                measurements(m) = funToInterpolate((parameters.block(0, 0, num_measures, num_params)).row(m), coeff_fun);
            }

            double scale_factor{sqrt(num_params)/10};
            RBFInterpolator interpolatorRBF(&RBFunctions::gaussian, scale_factor);

            Eigen::VectorXd RBF_points_interpolated = interpolatorRBF.interpolate(parametersFORinterp.block(0, 0, num_points, num_params), parameters.block(0, 0, num_measures, num_params), measurements);


            Eigen::VectorXd tmp = RBF_points_interpolated;

            Eigen::VectorXd points_real(num_points);
            for (size_t k = 0; k < num_points; ++k)
            {
                points_real(k) = funToInterpolate((parametersFORinterp.block(0, 0, num_points, num_params)).row(k), coeff_fun);
            }

            for (size_t l = 0; l < points_real.size(); ++l)
            {
                MSE(i,j) += (points_real(l) - RBF_points_interpolated(l))*(points_real(l) - RBF_points_interpolated(l));
                MAD(i,j) += abs(points_real(l) - RBF_points_interpolated(l));
                ME(i,j) += points_real(l) - RBF_points_interpolated(l);
            }
            MSE(i,j) /= num_points;
            MAD(i,j) /= num_points;
            RMSE(i,j) = sqrt(MSE(i,j));
            ME(i,j) /= num_points;
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
    plotData(DIMENSIONS, KNOWN_POINTS, MAD, error, EXPORT);

    error = "RMSE";
    plotData(DIMENSIONS, KNOWN_POINTS, RMSE, error, EXPORT);

    error = "MEAN";
    plotData(DIMENSIONS, KNOWN_POINTS, ME, error, EXPORT);
    
    return 0;
}