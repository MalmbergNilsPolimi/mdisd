#include "data.hpp"

Data::Data(const Eigen::VectorXd& measurements, const Eigen::MatrixXd& parameters)
    : measurements(measurements), parameters(parameters) {
}

const Eigen::VectorXd& Data::getMeasurements() const {
    return measurements;
}

const Eigen::MatrixXd& Data::getParameters() const {
    return parameters;
}