#include "data.hpp"

Data::Data(const std::vector<double>& measurements, const std::vector<std::vector<double>>& parameters)
    : measurements(measurements), parameters(parameters) {
}

const std::vector<double>& Data::getMeasurements() const {
    return measurements;
}

const std::vector<std::vector<double>>& Data::getParameters() const {
    return parameters;
}