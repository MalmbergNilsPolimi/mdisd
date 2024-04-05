#ifndef DATA_HPP
#define DATA_HPP

#include <vector>

class Data {
private:
    std::vector<double> measurements;
    std::vector<std::vector<double>> parameters;
public:
    Data(const std::vector<double>& measurements, const std::vector<std::vector<double>>& parameters);
    const std::vector<double>& getMeasurements() const;
    const std::vector<std::vector<double>>& getParameters() const;
};

#endif // DATA_HPP