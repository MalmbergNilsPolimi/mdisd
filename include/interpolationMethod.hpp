#ifndef INTERPOLATIONMETHOD_HPP
#define INTERPOLATIONMETHOD_HPP

#include <vector>

class InterpolationMethod {
public:
    virtual double interpolate(const std::vector<double>& parametersFORinterp, const std::vector<std::vector<double>>& parameters, const std::vector<double>& measurements) const = 0;
};

#endif // INTERPOLATIONMETHOD_HPP
