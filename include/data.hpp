#ifndef DATA_HPP
#define DATA_HPP

#include <Eigen/Core>

class Data {
private:
    Eigen::VectorXd measurements;
    Eigen::MatrixXd parameters;
public:
    Data(const Eigen::VectorXd& measurements, const Eigen::MatrixXd& parameters);
    const Eigen::VectorXd& getMeasurements() const;
    const Eigen::MatrixXd& getParameters() const;
};

#endif // DATA_HPP