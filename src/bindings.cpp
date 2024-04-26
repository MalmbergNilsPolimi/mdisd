#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "interpolator.hpp"
#include "OLSinterpolator.hpp"
#include "RBFinterpolator.hpp"
#include "RBFunctions.hpp"
#include "rescaling.hpp"


namespace py = pybind11;


template <typename T>
py::tuple interpolate_with_regression(T& self, const Eigen::MatrixXd& points_to_interpolate, const Eigen::MatrixXd& known_parameters, const Eigen::VectorXd& known_measurements) {
    Eigen::VectorXd regression;
    Eigen::VectorXd results = self.interpolate(points_to_interpolate, known_parameters, known_measurements, &regression);
    return py::make_tuple(results, regression);
}


PYBIND11_MODULE(mdisd_py, m) {

    py::class_<OLSInterpolator>(m, "OLSInterpolator")
        .def(py::init<>())
        .def("interpolate", static_cast<Eigen::VectorXd (OLSInterpolator::*)(const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::VectorXd&) const>(&OLSInterpolator::interpolate), py::arg("points_to_interpolate"), py::arg("known_parameters"), py::arg("known_measurements"))
        .def("interpolate_with_regression", &interpolate_with_regression<OLSInterpolator>);

    py::class_<RBFInterpolator>(m, "RBFInterpolator")
        .def(py::init<std::function<double(double, double)>, double, bool, bool>(), py::arg("rbf_function"), py::arg("scale_factor"), py::arg("flag1"), py::arg("flag2"))
        .def("interpolate", static_cast<Eigen::VectorXd (RBFInterpolator::*)(const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::VectorXd&) const>(&RBFInterpolator::interpolate), py::arg("points_to_interpolate"), py::arg("known_parameters"), py::arg("known_measurements"))
        .def("interpolate_with_regression", &interpolate_with_regression<RBFInterpolator>);

    
    py::class_<std::function<double(double, double)>>(m, "FunctionDD");

    py::module m_rbfunction = m.def_submodule("rbfunction", "RBFunctions submodule");

    auto multiquadratic_wrapper = []() -> std::function<double(double, double)> {
        return [](double r, double r0) -> double {
            return RBFunctions::multiquadratic(r, r0);
        };
    };
    m_rbfunction.def("Multiquadratic", multiquadratic_wrapper, "Multiquadratic radial basis function");

    m_rbfunction.def("call_Multiquadratic", [](double r, double r0) {
        return RBFunctions::multiquadratic(r, r0);
    }, "Call Multiquadratic radial basis function");


    auto invmultiquadratic_wrapper = []() -> std::function<double(double, double)> {
        return [](double r, double r0) -> double {
            return RBFunctions::inverseMultiquadratic(r, r0);
        };
    };
    m_rbfunction.def("invMultiquadratic", invmultiquadratic_wrapper, "Inverse multiquadratic radial basis function");

    m_rbfunction.def("call_invMultiquadratic", [](double r, double r0) {
        return RBFunctions::inverseMultiquadratic(r, r0);
    }, "Call Inverse multiquadratic radial basis function");


    auto gaussian_wrapper = []() -> std::function<double(double, double)> {
        return [](double r, double r0) -> double {
            return RBFunctions::gaussian(r, r0);
        };
    };
    m_rbfunction.def("Gaussian", gaussian_wrapper, "Gaussian radial basis function");

    m_rbfunction.def("call_Gaussian", [](double r, double r0) {
        return RBFunctions::gaussian(r, r0);
    }, "Call Gaussian radial basis function");


    auto thinPlateSpline_wrapper = []() -> std::function<double(double, double)> {
        return [](double r, double r0) -> double {
            return RBFunctions::thinPlateSpline(r, r0);
        };
    };
    m_rbfunction.def("ThinPlateSpline", thinPlateSpline_wrapper, "Thin Plate Spline radial basis function");

    m_rbfunction.def("call_ThinPlateSpline", [](double r, double r0) {
        return RBFunctions::thinPlateSpline(r, r0);
    }, "Call Thin Plate Spline radial basis function");


    py::class_<Rescaling>(m, "Rescaling")
        .def(py::init<>())
        .def("Mean", &Rescaling::meanNormalization, py::arg("data1"), py::arg("data2") = nullptr)
        .def("MinMax", &Rescaling::minMaxNormalization, py::arg("data1"), py::arg("data2") = nullptr)
        .def("Zscore", &Rescaling::zScoreNormalization, py::arg("data1"), py::arg("data2") = nullptr);
}
