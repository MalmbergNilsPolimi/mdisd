import sys
sys.path.append('../../build')

import mdisd_py.rbfunction as rbfunction
import mdisd_py

import numpy as np
import matplotlib.pyplot as plt



#Test of the different radial basis functions.
print("\n")

print("Multiquadratique(2.5 , 1.0) = ", rbfunction.call_Multiquadratic(2.5, 1.0))
print("invMultiquadratique(2.5 , 1.0) = ", rbfunction.call_invMultiquadratic(2.5, 1.0))
print("Gaussian(2.5, 1.0) = ", rbfunction.call_Gaussian(2.5, 1.0))
print("ThinPlateSpline(2.5, 1.0) = ", rbfunction.call_ThinPlateSpline(2.5, 1.0))

print("\n")

# Test of the different rescaling methods.
data1 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
data2 = np.array([[10, 20, 30], [40, 50, 60], [70, 80, 90]])

mean_normalized_data1, mean_normalized_data2 = mdisd_py.Rescaling().Mean(data1, data2)
print("Mean normalized data1:\n", mean_normalized_data1)
print("Mean normalized data2:\n", mean_normalized_data2)

print("\n")

minmax_normalized_data1, minmax_normalized_data2 = mdisd_py.Rescaling().MinMax(data1, data2)
print("minmax normalized data1:\n", minmax_normalized_data1)
print("minmax normalized data2:\n", minmax_normalized_data2)

print("\n")

zscore_normalized_data1, zscore_normalized_data2 = mdisd_py.Rescaling().Zscore(data1, data2)
print("zscore normalized data1:\n", zscore_normalized_data1)
print("zscore normalized data2:\n", zscore_normalized_data2)


# Test of interpolation.

def f(x):
    return 0.5 * x - 4.3

known_parameters = np.array([[-2], [3.7], [0.1], [-6], [18.2]])
known_measurements = f(known_parameters)

points_to_interpolate = []
inf , sup = -10 , 20
num_points = 10

for i in range(num_points):
    points_to_interpolate.append(inf + i * (sup - inf) / (num_points - 1))


# OLS interpolation.
ols_interpolator = mdisd_py.OLSInterpolator()

ols_results = ols_interpolator.interpolate(points_to_interpolate, known_parameters, known_measurements)
ols_results2, ols_regression = ols_interpolator.interpolate_with_regression(points_to_interpolate, known_parameters, known_measurements)


print("\n")

print("OLS weights:\n", ols_regression)


# RBF interpolation.
rbf_interpolator = mdisd_py.RBFInterpolator(rbfunction.Multiquadratic(), 0, False, False)

rbf_results = rbf_interpolator.interpolate(points_to_interpolate, known_parameters, known_measurements)
rbf_results2, rbf_regression = rbf_interpolator.interpolate_with_regression(points_to_interpolate, known_parameters, known_measurements)


print("RBF weights:\n", rbf_regression)

# NRBF interpolation.
nrbf_interpolator = mdisd_py.RBFInterpolator(rbfunction.Multiquadratic(), 0, True, False)

nrbf_results = nrbf_interpolator.interpolate(points_to_interpolate, known_parameters, known_measurements)
nrbf_results2, nrbf_regression = nrbf_interpolator.interpolate_with_regression(points_to_interpolate, known_parameters, known_measurements)


print("NRBF weights:\n", nrbf_regression)


# RBFP interpolation.
rbfp_interpolator = mdisd_py.RBFInterpolator(rbfunction.Multiquadratic(), 0, False, True)

rbfp_results = rbfp_interpolator.interpolate(points_to_interpolate, known_parameters, known_measurements)
rbfp_results2, rbfp_regression = rbfp_interpolator.interpolate_with_regression(points_to_interpolate, known_parameters, known_measurements)

print("RBFP weights:\n", rbfp_regression)




plt.scatter(known_parameters, known_measurements, label="regressors", color='green')
plt.plot(points_to_interpolate, ols_results, label="OLS interpolation", linestyle='dashed', color='red')
plt.plot(points_to_interpolate, rbf_results, label="RBF interpolation", color='blue')
plt.plot(points_to_interpolate, nrbf_results, label="NRBF interpolation", color='green')
plt.scatter(points_to_interpolate, rbfp_results, label="RBFP interpolation", color='orange')

plt.xlabel("x")
plt.ylabel("f(x)")
plt.legend()
plt.savefig('plot_test_python.png')