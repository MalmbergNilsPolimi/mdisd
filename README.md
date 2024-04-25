<!---
Copyright 2024 Nils Malmberg.
All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
-->

<h1 align="center">
    <p>mdisd</p>
    <p>Multidimensional Interpolation for Scattered Data Library</p>
</h1>

<p align="center">
    <a href="https://github.com/MalmbergNilsPolimi/mdisd/blob/main/LICENSE">
        <img alt="GitHub" src="https://img.shields.io/github/license/MalmbergNilsPolimi/mdisd?color=blue">
    </a>
    <a href="https://github.com/MalmbergNilsPolimi/mdisd/releases">
        <img alt="GitHub release" src="https://img.shields.io/github/release/MalmbergNilsPolimi/mdisd.svg">
    </a>
</p>



The mdisd C++ library is dedicated to the interpolation of scattered data with one or more dimensions. The library is also available for use with Python (IN COMING). It was created in 2024 by Nils Malmberg, then a student at  <a href="https://www.polimi.it/en" target="_blank">Politecnico di Milano (polimi)</a>. This library is a project developed as part of the "Advanced Programming for Scientific Computing" course given at polimi.

## Features

- Several interpolation methods are available: Radial Basis Function (RBF) and Ordinary Least Squares (OLS).
- Interpolation of multi-dimensional data.
- Offers the possibility to normalise RBF interpolation.
- Offers the possibility to add a linear polynomial to RBF interpolation.
- Offers pre-processing options (data rescaling and normalization).
- Efficient computation using <a href="https://eigen.tuxfamily.org/" target="_blank">Eigen</a> library for linear algebra.
- Use of the library in Python possible thanks to bindings made with <a href="https://github.com/pybind/pybind11" target="_blank">pybind11<\a> (IN COMING).

## How to test the library

The library presents 6 test cases numbered from 0 to 5 which can be launched by executing the following commands (example for the test case number 0):

```bash
git clone git@github.com:MalmbergNilsPolimi/mdisd.git
cd mdisd
cd test/case0
make
make run
```

For more information on the various test cases implemented, you can consult the commented code directly or consult the section dedicated to tests in the mdisdReport.pdf report located in the doc directory.

## Installation

Before building the library, make sure you have the following dependencies installed:
- CMake (version 3.0 or higher)
- Eigen3

To build and install the library:
```bash
mkdir build
cd build
cmake ..
make
sudo make install
```

## Usage

Once installed, you can use the library in your C++ projects by including the appropriate headers and linking against the `mdisd` library. Here's an example:

```cpp
#include "mdisd/RBFinterpolator.hpp"

int main() {
    // Your code here
    return 0;
}
```

## Directory Structure

The directory structure of the library is as follows:

- `include/`: Contains all the header files for the library.
- `src/`: Contains the source code for the library.
- `test/`: Contains example programs demonstrating the usage of the library.
- `doc/`: Contains Doxyfile of the library and a report written in 2024.

## Documentation

For detailed documentation on how to use the library, please refer to the comments in the code and the provided examples.

## Contributing

Contributions are welcome! If you find any issues or have suggestions for improvements, please open an issue or submit a pull request on GitHub.

## License

This library is licensed under the Apache License, Version 2.0. See the [LICENSE](LICENSE) file for details.
