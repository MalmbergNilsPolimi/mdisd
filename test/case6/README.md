To use the test case, first modify the CMakeLists.txt file to update the different path. If you don't have pybind, you can directly clone it in the mdisd repository, so you won't have to modify pybind11 paths. If you don't know how to update the path to Python you can refer to mdisd/doc/mdisdReport.pdf in the section 8.2.1.

To clone pybind11 (first go to mdisd/):
git clone https://github.com/pybind/pybind11.git

Then create a build directory in the mdisd repository:
mkdir build
cd build

Then in the build folder use the following commands:
cmake ..
make

Now you can go back in test/case6 and in the shell use the following commands:
python3 test_python.py

To use the library in your project, after using make, you can:
- install it using (in the build folder): sudo make install.
- or copy mdisd_py.cpython-39-x86_64-linux-gnu.so and paste it in your Python project.
- or use the same two first commands as in test_python.py to add the path to the .so file.