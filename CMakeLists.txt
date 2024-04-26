cmake_minimum_required(VERSION 3.0)

project(mdisd VERSION 1.1.0)

set(SUBDIRECTORIES "OLS" "RBF")

file(GLOB_RECURSE MAIN_SOURCE "src/*.cpp")
set(SOURCES ${MAIN_SOURCE})

# Retrieve all source files in the subdirectories under src.
foreach(SUBDIR ${SUBDIRECTORIES})
    file(GLOB_RECURSE ${SUBDIR}_SOURCES "src/${SUBDIR}/*.cpp")
    set(SOURCES ${SOURCES} ${${SUBDIR}_SOURCES})
endforeach()

# Include the main header directory and subdirectories.
include_directories(include)
foreach(SUBDIR ${SUBDIRECTORIES})
    include_directories(include/${SUBDIR})
endforeach()

# Find Python3 and set Python_INCLUDE_DIRS.
###### The user should update the required version if needed.
find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter Development)

# Set pybind11 directory.
###### The user should update the path to the pybind11 directory if it's located elsewhere.
set(pybind11_DIR "pybind11")

# Find Eigen3 using pkg-config.
find_package(PkgConfig REQUIRED)
pkg_search_module(EIGEN3 REQUIRED eigen3)

# Include directories.
include_directories(${Python_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${pybind11_INCLUDE_DIRS})

# Add pybind11.
add_subdirectory(${pybind11_DIR})

# Create the shared library.
add_library(mdisd SHARED ${SOURCES})

# Link necessary Eigen3 library.
target_link_libraries(mdisd ${EIGEN3_LIBRARIES})

# Create Python module with pybind11.
pybind11_add_module(mdisd_py MODULE src/bindings.cpp)

# Link the library with the Python module.
target_link_libraries(mdisd_py PRIVATE mdisd pybind11::module)

# Set C++ standard.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Install Python module.
###### The user should update the destination path if needed.
install(TARGETS mdisd_py DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.9/site-packages)

# Set library version.
set_target_properties(mdisd PROPERTIES VERSION ${PROJECT_VERSION})

# Set .so file version number.
set_target_properties(mdisd PROPERTIES SOVERSION 1)