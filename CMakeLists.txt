cmake_minimum_required(VERSION 3.0) # version of cmake to use.

project(mdisd VERSION 1.1.0) # Name and version of the library.

set(SUBDIRECTORIES "OLS" "RBF") # List of subdirectories to include.

file(GLOB_RECURSE MAIN_SOURCE "src/*.cpp")
list(APPEND SOURCES ${MAIN_SOURCE})

# Retrieve all source files in the subdirectories under src.
foreach(SUBDIR ${SUBDIRECTORIES})
    file(GLOB_RECURSE ${SUBDIR}_SOURCES "src/${SUBDIR}/*.cpp")
    list(APPEND SOURCES ${${SUBDIR}_SOURCES})
endforeach()

# Include the main header directory and subdirectories.
include_directories(include)
foreach(SUBDIR ${SUBDIRECTORIES})
    include_directories(include/${SUBDIR})
endforeach()

######################################################################################
######## THE USER NEED TO CHANGE HERE THE PATH TO PYTHON AND PYBIND11 LIBRARY ########
######################################################################################
include_directories(/u/sw/toolchains/gcc-glibc/11.2.0/base/include/python3.9)
include_directories(pybind11/include)
######################################################################################


######################################################################################
################### THE USER NEED TO CHECK THAT EIGEN3 IS INSTALLED ##################
######################################################################################
# Use pkg-config to find information about Eigen3.
find_package(PkgConfig REQUIRED)
pkg_search_module(EIGEN3 REQUIRED eigen3)
######################################################################################

# Include paths to Eigen3 header files.
include_directories(${EIGEN3_INCLUDE_DIRS})

######################################################################################
############### THE USER NEED TO CHANGE HERE THE PATH PYBIND11 LIBRARY ###############
######################################################################################
# Add pybind11.
add_subdirectory(pybind11)
######################################################################################

# Create the shared library.
add_library(mdisd SHARED ${SOURCES})

# Link necessary Eigen3 library.
target_link_libraries(mdisd ${EIGEN3_LIBRARIES})

# Create Python module with pybind11.
pybind11_add_module(mdisd_py MODULE src/bindings.cpp)

# Link the library with the Python module.
target_link_libraries(mdisd_py PRIVATE mdisd)

# Specify C++.
target_compile_features(mdisd_py PRIVATE cxx_std_11)

# Install Python module.
install(TARGETS mdisd_py DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.8/site-packages)

# Set library version.
set_target_properties(mdisd PROPERTIES VERSION ${PROJECT_VERSION})

# Set .so file version number.
set_target_properties(mdisd PROPERTIES SOVERSION 1) # To change with the version