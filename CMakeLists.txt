cmake_minimum_required(VERSION 3.0) # version of cmake to use.

project(mdsid VERSION 1.0.2) # Name and version of the library.

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

# Use pkg-config to find information about Eigen3.
find_package(PkgConfig REQUIRED)
pkg_search_module(EIGEN3 REQUIRED eigen3)

# Include paths to Eigen3 header files.
include_directories(${EIGEN3_INCLUDE_DIRS})

# Create the shared library.
add_library(mdsid SHARED ${SOURCES})

# Link necessary Eigen3 library.
target_link_libraries(mdsid ${EIGEN3_LIBRARIES})

# Set library version.
set_target_properties(mdsid PROPERTIES VERSION ${PROJECT_VERSION})

# Set .so file version number.
set_target_properties(mdsid PROPERTIES SOVERSION 1) # To change with the version