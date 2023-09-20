cmake_minimum_required(VERSION 3.12)
project(my_project)

# Specify the Python version we want to use.
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

# Run Python code to get the file extension for Python extensions on the current platform
execute_process(COMMAND "${Python_EXECUTABLE}" -c 
    "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'), end='')" 
    OUTPUT_VARIABLE PYTHON_EXTENSION_SUFFIX)

# Check if the compiler supports C++20
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++20" COMPILER_SUPPORTS_CXX20)
if(COMPILER_SUPPORTS_CXX20)
    set(CMAKE_CXX_STANDARD 20)
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++20 support. Please use a different C++ compiler.")
endif()

# Detect the compiler details
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message("You're using GCC.")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("You're using Clang.")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    message("You're using MSVC.")
else()
    message("You're using a compiler I don't recognize.")
endif()

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -Wall -fopenmp -fopenmp-simd -funsafe-math-optimizations -DNDEBUG -pipe -shared -fPIC")

add_library(dot MODULE pybind11.cpp)

# Remove the "lib" prefix from the output library
set_target_properties(dot PROPERTIES PREFIX "")

# Set the output library's suffix to the Python extension suffix
set_target_properties(dot PROPERTIES SUFFIX ${PYTHON_EXTENSION_SUFFIX})

# Link with Python libraries
target_include_directories(dot PRIVATE ${Python_INCLUDE_DIRS})
target_link_libraries(dot PRIVATE pybind11::module ${Python_LIBRARIES})
