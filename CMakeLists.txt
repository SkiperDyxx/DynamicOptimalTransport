cmake_minimum_required(VERSION 3.12)
project(dot)

find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

execute_process(COMMAND "${Python_EXECUTABLE}" -c
    "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'), end='')"
    OUTPUT_VARIABLE PYTHON_EXTENSION_SUFFIX)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++20" COMPILER_SUPPORTS_CXX20)
if(COMPILER_SUPPORTS_CXX20)
    set(CMAKE_CXX_STANDARD 20)
else()
    message(FATAL_ERROR "No C++20 support.")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message("Using GCC")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("Using Clang")
else()
	message(FATAL_ERROR "Unknown compiler. Please use GCC or Clang")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -Wall -fopenmp -fopenmp-simd -funsafe-math-optimizations -DNDEBUG -pipe -shared -fPIC")

# Update paths for source, include, and third-party directories
add_library(dot MODULE src/pybind11.cpp)
target_include_directories(dot PRIVATE include)
target_include_directories(dot PRIVATE third_party)

set_target_properties(dot PROPERTIES PREFIX "")
set_target_properties(dot PROPERTIES SUFFIX ${PYTHON_EXTENSION_SUFFIX})
target_include_directories(dot PRIVATE ${Python_INCLUDE_DIRS})
target_link_libraries(dot PRIVATE pybind11::module ${Python_LIBRARIES})

