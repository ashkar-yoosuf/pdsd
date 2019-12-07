cmake_minimum_required(VERSION 3.15)
project(ParAlgo)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP REQUIRED)
add_executable(ParAlgo main.cpp)

target_link_libraries(ParAlgo PRIVATE OpenMP::OpenMP_CXX)