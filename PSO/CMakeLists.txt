cmake_minimum_required(VERSION 3.17)
project(untitled1)
include_directories("C:/Users/ASUS/Desktop/eigen3")
set(CMAKE_CXX_STANDARD 14)
find_package(Eigen3 3.3.4 REQUIRED)

add_executable(untitled1 main.cpp)