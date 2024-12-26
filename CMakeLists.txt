cmake_minimum_required(VERSION 3.20)
project(fluid)

set(CMAKE_CXX_STANDARD 20)

add_executable(fluid Fixed.h main.cpp Procesing.h ProcesingDynamic.h RandomGenerator.h)

