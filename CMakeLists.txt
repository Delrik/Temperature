cmake_minimum_required(VERSION 3.17)
project(Temperature)

set(CMAKE_CXX_STANDARD 14)

add_executable(Temperature main.cpp fGNV.h Integral.h)