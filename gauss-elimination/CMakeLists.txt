cmake_minimum_required(VERSION 3.14)
project(gauss_elimination C)

set(CMAKE_C_STANDARD 99)

add_executable(gauss_elimination main.c libs/matrix.c)