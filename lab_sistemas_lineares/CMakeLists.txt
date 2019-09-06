cmake_minimum_required(VERSION 3.14)
project(lab_sistemas_lineares C)

set(CMAKE_C_STANDARD 99)

add_executable(lab_sistemas_lineares SistemasLineares.c utils.c labSisLin.c main.c)