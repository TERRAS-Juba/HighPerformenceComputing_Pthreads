cmake_minimum_required(VERSION 3.20)
project(Projet_HPC C)

set(CMAKE_C_STANDARD 99)

add_executable(Projet_HPC main.c)
# openMP configuration
FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()
