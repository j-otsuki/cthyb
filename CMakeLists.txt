cmake_minimum_required(VERSION 3.2)
project (hybqmc)

find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

find_package(GSL REQUIRED)

add_subdirectory(src)
