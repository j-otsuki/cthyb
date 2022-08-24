cmake_minimum_required(VERSION 3.2)
project (hybqmc)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

find_package(GSL REQUIRED)  # cmake>=3.2

find_package(Boost REQUIRED)

# Set dafault CMAKE_BUILD_TYPE to Release
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()

add_subdirectory(src)
