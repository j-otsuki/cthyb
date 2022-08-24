cmake_minimum_required(VERSION 3.2)
project (hybqmc)

find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

find_package(GSL REQUIRED)

# CMAKE_BUILD_TYPE default to Release
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()

add_subdirectory(src)
