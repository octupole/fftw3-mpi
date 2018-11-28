# Policy
cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(SOURCE Fftw3pp.cpp)

add_library(fftw3pp STATIC ${SOURCE})
