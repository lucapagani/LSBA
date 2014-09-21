cmake_minimum_required(VERSION 2.6)
project(lsba)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH})

# add_definitions("-std=c++11 -Wall -Wno-sign-compare -g -O0")
add_definitions("-std=c++11 -Wall -Wno-sign-compare -O3 -DEIGEN_NO_DEBUG -fopenmp")

find_package(Eigen REQUIRED)
include_directories(SYSTEM ${EIGEN_INCLUDE_DIRS})
add_definitions(-DEIGEN_USE_NEW_STDVECTOR -DEIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET)

include_directories(/usr/local/include/mba /usr/local/include/lsmg/include)

INCLUDE_DIRECTORIES(
  ${lsba_SOURCE_DIR}/include
  ${mba_SOURCE_DIR}/include
  ${Boost_INCLUDE_DIRS}
  )


# Linked in libraries
SET(DEPLIBS
  mba
  )


# Make the lsba library
FILE(GLOB_RECURSE lsba_SRCS src/*.C src/*.cpp include/*.hpp)
ADD_LIBRARY(lsba ${lsba_SRCS})
