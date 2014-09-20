cmake_minimum_required(VERSION 2.6)
project(lsba)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# add_definitions("-std=c++11 -Wall -Wno-sign-compare -g -O0 -pg")
add_definitions("-std=c++11 -Wall -Wno-sign-compare -O3 -DEIGEN_NO_DEBUG -pg")

include_directories(/usr/local/include/mba/include /usr/local/include/lsmg/include)

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
