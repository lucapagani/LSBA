cmake_minimum_required(VERSION 2.6)
project(example_lsba)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH})

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "possible configurations" FORCE)

# In case the user does not setup CMAKE_BUILD_TYPE, assume it's RelWithDebInfo
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
set(CMAKE_BUILD_TYPE Release CACHE STRING "build type default to Release" FORCE)
endif("${CMAKE_BUILD_TYPE}" STREQUAL "")


IF(CMAKE_BUILD_TYPE MATCHES "Release")
  add_definitions("-std=c++11 -Wall -Wno-sign-compare -O3")
ELSE(CMAKE_BUILD_TYPE MATCHES "Release")
  add_definitions("-std=c++11 -Wall -Wno-sign-compare -g -O0")
ENDIF(CMAKE_BUILD_TYPE MATCHES "Release")

find_package(Eigen REQUIRED)
include_directories(SYSTEM ${EIGEN_INCLUDE_DIRS})
add_definitions("-DEIGEN_NO_DEBUG")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(/usr/local/include/mba/include /usr/local/include/lsmg/include /usr/local/include/mba /usr/local/include/lsmg /usr/local/include/lsba)

add_executable(example_lsba src/main.cpp)

target_link_libraries(example_lsba liblsba.a liblsmg.a libmba.a libcubature.a)