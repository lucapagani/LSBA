cmake_minimum_required(VERSION 3.0)
project(lsba)

set(LSBA_MAJOR_VERSION 1)
set(LSBA_MINOR_VERSION 0)
set(LSBA_VERSION ${LSBA_MAJOR_VERSION}.${LSBA_MINOR_VERSION})

# Set policy
cmake_policy(SET CMP0012 NEW)

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR cmake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/cmake/lsba)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

# set output paths
set(PROJECT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

set(CMAKE_CONFIGURATION_TYPES "Debug;RelWithDebInfo;Release" CACHE STRING "possible configurations" FORCE)

# In case the user does not setup CMAKE_BUILD_TYPE, assume it's Release
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
set(CMAKE_BUILD_TYPE Release CACHE STRING "build type default to Release" FORCE)
endif("${CMAKE_BUILD_TYPE}" STREQUAL "")

if (CMAKE_BUILD_TYPE STREQUAL "Release")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_NO_DEBUG")
endif ()
if (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_NO_DEBUG")
endif ()
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH})

# # Set shared or static library
# option ( BUILD_SHARED_LIBS "Build libraries shared" ON )

# if(BUILD_SHARED_LIBS)
  # set ( LIBRARY_TYPE SHARED )
# else ()
  # set ( LIBRARY_TYPE STATIC )
# endif ()

# Set include directory
include_directories( ${CMAKE_SOURCE_DIR}/include/ ${CMAKE_BINARY_DIR} )

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package ( Eigen3 3.3 REQUIRED NO_MODULE )

find_package ( Boost REQUIRED )
include_directories ( ${Boost_INCLUDE_DIRS} )

find_package ( MBA REQUIRED )
include_directories ( ${MBA_INCLUDE_DIRS} )

find_package ( LSMG REQUIRED )
include_directories ( ${LSMG_INCLUDE_DIRS} )

find_package ( Cubature REQUIRED )
include_directories ( ${Cubature_INCLUDE_DIRS} )

include(CMakeToolsHelpers OPTIONAL)

file ( GLOB src_list "./src/*.cpp" )
file ( GLOB include_list "./include/*.hpp" )
list ( APPEND include_list "${CMAKE_BINARY_DIR}/lsba_EXPORTS.h" )

# add_library(lsba ${LIBRARY_TYPE} ${src_list})
add_library(lsba STATIC ${src_list})
target_link_libraries(lsba Eigen3::Eigen ${Boost_LIBRARIES} ${MBA_LIBRARIES} ${LSMG_LIBRARIES} ${Cubature_LIBRARIES})

include ( GenerateExportHeader )

if ( WIN32 OR MINGW ) 
  set(CMAKE_C_VISIBILITY_PRESET hidden)
  set(CMAKE_CXX_VISIBILITY_PRESET hidden)
endif ()

generate_export_header ( lsba
  BASE_NAME lsba
  EXPORT_MACRO_NAME LSBA_EXPORT
  EXPORT_FILE_NAME lsba_EXPORTS.h
  STATIC_DEFINE SHARED_EXPORTS_BUILT_AS_STATIC
  )

install(TARGETS lsba 
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  )
install(FILES ${include_list} DESTINATION  include/lsba)

option(BUILD_BIN "Build bin" OFF)

if ( BUILD_BIN )
  # add_subdirectory ( bin/lsba_form_remove )
endif ()

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE lsba)

# Create the LSBAConfig and LSBAConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(LSBAConfig.cmake.in "${PROJECT_BINARY_DIR}/LSBAConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${LSBA_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(LSBAConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/LSBAConfig.cmake" @ONLY)
# ... for both
configure_file(LSBAConfigVersion.cmake.in "${PROJECT_BINARY_DIR}/LSBAConfigVersion.cmake" @ONLY)
# Install the LSBAConfig.cmake and LSBAConfigVersion.cmake
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/LSBAConfig.cmake" "${PROJECT_BINARY_DIR}/LSBAConfigVersion.cmake" DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
# Install FindPackage files
install(FILES "${PROJECT_SOURCE_DIR}/cmake/FindCubature.cmake" "${PROJECT_SOURCE_DIR}/cmake/FindLSMG.cmake" "${PROJECT_SOURCE_DIR}/cmake/FindMBA.cmake" DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

# CPack packaging
include(InstallRequiredSystemLibraries)
# set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/License.txt")
set(CPACK_PACKAGE_VERSION_MAJOR "${LSBA_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${LSBA_VERSION_MINOR}")
include(CPack)
