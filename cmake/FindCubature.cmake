# Try to find Cubature library
# Once done this will define
#  Cubature_FOUND - if system found Cubature library
#  Cubature_INCLUDE_DIRS - The Cubature include directories
#  Cubature_LIBRARIES - The libraries needed to use Cubature

find_package(PkgConfig)
pkg_check_modules(PC_Cubature QUIET libcubature)
set(Cubature_DEFINITIONS ${PC_Cubature_CFLAGS_OTHER})

find_path(Cubature_INCLUDE_DIR cubature.h
          HINTS ${PC_Cubature_INCLUDEDIR} ${PC_Cubature_INCLUDE_DIRS})

find_library(Cubature_LIBRARY NAMES cubature libcubature
             HINTS ${PC_Cubature_LIBDIR} ${PC_Cubature_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set Cubature_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Cubature DEFAULT_MSG
                                  Cubature_LIBRARY Cubature_INCLUDE_DIR)

mark_as_advanced(Cubature_INCLUDE_DIR Cubature_LIBRARY)

set(Cubature_LIBRARIES ${Cubature_LIBRARY} )
set(Cubature_INCLUDE_DIRS ${Cubature_INCLUDE_DIR} )

