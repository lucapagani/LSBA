# Try to find LSMG library
# Once done this will define
#  LSMG_FOUND - if system found LSMG library
#  LSMG_INCLUDE_DIRS - The LSMG include directories
#  LSMG_LIBRARIES - The libraries needed to use LSMG

find_package(PkgConfig)
pkg_check_modules(PC_LSMG QUIET liblsmg)
set(LSMG_DEFINITIONS ${PC_LSMG_CFLAGS_OTHER})

find_path(LSMG_INCLUDE_DIR LSsystem.h
          HINTS ${PC_LSMG_INCLUDEDIR} ${PC_LSMG_INCLUDE_DIRS}
          PATH_SUFFIXES lsmg)

find_library(LSMG_LIBRARY NAMES lsmg liblsmg
             HINTS ${PC_LSMG_LIBDIR} ${PC_LSMG_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LSMG_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LSMG DEFAULT_MSG
                                  LSMG_LIBRARY LSMG_INCLUDE_DIR)

mark_as_advanced(LSMG_INCLUDE_DIR LSMG_LIBRARY)

set(LSMG_LIBRARIES ${LSMG_LIBRARY} )
set(LSMG_INCLUDE_DIRS ${LSMG_INCLUDE_DIR} )
