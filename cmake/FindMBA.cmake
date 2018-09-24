# Try to find MBA library
# Once done this will define
#  MBA_FOUND - if system found MBA library
#  MBA_INCLUDE_DIRS - The MBA include directories
#  MBA_LIBRARIES - The libraries needed to use MBA

find_package(PkgConfig)
pkg_check_modules(PC_MBA QUIET libmba)
set(MBA_DEFINITIONS ${PC_MBA_CFLAGS_OTHER})

find_path(MBA_INCLUDE_DIR MBA.h
          HINTS ${PC_MBA_INCLUDEDIR} ${PC_MBA_INCLUDE_DIRS}
          PATH_SUFFIXES mba)

find_library(MBA_LIBRARY NAMES mba libmba
             HINTS ${PC_MBA_LIBDIR} ${PC_MBA_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MBA_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(MBA DEFAULT_MSG
                                  MBA_LIBRARY MBA_INCLUDE_DIR)

mark_as_advanced(MBA_INCLUDE_DIR MBA_LIBRARY)

set(MBA_LIBRARIES ${MBA_LIBRARY} )
set(MBA_INCLUDE_DIRS ${MBA_INCLUDE_DIR} )
