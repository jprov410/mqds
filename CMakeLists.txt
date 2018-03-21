cmake_minimum_required(VERSION 3.2)
project(MQDS
	LANGUAGES Fortran
	VERSION 0.0.1
)


include(GNUInstallDirs)

find_package(LAPACK)
find_package(MPI)
find_package(Doxygen)
find_package(Latex)

if(DOXYGEN_FOUND)
    # Set input and output Doxygen files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/dox/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/dox/Doxyfile)

    #make sure this matches the version number above (for documentation)
    set(VERSION_MAJOR 0 CACHE STRING "Project major version number.")
    set(VERSION_MINOR 0 CACHE STRING "Project minor version number.")
    set(VERSION_PATCH 1 CACHE STRING "Project patch version number.")
    mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    # Build the docs with the application
    add_custom_target( docs
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
else()
    message(STATUS "Doxygen needs to be installed to generate the doxygen documentation")
endif(DOXYGEN_FOUND)


if(MPI_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
endif()

message(STATUS "Fortran Flags: ${CMAKE_Fortran_FLAGS}")

add_subdirectory(mqds)

enable_testing()
add_subdirectory(test)
