cmake_minimum_required(VERSION 3.2)
project(MQDS
	LANGUAGES Fortran
	VERSION 1.0
)


include(GNUInstallDirs)

find_package(LAPACK)
find_package(MPI)

message(STATUS "Fortran Flags: ${CMAKE_Fortran_FLAGS}")

add_subdirectory(src)
