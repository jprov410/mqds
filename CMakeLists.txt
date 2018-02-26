cmake_minimum_required(VERSION 3.2)
project(MQDS
	LANGUAGES Fortran
	VERSION 0.0.1
)


include(GNUInstallDirs)

find_package(LAPACK)
find_package(MPI)

if(MPI_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
endif()

message(STATUS "Fortran Flags: ${CMAKE_Fortran_FLAGS}")

add_subdirectory(mqds)

enable_testing()
add_subdirectory(test)
