set(MQDS_SURFACE_HOPPING_LIB
        surface_hopping_master.f90
        )

### Link the sqc library with general library
include_directories(${MQDS_BINARY_DIR}/mqds/src/general_src)

add_library(surfacehoppinglib ${MQDS_SURFACE_HOPPING_LIB})

### If MPI was found, link the sqc library with MPI
if(MPI_FOUND)
    target_link_libraries(surfacehoppinglib
            generallib
            MPI::MPI_Fortran
            )
else()
    target_link_libraries(surfacehoppinglib
            generallib
            )
endif()
