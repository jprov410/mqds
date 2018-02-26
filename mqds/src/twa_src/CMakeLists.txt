set(MQDS_TWA_LIB
        twa_master.f90
        calculate_twa_redmat.f90
        )

### Link the twa library with general library
include_directories(${MQDS_BINARY_DIR}/mqds/src/general_src)

### add a library containing all of the source other than main.f90
add_library(twalib ${MQDS_TWA_LIB})

### If MPI was found, link the general library with MPI
if(MPI_FOUND)
    target_link_libraries(twalib
            generallib
            MPI::MPI_Fortran
            )
else()
    target_link_libraries(twalib
            generallib
            )
endif()
