set(MQDS_EQUILIBRIUM_LIB
        equilibrium_master.f90
        calculate_equilibrium_site.f90
        )


include_directories(${MQDS_BINARY_DIR}/mqds/src/general_src)

add_library(equilibriumlib ${MQDS_EQUILIBRIUM_LIB})

### If MPI was found, link the general library with MPI
if(MPI_FOUND)
    target_link_libraries(equilibriumlib
            generallib
            MPI::MPI_Fortran
            )
else()
    target_link_libraries(equilibriumlib
            generallib
            )
endif()
