set(MQDS_IPLDM_LIB
        ipldm_master.f90
        calculate_ipldm_redmat.f90
        )


include_directories(${MQDS_BINARY_DIR}/mqds/src/general_src)

add_library(ipldmlib ${MQDS_IPLDM_LIB})

### If MPI was found, link the general library with MPI
if(MPI_FOUND)
    target_link_libraries(ipldmlib
            generallib
            MPI::MPI_Fortran
            )
else()
    target_link_libraries(ipldmlib
            generallib
            )
endif()
