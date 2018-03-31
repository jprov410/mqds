set(MQDS_PLDM_LIB
        pldm_master.f90
        calculate_pldm_redmat.f90
        calculate_pldm_absorption.f90
        calculate_pldm_nonlinear.f90
        )


include_directories(${MQDS_BINARY_DIR}/mqds/src/general_src)

add_library(pldmlib ${MQDS_PLDM_LIB})

### If MPI was found, link the general library with MPI
if(MPI_FOUND)
    target_link_libraries(pldmlib
            generallib
            MPI::MPI_Fortran
            )
else()
    target_link_libraries(pldmlib
            generallib
            )
endif()
