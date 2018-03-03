set(MQDS_PLDM_LIB_MPI
        pldm_master_mpi.f90
        calculate_pldm_redmat_mpi.f90
        calculate_pldm_absorption_mpi.f90
        )


include_directories(${MQDS_BINARY_DIR}/mqds/src/general_src)
include_directories(${MQDS_BINARY_DIR}/mqds/src/mpi_src)

add_library(pldmlib_mpi ${MQDS_PLDM_LIB_MPI})

### If MPI was found, link the general library with MPI
target_link_libraries(pldmlib_mpi
        mpilib
        generallib
        MPI::MPI_Fortran
        )
