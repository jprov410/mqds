set(MQDS_IPLDM_LIB_MPI
        ipldm_master_mpi.f90
        calculate_ipldm_redmat_mpi.f90
        )


include_directories(${MQDS_BINARY_DIR}/mqds/src/general_src)
include_directories(${MQDS_BINARY_DIR}/mqds/src/mpi_src)

add_library(ipldmlib_mpi ${MQDS_IPLDM_LIB_MPI})

### If MPI was found, link the general library with MPI
target_link_libraries(ipldmlib_mpi
        mpilib
        generallib
        MPI::MPI_Fortran
        )

