set(MQDS_REDFIELD_LIB_MPI
        redfield_master_mpi.f90
        calculate_redfield_rates_mpi.f90
        )


include_directories(${MQDS_BINARY_DIR}/mqds/src/general_src)
include_directories(${MQDS_BINARY_DIR}/mqds/src/mpi_src)

add_library(redfieldlib_mpi ${MQDS_REDFIELD_LIB_MPI})

### If MPI was found, link the general library with MPI
target_link_libraries(redfieldlib_mpi
        mpilib
        generallib
        MPI::MPI_Fortran
        )