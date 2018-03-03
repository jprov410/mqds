set(MQDS_TWA_LIB_MPI
        twa_master_mpi.f90
        calculate_twa_redmat_mpi.f90
        )

### Link the twa library with general library
include_directories(${MQDS_BINARY_DIR}/mqds/src/general_src)
include_directories(${MQDS_BINARY_DIR}/mqds/src/mpi_src)

### add a library containing all of the source other than main.f90
add_library(twalib_mpi ${MQDS_TWA_LIB_MPI})

target_link_libraries(twalib_mpi
        generallib
        mpilib
        MPI::MPI_Fortran
        )
