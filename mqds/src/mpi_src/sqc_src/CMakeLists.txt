set(MQDS_SQC_LIB_MPI
        sqc_master_mpi.f90
        calculate_sqc_redmat_mpi.f90
        calculate_sqc_absorption_mpi.f90
        )

### Link the sqc library with general library
include_directories(${MQDS_BINARY_DIR}/mqds/src/general_src)
include_directories(${MQDS_BINARY_DIR}/mqds/src/mpi_src)

add_library(sqclib_mpi ${MQDS_SQC_LIB_MPI})

### If MPI was found, link the sqc library with MPI
target_link_libraries(sqclib_mpi
        generallib
        mpilib
        MPI::MPI_Fortran
        )
