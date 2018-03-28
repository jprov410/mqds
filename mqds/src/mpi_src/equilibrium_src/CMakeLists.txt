set(MQDS_EQUILIBRIUM_LIB_MPI
        equilibrium_master_mpi.f90
        calculate_equilibrium_site_mpi.f90
        )


include_directories(${MQDS_BINARY_DIR}/mqds/src/general_src)
include_directories(${MQDS_BINARY_DIR}/mqds/src/mpi_src)

add_library(equilibriumlib_mpi ${MQDS_EQUILIBRIUM_LIB_MPI})

### If MPI was found, link the general library with MPI
target_link_libraries(equilibriumlib_mpi
        mpilib
        generallib
        MPI::MPI_Fortran
        )