set(MQDS_SQC_LIB
        sqc_master.f90
        calculate_sqc_redmat.f90
        calculate_sqc_absorption.f90
        )

### Link the sqc library with general library
include_directories(${MQDS_BINARY_DIR}/mqds/src/general_src)

add_library(sqclib ${MQDS_SQC_LIB})

### If MPI was found, link the sqc library with MPI
if(MPI_FOUND)
    target_link_libraries(sqclib
            generallib
            MPI::MPI_Fortran
            )
else()
    target_link_libraries(sqclib
            generallib
            )
endif()
