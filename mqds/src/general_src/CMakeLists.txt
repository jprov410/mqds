set(MQDS_GENERAL_LIB
        input_output.f90
        unit_conversions.f90
        random_numbers.f90
        harmonic_bath.f90
        mapping_variables.f90
        hamiltonians.f90
        kinds.f90
        parameters.f90
        spectroscopy.f90
        focusing.f90
        windows.f90
        linear_algebra.f90
        imaginary_time.f90)

### add a library containing all of the source other than main.f90
add_library(generallib
        ${MQDS_GENERAL_LIB}
        )

### If MPI was found, link the general library with MPI
if(MPI_Fortran_FOUND)
    target_link_libraries(generallib
            MPI::MPI_Fortran
            ${LAPACK_LIBRARIES}
            ${BLAS_LIBRARIES}
            )
else()
    target_link_libraries(generallib
            ${LAPACK_LIBRARIES}
            ${BLAS_LIBRARIES}
            )
endif()
