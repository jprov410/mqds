set(MQDS_GENERAL_LIB
        input_output.f90
        unit_conversions.f90
        random_numbers.f90
        harmonic_bath.f90
        mapping_variables.f90
        pldm_master.f90
        hamiltonians.f90
        calculate_pldm_redmat.f90
        twa_master.f90
        calculate_twa_redmat.f90
        kinds.f90
        windows.f90
        parameters.f90
        spectroscopy.f90
        calculate_pldm_absorption.f90
        sqc_master.f90
        calculate_twa_redmat.f90
        calculate_sqc_redmat.f90
        sqc_master.f90
        calculate_sqc_absorption.f90
        focusing.f90)

### add a library containing all of the source other than main.f90
add_library(generallib ${MQDS_GENERAL_LIB})

### If MPI was found, link the general library with MPI
if(MPI_FOUND)
    target_link_libraries(generallib MPI::MPI_Fortran)
endif()

add_executable(mqds.exe main.f90)

### If MPI was found, link the executable with generallib and MPI
if(MPI_FOUND)
    target_link_libraries(mqds.exe generallib MPI::MPI_Fortran)
    message(STATUS "MPI has been linked to executable")
### Otherwise, just link it with generallib
else()
    target_link_libraries(mqds.exe generallib)
endif()