set(MQDS_GENERAL_LIB
        general_src/input_output.f90
        general_src/unit_conversions.f90
        general_src/random_numbers.f90
        general_src/harmonic_bath.f90
        general_src/mapping_variables.f90
        pldm_src/pldm_master.f90
        general_src/hamiltonians.f90
        pldm_src/calculate_pldm_redmat.f90
        twa_src/twa_master.f90
        twa_src/calculate_twa_redmat.f90
        general_src/kinds.f90
        pldm_src/windows.f90
        general_src/parameters.f90
        general_src/spectroscopy.f90
        pldm_src/calculate_pldm_absorption.f90
        sqc_src/sqc_master.f90
        twa_src/calculate_twa_redmat.f90
        sqc_src/calculate_sqc_redmat.f90
        sqc_src/sqc_master.f90
        sqc_src/calculate_sqc_absorption.f90
        general_src/focusing.f90)

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