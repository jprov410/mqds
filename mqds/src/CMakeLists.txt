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
        sqc_master.f90)


add_library(generallib ${MQDS_GENERAL_LIB})

add_executable(mqds.exe main.f90)
target_link_libraries(mqds.exe generallib)
