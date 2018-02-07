set(MQDS_GENERAL_LIB	kinds.f90
			parameters.f90
			input_output.f90
			unit_conversions.f90
			random_numbers.f90
			harmonic_bath.f90
			mapping_variables.f90
			pldm_master.f90
			hamiltonians.f90
			calculate_pldm_redmat.f90)
			

add_library(generallib ${MQDS_GENERAL_LIB})
target_link_libraries(generallib)

add_executable(mqds main.f90)
target_link_libraries(mqds generallib)
			

