### Include the directories with .o .mod files
message(STATUS "building unit tests")

INCLUDE_DIRECTORIES(
        ${MQDS_BINARY_DIR}/mqds/src/general_src
        ${MQDS_BINARY_DIR}/mqds/src/pldm_src
        ${MQDS_BINARY_DIR}/mqds/src/sqc_src
        ${MQDS_BINARY_DIR}/mqds/src/twa_src)

add_executable(test_unit_conversions test_unit_conversions.f90)
target_link_libraries(test_unit_conversions generallib)
add_test(UNIT_CONVERSIONS test_unit_conversions)

add_executable(test_random_numbers test_random_numbers.f90)
target_link_libraries(test_random_numbers generallib)
add_test(RANDOM_NUMBERS test_random_numbers)

add_executable(test_focusing test_focusing.f90)
target_link_libraries(test_focusing generallib)
add_test(FOCUSING test_focusing)

add_executable(test_hamiltonians test_hamiltonians.f90)
target_link_libraries(test_hamiltonians generallib)
add_test(HAMILTONIANS test_hamiltonians)

add_executable(test_kinds test_kinds.f90)
target_link_libraries(test_kinds generallib)
add_test(KINDS test_kinds)

add_executable(test_harmonic_bath test_harmonic_bath.f90)
target_link_libraries(test_harmonic_bath generallib)
add_test(HARMONIC_BATH test_harmonic_bath)

add_executable(test_mapping_variables test_mapping_variables.f90)
target_link_libraries(test_mapping_variables generallib)
add_test(MAPPING_VARIABLES test_mapping_variables)

add_executable(test_windows test_windows.f90)
target_link_libraries(test_windows generallib)
add_test(WINDOWS test_windows)

add_executable(test_spectroscopy test_spectroscopy.f90)
target_link_libraries(test_spectroscopy generallib)
add_test(SPECTROSCOPY test_spectroscopy)


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/hel.in
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/couplings.in
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/continuumsd.in
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/dipole.in
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})