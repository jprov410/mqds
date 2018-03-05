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

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/hel.in
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/couplings.in
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/continuumsd.in
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})