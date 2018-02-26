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