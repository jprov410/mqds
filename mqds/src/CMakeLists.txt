add_subdirectory(general_src)
add_subdirectory(pldm_src)
add_subdirectory(sqc_src)
add_subdirectory(twa_src)

### Include the dircectories that contain libraries
include_directories(
        ${MQDS_BINARY_DIR}/mqds/src/general_src
        ${MQDS_BINARY_DIR}/mqds/src/pldm_src
        ${MQDS_BINARY_DIR}/mqds/src/sqc_src
        ${MQDS_BINARY_DIR}/mqds/src/twa_src
)
set(MQDS_LIBRARIES
        generallib
        pldmlib
        sqclib
        twalib)

add_executable(mqds.exe main.f90)

### If MPI was found, link the executable with generallib and MPI
if(MPI_FOUND)
    target_link_libraries(mqds.exe
            ${MQDS_LIBRARIES}
            MPI::MPI_Fortran
            )
    message(STATUS "MPI has been linked to executable")

    ### Otherwise, just link it with generallib
else()
    target_link_libraries(mqds.exe
            ${MQDS_LIBRARIES}
            )
endif()

# Copy executable to directory with submission/analysis scripts
set(MQDS_RUN_PATH ${CMAKE_HOME_DIRECTORY}/run/mqds_run/.)
add_custom_command(TARGET mqds.exe
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:mqds.exe> ${MQDS_RUN_PATH})
MESSAGE(STATUS "executable has been copied to ${MQDS_RUN_PATH}")