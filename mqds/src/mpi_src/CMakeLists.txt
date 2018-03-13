if(MPI_Fortran_FOUND)
    MESSAGE(STATUS "MPI was found, building mpi executable")
    add_subdirectory(pldm_src)
    add_subdirectory(sqc_src)
    add_subdirectory(twa_src)

### Include the dircectories that contain libraries
    include_directories(
            ${MQDS_BINARY_DIR}/mqds/src/general_src
            ${MQDS_BINARY_DIR}/mqds/src/mpi_src/pldm_src
            ${MQDS_BINARY_DIR}/mqds/src/mpi_src/sqc_src
            ${MQDS_BINARY_DIR}/mqds/src/mpi_src/twa_src
    )
    add_library(mpilib mpi_variables.f90)

    set(MQDS_LIBRARIES_MPI
            generallib
            pldmlib_mpi
            sqclib_mpi
            twalib_mpi
            mpilib)

add_executable(mqds_mpi.exe main_mpi.f90)

    target_link_libraries(mqds_mpi.exe
            ${MQDS_LIBRARIES_MPI}
            MPI::MPI_Fortran
            )
    message(STATUS "MPI has been linked to executable")

# Copy executable to directory with submission/analysis scripts
    set(MQDS_RUN_PATH ${CMAKE_HOME_DIRECTORY}/run/mqds_run/.)
    add_custom_command(TARGET mqds_mpi.exe
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:mqds_mpi.exe> ${MQDS_RUN_PATH})
    MESSAGE(STATUS "MPI executable has been copied to ${MQDS_RUN_PATH}")
endif()