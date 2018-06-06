add_executable(fourier_transform.exe fourier_transform.f90)

# Copy executable to directory with submimssion/analysis scripts
set(MQDS_RUN_PATH ${CMAKE_HOME_DIRECTORY}/run/analysis/.)
add_custom_command(TARGET fourier_transform.exe
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:fourier_transform.exe> ${MQDS_RUN_PATH})
message(STATUS "fourier_transform.exe executable will be copied to ${MQDS_RUN_PATH}")