!> Main executable for the Modular Quantum Dynamics and Spectroscopy
!! package that reads the input file called "processed_run.in" and
!! calls the necessary based on the specified dynamics method.
PROGRAM mqds
  USE input_output
  IMPLICIT NONE
  
  CALL read_input
  
  IF ( MOD( nbstep, dump ) /= 0 ) THEN
     OPEN(UNIT=10, FILE=ERRORLOG)
     WRITE(10,*) 'Error: nbstep must be divisible by dump!'
     CLOSE(10)
     STOP
  END IF

  ! Mapping Variable Dynamics methods
  IF ( method == 'pldm' ) CALL pldm_master
  
  IF ( method == 'twa' ) CALL twa_master

  IF ( method == 'sqc' ) CALL sqc_master

  IF ( method == 'ipldm' ) CALL ipldm_master

  ! Equilibrium Reduced Density Matrix method


END PROGRAM mqds

