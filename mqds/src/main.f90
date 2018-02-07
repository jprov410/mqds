!
! This is the main program that interprets the input file
! and decides which methods to execute
!
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

  IF ( method .EQ. 'pldm' ) CALL pldm_master
  

END PROGRAM mqds
