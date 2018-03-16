!> This is the MPI version of the main program that interprets the
!! input file and decides which methods to execute with mpi implemented.
PROGRAM mqds_mpi
  USE input_output
  USE mpi_variables
  IMPLICIT NONE

  ! Initialize mpi & communication
  CALL MPI_INIT(ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, npes, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, mype, ierr)

  CALL read_input
  
  IF ( MOD( nbstep, dump ) /= 0 ) THEN
      OPEN(UNIT=10, FILE=ERRORLOG)
      WRITE(10,*) 'Error: nbstep must be divisible by dump!'
      CLOSE(10)
      STOP
  END IF

  ! mapping variable methods
  IF ( method == 'pldm' ) THEN
      CALL pldm_master_mpi
  END IF

  IF ( method == 'twa' ) THEN
      CALL twa_master_mpi
  END IF

  IF ( method == 'sqc' ) THEN
      CALL sqc_master_mpi
  END IF

  IF ( method == 'ipldm' ) THEN
      CALL ipldm_master_mpi
  END IF

  CALL MPI_FINALIZE(ierr)

END PROGRAM mqds_mpi

