!
! This is the TWA master subroutine that determines
!      which TWA-based calculation to perform
!
SUBROUTINE sqc_master
  USE input_output
  USE harmonic_bath
  USE mapping_variables
  USE hamiltonians
  USE random_numbers
  USE windows
  IMPLICIT NONE

  ! Initialize random seed (with mype argument for MPI)
  CALL initialize_rn

  ! Allocate arrays and initialize mapping variables
  ! sqc and twa use the same initialization of mapping
  ! variable arrays
  CALL initialize_twa_map
  CALL initialize_windows

  ! Do calculations ins the diabatic basis with harmonic bath
  IF ( bath == 'harmonic' ) THEN

     CALL initialize_bath

     IF ( basis == 'diabatic' ) THEN

        CALL initialize_hel

        IF ( calculation == 'redmat' ) CALL calculate_sqc_redmat

        IF ( calculation == 'absorption' ) CALL calculate_sqc_absorption

        CALL finalize_hel

     END IF

     CALL finalize_bath

  END IF

  ! Deallocate mapping variables
  CALL finalize_twa_map
  CALL finalize_windows

END SUBROUTINE sqc_master
