!> This is the TWA master subroutine that determines
!! which TWA-based calculation to perform.
SUBROUTINE surface_hopping_master
  USE input_output
  USE harmonic_bath
  USE mapping_variables
  USE hamiltonians
  USE random_numbers
  USE windows
  IMPLICIT NONE

  ! Initialize random seed (with mype argument for MPI)
  CALL initialize_rn

  ! Do calculations in the adiabatic basis with harmonic bath
  IF ( bath == 'harmonic' ) THEN

     CALL initialize_bath

     IF ( basis == 'adiabatic' ) THEN

        CALL initialize_hel

!        IF ( calculation == 'redmat' ) CALL calculate_surface_hopping_redmat

        CALL finalize_hel

     END IF

     CALL finalize_bath

  END IF

END SUBROUTINE surface_hopping_master
