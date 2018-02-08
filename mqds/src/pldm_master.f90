!
! This is the PLDM master subroutine that determines
!      which PLDM-based calculation to perform
!
SUBROUTINE pldm_master
  USE input_output
  USE harmonic_bath
  USE mapping_variables
  USE hamiltonians
  USE random_numbers
  IMPLICIT NONE
  INTEGER :: time

  ! Initialize random seed (with mype argument for MPI)
  CALL initialize_rn

  ! Allocate arrays and initialize mapping variables
  CALL initialize_pldm_map

  ! Do calculations ins the diabatic basis with harmonic bath
  IF ( bath .EQ. 'harmonic' ) THEN

     CALL initialize_bath

     IF ( basis .EQ. 'diabatic' ) THEN        

        CALL initialize_hel

        IF ( calculation .EQ. 'redmat' ) CALL calculate_pldm_redmat

        !IF (calculation .EQ. 'absorption') CALL calculate_pldm_absorption

        !IF (calculation .EQ. 'nonlinear') CALL calculate_pldm_nonlinear

        CALL finalize_hel

     END IF

     CALL finalize_bath

  END IF


  ! Deallocate mapping variables
  CALL finalize_pldm_map

END SUBROUTINE pldm_master
