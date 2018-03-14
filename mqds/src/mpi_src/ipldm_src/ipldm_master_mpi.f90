!
! This is the PLDM master subroutine that determines
!      which PLDM-based calculation to perform
!
SUBROUTINE ipldm_master
  USE input_output
  USE harmonic_bath
  USE mapping_variables
  USE hamiltonians
  USE random_numbers
  USE mpi_variables
  IMPLICIT NONE

  ! Initialize random seed (with mype argument for MPI)
  CALL initialize_rn(mype)

  ! Allocate arrays and initialize mapping variables
  CALL initialize_pldm_map

  ! Do calculations ins the diabatic basis with harmonic bath
  IF ( bath == 'harmonic' ) THEN

     CALL initialize_bath

     IF ( basis == 'diabatic' ) THEN

        CALL initialize_hel

        IF ( calculation == 'redmat' ) CALL calculate_ipldm_redmat_mpi

        !IF (calculation == 'absorption' ) CALL calculate_ipldm_absorption_mpi

        !IF (calculation == 'nonlinear') CALL calculate_pldm_nonlinear_mpi

        CALL finalize_hel

     END IF

     CALL finalize_bath

  END IF


  ! Deallocate mapping variables
  CALL finalize_pldm_map

END SUBROUTINE ipldm_master
