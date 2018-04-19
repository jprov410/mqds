!> This is the PLDM master subroutine that determines
!! which PLDM-based calculation to perform.

!> The Partially Linearized Density Matrix Method is a
!! mixed quantum-classical method, although the "quantum"
!! aspect comes at the cost of a semiclassical treatment.
!! See ( The Journal of Chemical Physics 135, 201101 (2011) )
!! for more information.

!> @image html pldm_schematic.pdf
!! @image latex pldm_schematic.pdf "Importance sampling scheme" width=0.5/linewidth
SUBROUTINE pldm_master_mpi
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

          IF ( calculation == 'redmat' ) THEN
              CALL calculate_pldm_redmat_mpi
          END IF

          IF ( calculation == 'absorption' ) THEN
              CALL calculate_pldm_absorption_mpi
          END IF

          IF ( calculation == 'nonlinear') THEN
              CALL calculate_pldm_nonlinear_mpi
          END IF

          CALL finalize_hel

      END IF

      CALL finalize_bath

  END IF


  ! Deallocate mapping variables
  CALL finalize_pldm_map

END SUBROUTINE pldm_master_mpi
