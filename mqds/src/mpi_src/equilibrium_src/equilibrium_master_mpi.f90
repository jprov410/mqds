!> This is the equilibrium master subroutine that
!! uses input to determine which type of equilibrium
!! reduced density matrix calculation to do.
!! Only currently calculates equilibrium site populations.
SUBROUTINE equilibrium_master_mpi
    USE kinds
    USE input_output
    USE random_numbers
    USE harmonic_bath
    USE hamiltonians
    USE mpi_variables
    IMPLICIT NONE

    CALL initialize_rn(mype)
    CALL initialize_bath
    CALL initialize_hel



    CALL calculate_equilibrium_site_mpi



    CALL finalize_hel
    CALL finalize_bath

END SUBROUTINE equilibrium_master_mpi