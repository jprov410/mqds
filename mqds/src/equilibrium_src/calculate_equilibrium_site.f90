!> This subroutine calculates the exact equilibrium
!! reduced density matrix in the site basis using the
!! methodology described in
!! PHYSICAL REVIEW B 85, 115412 (2012)
SUBROUTINE calculate_equilibrium_site
    USE kinds
    USE random_numbers
    USE harmonic_bath
    USE hamiltonians
    USE linear_algebra
    USE unit_conversions
    USE imaginary_time
    USE input_output
    IMPLICIT NONE
    INTEGER :: istate
    REAL(dp) :: beta
    REAL(dp) :: eigenvectors( nstate, nstate )
    REAL(dp) :: eigenvalues( nstate )
    REAL(dp) :: covariance_matrix( nbstep, nbstep, nstate )
    REAL(dp) :: covariance_eigenvectors_inverse( nbstep, nbstep, nstate )
    REAL(dp) :: covariance_eigenvalues( nbstep, nstate )
    REAL(dp) :: equilibrium_populations( nstate )
    REAL(dp) :: redmat( nstate, nstate )
    REAL(dp) :: redmat_init( nstate, nstate )
    eigenvectors = hel ; eigenvalues = 0.0_dp ; covariance_matrix = 0.0_dp
    covariance_eigenvectors_inverse = 0.0_dp ; covariance_eigenvalues = 0.0_dp
    equilibrium_populations = 0.0_dp ; redmat = 0.0_dp ; redmat_init = 0.0_dp

    ! Set Boltzmann according to input temperature
    beta = 1.0_dp / ( temperature * convert('kelvin','au_energy') )

    ! Read the input electonic Hamiltonian
    CALL read_hel ; hel = hel * convert('wvnbr','au_energy')

    ! Initialize in the high temperature limit, all states equally populated
    DO istate = 1, nstate
        redmat_init( istate, istate ) = 1.0_dp / nstate
    END DO

    CALL get_covariance_matrix( covariance_matrix, beta )


END SUBROUTINE calculate_equilibrium_site