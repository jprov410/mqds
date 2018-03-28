!> This subroutine calculates the exact equilibrium
!! reduced density matrix in the site basis using the
!! methodology described in
!! PHYSICAL REVIEW B 85, 115412 (2012)
SUBROUTINE calculate_equilibrium_site_mpi
    USE kinds
    USE random_numbers
    USE harmonic_bath
    USE hamiltonians
    USE linear_algebra
    USE unit_conversions
    USE imaginary_time
    USE input_output
    USE mpi_variables
    IMPLICIT NONE
    INTEGER :: istate, itraj, istep
    REAL(dp) :: beta, dt, trace
    REAL(dp) :: hel_eigenvectors( nstate, nstate )
    REAL(dp) :: hel_eigenvalues( nstate )
    REAL(dp) :: covariance_matrix( 0 : nbstep, 0 : nbstep, nstate )
    REAL(dp) :: covariance_eigenvectors_inverse( 0 : nbstep, 0 : nbstep, nstate )
    REAL(dp) :: covariance_eigenvalues( 0 : nbstep, nstate )
    REAL(dp) :: temporary( 0 : nbstep, 0 : nbstep )
    REAL(dp) :: equilibrium_site_populations( nstate )
    REAL(dp) :: sum_equilibrium_site_populations( nstate )
    REAL(dp) :: redmat( nstate, nstate )
    REAL(dp) :: redmat_init( nstate, nstate )
    REAL(dp) :: system_propagator( nstate, nstate ), bath_propagator( nstate, nstate )
    REAL(dp) :: xi( 0 : nbstep, nstate )
    hel_eigenvectors = 0.0_dp ; hel_eigenvalues = 0.0_dp ; covariance_matrix = 0.0_dp
    covariance_eigenvectors_inverse = 0.0_dp ; covariance_eigenvalues = 0.0_dp
    equilibrium_site_populations = 0.0_dp ; redmat = 0.0_dp ; redmat_init = 0.0_dp
    system_propagator = 0.0_dp ; bath_propagator = 0.0_dp ; dt = 0.0_dp ; trace = 0.0_dp

    ! Set Boltzmann according to input temperature
    beta = 1.0_dp / ( temperature * convert('kelvin','au_energy') )
    dt = beta / nbstep

    ! Read the input electonic Hamiltonian
    CALL read_hel ; hel = hel * convert('wvnbr','au_energy')

    ! Initialize in the high temperature limit, all states equally populated
    DO istate = 1, nstate
        redmat_init( istate, istate ) = 1.0_dp / nstate
    END DO

    ! Initialize the imaginary time covariance matrix using the input continuum spectral density
    CALL get_covariance_matrix( covariance_matrix, beta )

    ! Get the electronic (system) eigenvectors and eigenvalues
    hel_eigenvectors = hel
    CALL diagonalize_real_symmetric( hel_eigenvectors, hel_eigenvalues )

    ! Use the system eigenstuff to calculate the system portion of the propagator
    ! to use with the Trotter factorization of the imaginary time propagator
    DO istate = 1, nstate
        system_propagator( istate, istate ) = EXP( -dt * hel_eigenvalues( istate ) )
    END DO
    system_propagator = MATMUL( system_propagator, TRANSPOSE( hel_eigenvectors ) )
    system_propagator = MATMUL( hel_eigenvectors, system_propagator )

    ! Get the eigenvectors and eigenvalues for the inverse of the covariance matrix
    covariance_eigenvectors_inverse = covariance_matrix
    DO istate = 1, nstate
        CALL diagonalize_real_symmetric( covariance_eigenvectors_inverse( :, :, istate ), &
                covariance_eigenvalues( :, istate ) )
        covariance_eigenvectors_inverse( :, :, istate ) =  TRANSPOSE( covariance_eigenvectors_inverse( :, :, istate ) )
    END DO

    ! Run Imaginary time trajectories
    DO itraj = 1, INT( ntraj / npes )
        redmat = redmat_init


        DO istate = 1, nstate
            xi( :, istate ) = gaussian_rn( xi( :, istate ) ) * DSQRT( covariance_eigenvalues( :, istate ) ) !/ dt
            xi( :, istate ) = MATMUL( covariance_eigenvectors_inverse( :, :, istate ), xi( :, istate ) )
        END DO

        ! Propagate the reduced density matrix with trotter in imaginary time
        DO istep = 0, nbstep

            ! Build current bath propagator
            DO istate = 1, nstate
                bath_propagator( istate, istate ) = EXP( -0.5_dp * dt * xi( istep, istate ) )
            END DO

            ! Trotter propagation in imaginary time
            redmat = MATMUL( bath_propagator, redmat )
            redmat = MATMUL( system_propagator, redmat )
            redmat = MATMUL( bath_propagator, redmat )

        END DO

        DO istate = 1, nstate
            equilibrium_site_populations( istate ) = equilibrium_site_populations( istate ) &
                    + redmat( istate, istate )
        END DO

    END DO

    trace = SUM( equilibrium_site_populations )
    equilibrium_site_populations = equilibrium_site_populations / trace

    result_size = SIZE( equilibrium_site_populations )

    CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
    IF ( mype /= 0 ) THEN
        CALL MPI_Send(equilibrium_site_populations, result_size, MPI_DOUBLE, &
                0, tag, MPI_COMM_WORLD, ierr)
    ELSE
        sum_equilibrium_site_populations = equilibrium_site_populations
        DO ipe=1, npes-1
            CALL MPI_recv(equilibrium_site_populations, result_size, MPI_DOUBLE,&
                    ipe, tag, MPI_COMM_WORLD, status,ierr)
            sum_equilibrium_site_populations = &
                    sum_equilibrium_site_populations + equilibrium_site_populations
        END DO
        sum_equilibrium_site_populations = sum_equilibrium_site_populations / DBLE(npes)
        CALL write_equilibrium_site_populations( sum_equilibrium_site_populations )
    END IF

END SUBROUTINE calculate_equilibrium_site_mpi