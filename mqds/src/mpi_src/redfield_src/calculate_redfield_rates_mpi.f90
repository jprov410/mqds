!> subroutine for the calculation of rates in the Redfield approximation
SUBROUTINE calculate_redfield_rates_mpi
    USE harmonic_bath
    USE hamiltonians
    USE input_output
    USE kinds
    USE linear_algebra
    USE parameters
    USE fourier
    USE unit_conversions
    IMPLICIT NONE
    INTEGER :: i,j,k
    REAL(dp) :: eigenvectors(nstate,nstate), exciton_energies(nstate)
    REAL(dp) :: exciton_reorg(nstate)
    REAL(dp) :: site_bath_cf(nstate, 0:nbstep)
    REAL(dp) :: exciton_bath_cf(nstate, nstate, 0:nbstep)
    REAL(dp) :: ft_exciton_bath_cf(nstate,nstate)
    REAL(dp) :: beta, dt, dw_electronic, freq
    exciton_energies = 0.0_dp ; eigenvectors = 0.0_dp
    exciton_reorg = 0.0_dp ; dw_electronic = 0.0_dp
    site_bath_cf = 0.0_dp ; ft_exciton_bath_cf = 0.0_dp
    exciton_bath_cf = 0.0_dp ; freq = 0.0_dp

    beta = 1.0_dp / ( temperature * convert('kelvin','au_energy') )
    runtime = runtime * convert('fs','au_time')
    dt = runtime / DBLE(nbstep)

    !> Read in the system Hamiltonian
    CALL read_hel
    hel = hel * convert('wvnbr','au_energy')

    !> Diagonalize to get Eigenvectors and Energy Eigenvalues
    eigenvectors = hel
    CALL diagonalize_real_symmetric(eigenvectors, exciton_energies)

    !> Calculate the exciton basis reorganization energy
    DO i = 1, nstate
        DO j = 1, nstate
            exciton_reorg(i) = exciton_reorg(i) + &
                    site_reorg(j) * eigenvectors(j,i) ** 4
        END DO
    END DO

!    !> Calculate size of electronic bath frequency grid
!    dw_electronic = electronic_bath_freq(2) - electronic_bath_freq(1)

    !> Calculate the bath correlation functions for the modes coupled to sites
    DO i = 1, nstate
        site_bath_cf(i,:) = redfield_bath_correlation_function(site_sd(i,:), &
                site_bath_freq,beta)
    END DO

    !> Linearly combine the site bath correlation functions to get the exciton ones
    DO i=1,nstate
        DO j=1,nstate
            DO k=1,nstate
                exciton_bath_cf(i,j,:) = exciton_bath_cf(i,j,:) &
                        + eigenvectors(k,i)**2 * site_bath_cf(k,:) * eigenvectors(k,j)**2
            END DO
        END DO
    END DO


    !> Fourier Transform the bath correlation function in exciton basis
    DO i=1, nstate
        DO j=1, nstate
                freq = exciton_energies(i) - exciton_energies(j)
                ft_exciton_bath_cf(i,j) = ft_at_omega( exciton_bath_cf(i,j,:), freq, dt )
        END DO
    END DO




END SUBROUTINE calculate_redfield_rates_mpi