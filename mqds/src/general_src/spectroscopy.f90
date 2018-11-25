!> This is a module that contains the necessary variables,
!! subroutines and functions to perform linear (and soon nonlinear)
!! spectroscopy calculations.
MODULE spectroscopy
    USE kinds
    USE input_output
    IMPLICIT NONE
    REAL(dp), ALLOCATABLE, PRIVATE :: dipole(:,:)
    REAL(dp), PUBLIC :: orientation(40,3)
    REAL(dp), ALLOCATABLE, PRIVATE :: dip_xyz(:, :)
    REAL(dp), ALLOCATABLE, PUBLIC :: dip_weight(:, :)

CONTAINS

    !> Allocate the necessary arrays for performing
    !! a spectroscopy calculation and read input from
    !! dipole.in
    SUBROUTINE initialize_spectroscopy
        USE input_output
        IMPLICIT NONE
        INTEGER :: i, j

        ALLOCATE( dipole(nstate, nstate) )

        !read input dipole matrix elements
        OPEN(UNIT=10, FILE=DIPOLEIN)
        DO i=1, nstate
            READ(10,*) ( dipole(j,i), j=1, nstate )
        END DO
        CLOSE(10)

    END SUBROUTINE initialize_spectroscopy

    SUBROUTINE initialize_3d_dipole_mag
        USE input_output
        IMPLICIT NONE
        INTEGER i

        ALLOCATE(dip_xyz(nstate-1,3))
        allocate(dipole(nstate,nstate))

        open(unit = 10 , file = 'dip_xyz.in')
        DO i=1,nstate-1
            read(10,*) dip_xyz(i,1),dip_xyz(i,2),dip_xyz(i,3)
        END DO
        close(10)

        dipole = 0.d0

        DO i=1,nstate-1
            dipole(1,i+1) = sqrt( dip_xyz(i,1)**2 + dip_xyz(i,2)**2 + dip_xyz(i,3)**2 )
            dipole(i+1,1) = dipole(1,i+1)
        END DO

    END SUBROUTINE initialize_3d_dipole_mag

    SUBROUTINE read_3d_dipole
        USE input_output
        IMPLICIT NONE
        INTEGER :: i

        ALLOCATE(dip_xyz(nstate-1,3))
        allocate(dipole(nstate,nstate))

        open(unit = 10 , file = 'dip_xyz.in')
        DO i=1,nstate-1
            read(10,*) dip_xyz(i,1),dip_xyz(i,2),dip_xyz(i,3)
        END DO
        close(10)

    END SUBROUTINE read_3d_dipole

    SUBROUTINE calculate_oriented_dipole( iori )
        IMPLICIT NONE
        INTEGER :: iori
        INTEGER :: i

        dipole = 0.d0
        DO i = 1, nstate-1
            dipole(1,i+1) = dot_product( orientation(iori,:), dip_xyz(i,:) )
            dipole(i+1,1) = dipole(1,i+1)
        END DO

    END SUBROUTINE calculate_oriented_dipole

    SUBROUTINE generate_orientations
        USE parameters
        IMPLICIT NONE
        INTEGER :: i,j
        REAL(dp) :: offset, phi, increment, radius,y
        offset = 2.0_dp / SIZE(orientation(:,1))
        increment = pi * (3.0_dp - sqrt(5.0_dp))
        radius = 0.0_dp ; y = 0.0_dp


        open(unit=10, file='test_orientation.dat')

        DO i =1, SIZE(orientation(:,1))
            y = ( (i-1) * offset - 1.0_dp ) + 0.5_dp * offset
            radius = sqrt( 1 - y**2 )
            phi =  MOD(i, SIZE(orientation(:,1))) * increment

            orientation(i,1) = COS(phi) * radius
            orientation(i,2) = y
            orientation(i,3) = SIN(phi) * radius

            write(10,*) (orientation(i,j),j=1,3)
        END DO
        close(10)
    END SUBROUTINE generate_orientations

    SUBROUTINE get_dip_weight(istate, istatet)
        IMPLICIT NONE
        INTEGER :: i,j
        integer, intent(in) :: istate,istatet
        dipole=0.0

        !DO i=1, nstate-1
        !    dipole(1,i+1) = 1.0_dp
        !    dipole(i+1,1) = 1.0_dp
        !END DO

        !print*, istate, istatet

        IF (allocated(dip_weight) .eqv. .false.) allocate(dip_weight(nstate,nstate))
        dip_weight = 0.0_dp
        !DO i=1,nstate
        DO j=1,nstate-1
            dip_weight(j+1,1) = DOT_PRODUCT(dip_xyz(istate-1,:),dip_xyz(j,:))
            dip_weight(1,j+1) = DOT_PRODUCT(dip_xyz(istate-1,:),dip_xyz(j,:))
        END DO
        !END DO

    END SUBROUTINE get_dip_weight

    SUBROUTINE get_cd_weight(istate,istatet)
        USE linear_algebra
        IMPLICIT NONE
        INTEGER :: i,j
        INTEGER, INTENT(in) :: istate, istatet
        REAL(dp) :: centmass(nstate-1,3)
        REAL(dp) :: cent_mass_diff(3), trans_dip_cross(3)
        dipole=0.0
        centmass(1,1) = 29.550079; centmass(1,2) = 11.225173; centmass(1,3) = 8.908127
        centmass(2,1) = 23.205517; centmass(2,2) = 9.151394; centmass(2,3) = -2.471576
        centmass(3,1) = 14.429798; centmass(3,2) = -4.427753; centmass(3,3) = -8.275462
        centmass(4,1) = 16.496734; centmass(4,2) = -13.705020; centmass(4,3) = -2.681825
        centmass(5,1) = 24.547451; centmass(5,2) = -12.347156; centmass(5,3) = 7.926116
        centmass(6,1) = 22.627356; centmass(6,2) = -2.492773; centmass(6,3) = 13.961508
        centmass(7,1) = 15.863782; centmass(7,2) = -1.605234; centmass(7,3) = 2.861800
        centmass(8,1) = 20.490587; centmass(8,2) = 20.137727; centmass(8,3) = 13.569722
        !centmass = centmass*1.E-10_dp

        IF (allocated(dip_weight) .eqv. .false.) then
            allocate(dip_weight(nstate,nstate))
            !dip_xyz = dip_xyz * 3.33564E-30_dp
        END IF

        dip_weight = 0.0_dp
        DO j=1,nstate-1
            trans_dip_cross(:) = cross_product_3d(dip_xyz(j,:),dip_xyz(istate-1,:))
            cent_mass_diff(:) = centmass(j,:) - centmass(istate-1,:)
            dip_weight(j+1,1) = DOT_PRODUCT(cent_mass_diff(:),trans_dip_cross(:))
            dip_weight(1,j+1) = DOT_PRODUCT(cent_mass_diff(:),trans_dip_cross(:))
        END DO

    END SUBROUTINE get_cd_weight

    SUBROUTINE unit_dipoles()
        IMPLICIT NONE
        INTEGER :: i
        dipole = 0.0_dp
        do i = 2, nstate
            dipole(1,i) = 1.0_dp
            dipole(i,1) = 1.0_dp
        end do

    END SUBROUTINE unit_dipoles

    SUBROUTINE gaussian_fluctuating_site_energy(hel)
        use unit_conversions
        use random_numbers
        implicit none
        INTEGER :: i
        REAL(dp) :: fwhm(nstate) !in wavenumbers
        REAL(dp) :: hel(nstate,nstate)
        REAL(dp) :: gran(nstate)
        fwhm(1) = 0.0_dp !ground state
        fwhm(2) = 60.0_dp * convert('wvnbr','au_energy')
        fwhm(3) = 100.0_dp * convert('wvnbr','au_energy')
        fwhm(4) = 60.0_dp * convert('wvnbr','au_energy')
        fwhm(5) = 60.0_dp * convert('wvnbr','au_energy')
        fwhm(6) = 120.0_dp * convert('wvnbr','au_energy')
        fwhm(7) = 120.0_dp * convert('wvnbr','au_energy')
        fwhm(8) = 120.0_dp * convert('wvnbr','au_energy')
        fwhm(9) = 120.0_dp * convert('wvnbr','au_energy')

        gran = gaussian_rn(gran)
        DO i=1,nstate
            hel(i,i) = hel(i,i) + gran(i)*fwhm(i)/sqrt(8.0_dp*log(2.))
        END DO

    END SUBROUTINE gaussian_fluctuating_site_energy

    !> Deallocate the arrays from spectroscopy calculation
    SUBROUTINE finalize_spectroscopy
        IMPLICIT NONE

        IF (allocated(dipole) .eqv. .true.) then
            DEALLOCATE( dipole )
        END IF

        IF (allocated(dip_xyz) .eqv. .true.) then
            DEALLOCATE( dip_xyz )
        END IF

        IF (allocated(dip_weight) .eqv. .true.) then
            DEALLOCATE( dip_weight )
        END IF

        IF (allocated(dip_weight) .eqv. .true.) deallocate(dip_weight)
    END SUBROUTINE finalize_spectroscopy

    !> Function to act dipole operator on the reduced
    !! density matrix from the left
    FUNCTION dipole_operator( redmat ) RESULT( res )
        IMPLICIT NONE
        COMPLEX(dp) :: redmat(nstate, nstate)
        COMPLEX(dp) :: res(nstate, nstate)

        res = MATMUL( dipole, redmat )

    END FUNCTION dipole_operator

    !> Function to act dipole commutator on the reduced
    !! density matrix from the left
    FUNCTION dipole_commutator( redmat ) RESULT( res )
        IMPLICIT NONE
        COMPLEX(dp) :: redmat(nstate, nstate)
        COMPLEX(dp) :: res(nstate, nstate)

        res = MATMUL( dipole, redmat ) - MATMUL( redmat, dipole )

    END FUNCTION dipole_commutator

    !> Function to find the linear response function for
    !! calculation of the absorption spectrum
    FUNCTION system_trace( redmat ) RESULT( res )
        IMPLICIT NONE
        COMPLEX(dp) :: redmat(nstate, nstate, 0 : nbstep / dump)
        COMPLEX(dp) :: res(0 : nbstep / dump)
        INTEGER :: istate, itime
        res = (0.0_dp, 0.0_dp)

        DO itime=0, ( nbstep / dump )
            DO istate=1, nstate
                res(itime) = res(itime) + redmat(istate, istate, itime)
            END DO
        END DO

    END FUNCTION system_trace

    !> Subroutine that applies the dipole commutator the current reduced
    !! density matrix and then performs the focusing procedure on the
    !! resulting matrix to use as the sole initially-occupied states
    !! for subsequent propagation
    SUBROUTINE dipole_focus( x, p, xt, pt, weight, k_sign )
        USE focusing
        USE mapping_variables
        IMPLICIT NONE
        INTEGER :: i, j, minf, maxf, minb, maxb
        INTEGER, INTENT(out) :: k_sign
        REAL(dp), INTENT(in) :: x(nstate), p(nstate)
        REAL(dp), INTENT(in) :: xt(nstate), pt(nstate)
        COMPLEX(dp), INTENT(inout) :: weight
        COMPLEX(dp) :: matrix( nstate, nstate )

        ! Calculate \mu^{\times} \rho
        matrix = pldm_redmat( x, p, xt, pt )
        matrix = dipole_commutator( matrix )

        ! Focus and assign new values to initstate
        ! and initstatet
        CALL focus_redmat( matrix, weight )


        !FIND MIN AND MAX STATE FORWARD STATE BEFORE t=t1 DIPOLE
        DO minf = 1, nstate, 1
            IF ( x(minf) /= 0.0_dp .OR. p(minf) /= 0.0_dp ) THEN
                EXIT
            END IF
        END DO

        DO maxf = nstate, 1, -1
            IF ( x(maxf) /= 0.0_dp .OR. p(maxf) /= 0.0_dp ) THEN
                EXIT
            END IF
        END DO

        !FIND MIN AND MAX STATE BACKWARD BEFORE t=t1 DIPOLE
        DO minb = 1, nstate
            IF ( xt(minb) /= 0.0_dp .OR. pt(minb) /= 0.0_dp ) THEN
                EXIT
            END IF
        END DO

        DO maxb = nstate, 1, -1
            IF ( xt(maxb) /= 0.0_dp .OR. pt(maxb) /= 0.d0 ) THEN
                EXIT
            END IF
        END DO

        !CHECK PLUS OR MINUS k
        IF ( initstate > maxf .OR. initstatet < minb ) THEN
            k_sign = 1
        END IF

        IF ( initstate < minf .OR. initstatet > maxb ) THEN
            k_sign = -1
        END IF

    END SUBROUTINE dipole_focus

    !> Function that uses the sign of incoming wavevectors along with the current
    !! mapping variable values to calculate the different contributions to the
    !! total 3rd order nonlinear response function
    FUNCTION nonlinear_response( x, p, xt, pt, wt, k1, k2, k3 ) RESULT( res )
        USE mapping_variables
        USE parameters
        IMPLICIT NONE
        INTEGER :: k1, k2, k3, diagram, i, j
        REAL(dp) :: x( nstate ), p( nstate )
        REAL(dp) :: xt( nstate ), pt( nstate )
        COMPLEX(dp) :: wt
        COMPLEX(dp) :: res( 8 )
        COMPLEX(dp) :: rhotilde( nstate, nstate )
        res = 0.0_dp ; diagram = 1

        ! Figure out which of the responses we are recording
        ! ( K1 = -1,+1,+1 ; +1,-1,-1 ) REPHASING
        ! ( K2 = +1,-1,+1 ; -1,+1,-1 ) NONREPHASING
        ! ( K3 = +1,+1,-1 ; -1,-1,+1 ) DOUBLE QUANTUM
        ! ( K4 = +1,+1,+1 ; -1,-1,-1 ) NOT SURE WHAT THIS IS CALLED

        ! REPHASING SIGNALS
        IF ( k1 == -1 .AND. k2 == +1 .AND. k3 == +1 ) THEN
            diagram = 1
        END IF
        IF ( k1 == +1 .AND. k2 == -1 .AND. k3 == -1 ) THEN
            diagram = 2
        END IF

        ! NONREPHASING SIGNALS
        IF ( k1 == +1 .AND. k2 == -1 .AND. k3 == +1 ) THEN
            diagram = 3
        END IF
        IF ( k1 == -1 .AND. k2 == +1 .AND. k3 == -1 ) THEN
            diagram = 4
        END IF

        ! DOUBLE QUANTUM SIGNALS
        IF ( k1 == +1 .AND. k2 == +1 .AND. k3 == -1 ) THEN
            diagram = 5
        END IF
        IF ( k1 == -1 .AND. k2 == -1 .AND. k3 == +1 ) THEN
            diagram = 6
        END IF

        ! OTHER SIGNALS (K4)
        IF ( k1 == +1 .AND. k2 == +1 .AND. k3 == +1 ) THEN
            diagram = 7
        END IF
        IF ( k1 == -1 .AND. k2 == -1 .AND. k3 == -1 ) THEN
            diagram = 8
        END IF

        ! Get matrix
        rhotilde = pldm_redmat( x, p, xt, pt )

        ! Dipole operate on matrix and assign Monte-Carlo weight
        rhotilde = wt * dipole_operator( rhotilde )

        ! Find response function at (t1, t2, t3)
        DO i = 1, nstate
            res( diagram ) = res( diagram ) + rhotilde( i, i )
        END DO

        ! Multiply by (i / hbar)^3
        res = - eye * res

    END FUNCTION nonlinear_response

END MODULE spectroscopy