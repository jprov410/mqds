!> This is a module that contains the necessary variables,
!! subroutines and functions to perform linear (and soon nonlinear)
!! spectroscopy calculations.
MODULE spectroscopy
    USE kinds
    USE input_output
    IMPLICIT NONE
    REAL(dp), ALLOCATABLE, PRIVATE :: dipole(:,:)

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

    !> Deallocate the arrays from spectroscopy calculation
    SUBROUTINE finalize_spectroscopy
        IMPLICIT NONE

        DEALLOCATE( dipole )

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

        !CHECK PLUS OR MINUS k2
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