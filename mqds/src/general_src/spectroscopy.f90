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


END MODULE spectroscopy