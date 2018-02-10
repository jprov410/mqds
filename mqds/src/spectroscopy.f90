!
! This is a module that contains the necessary variables,
! subroutines and functions to perform linear and nonlinear
! spectroscopy calculations.
!
MODULE spectroscopy
    USE kinds
    USE input_output
    IMPLICIT NONE
    REAL(dp), ALLOCATABLE, PRIVATE :: dipole(:,:)

    !strange that wont compile if I USE  input_output in the
    ! dipole operator and dipole commutator functions defined below


CONTAINS

    ! Allocate the necessary arrays for performing
    ! a spectroscopy calculation and read input
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

    ! Deallocate the arrays from spectroscopy calculation
    SUBROUTINE finalize_spectroscopy
        IMPLICIT NONE

        DEALLOCATE( dipole )

    END SUBROUTINE finalize_spectroscopy

    ! Function to act dipole operator on the reduced
    ! density matrix
    FUNCTION dipole_operator(input) RESULT( res )
!        USE input_output
        IMPLICIT NONE
        COMPLEX(dp) :: input(nstate, nstate)
        COMPLEX(dp) :: res(nstate, nstate)

        res = MATMUL( dipole, input )

    END FUNCTION dipole_operator

    ! Function to act dipole commutator on the reduced
    ! density matrix
    FUNCTION dipole_commutator(input) RESULT( res )
!        USE input_output
        IMPLICIT NONE
        COMPLEX(dp) :: input(nstate, nstate)
        COMPLEX(dp) :: res(nstate, nstate)

        res = MATMUL( dipole, input ) - MATMUL( input, dipole )

    END FUNCTION dipole_commutator

END MODULE spectroscopy