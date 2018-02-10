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
    FUNCTION dipole_operator( redmat ) RESULT( res )
!        USE input_output
        IMPLICIT NONE
        COMPLEX(dp) :: redmat(nstate, nstate)
        COMPLEX(dp) :: res(nstate, nstate)

        res = MATMUL( dipole, redmat )

    END FUNCTION dipole_operator

    ! Function to act dipole commutator on the reduced
    ! density matrix
    FUNCTION dipole_commutator( redmat ) RESULT( res )
!        USE input_output
        IMPLICIT NONE
        COMPLEX(dp) :: redmat(nstate, nstate)
        COMPLEX(dp) :: res(nstate, nstate)

        res = MATMUL( dipole, redmat ) - MATMUL( redmat, dipole )

    END FUNCTION dipole_commutator

    ! Function to find the linear response function for
    ! the absorption spectrum
    FUNCTION linear_response( redmat ) RESULT( res )
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

    END FUNCTION linear_response

    ! Subroutine to write the linear response function
    SUBROUTINE write_linear_response(method, resp_func, printstep)
        USE kinds
        IMPLICIT NONE
        INTEGER :: i, j, itime
        INTEGER :: tdim, namesize
        COMPLEX(dp), INTENT(in) :: resp_func(:)
        REAL(dp), INTENT(in) :: printstep
        CHARACTER(20) :: method
        CHARACTER(17) :: writeformat
        CHARACTER(50) :: filename
        tdim = SIZE(resp_func(:))

        ! Figure out how long the name is going to be
        namesize=LEN_TRIM(method)+16
        writeformat='(A,A)'

        ! create file for each redmat element and write
        WRITE(filename, writeformat) TRIM(method),'_linrespfunc.out'
        OPEN(UNIT=10, FILE=filename)
        DO itime=1, tdim
            WRITE(10,'(3(E13.6,2X))') printstep * (itime-1), resp_func(itime)
        END DO
        CLOSE(10)

    END SUBROUTINE write_linear_response



END MODULE spectroscopy