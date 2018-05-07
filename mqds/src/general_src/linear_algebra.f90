!> This module containes linear algebra subroutines
!! such as matrix diagonalization and inversion
MODULE linear_algebra
    USE kinds
    IMPLICIT NONE

CONTAINS

    !> Diagonalized a real symmetric matrix and returns the eigenvalues and the
    !! eigenvectors associated with that matrix
    SUBROUTINE diagonalize_real_symmetric( eigenvectors, eigenvalues )
        USE kinds
        IMPLICIT NONE
        REAL(dp), INTENT(out) :: eigenvalues( : )
        REAL(dp), INTENT(inout) :: eigenvectors( : , : )
        REAL(dp), ALLOCATABLE :: work( : )
        INTEGER :: info, lwork
        CHARACTER(1) :: jobz='V', uplo='U'
        lwork = 3 * SIZE(eigenvalues) - 1

        IF ( ALLOCATED( work ) .EQV. .TRUE. ) THEN
            DEALLOCATE( work )
        END IF
        ALLOCATE( work(lwork) )

        CALL dsyev( jobz, uplo, SIZE(eigenvalues), eigenvectors,&
                SIZE(eigenvalues), eigenvalues, work, lwork, info)

        IF ( ALLOCATED(work) .EQV. .TRUE. ) THEN
            DEALLOCATE( work )
        END IF

    END SUBROUTINE diagonalize_real_symmetric

    !> Function to perform the trace of a real square matrix
    FUNCTION trace_real_matrix( matrix ) RESULT( res )
        USE kinds
        IMPLICIT NONE
        REAL(dp) :: matrix(:,:), res
        INTEGER :: i
        res = 0.0_dp
        IF ( SIZE( matrix( : ,1 ) ) /= SIZE( matrix( 1 ,: ) ) ) THEN
            WRITE(*,*) 'matrix sent to trace function is not square'
            STOP 1
        END IF

        DO i = 1 , SIZE( matrix(:,1) )
            res = res + matrix(i,i)
        END DO

    END FUNCTION trace_real_matrix

    !> Function to perform the trace of a real square matrix
    FUNCTION trace_complex_matrix( matrix ) RESULT( res )
        USE kinds
        IMPLICIT NONE
        COMPLEX(dp) :: matrix(:,:), res
        INTEGER :: i
        res = 0.0_dp
        IF ( SIZE( matrix( : ,1 ) ) /= SIZE( matrix( 1 ,: ) ) ) THEN
            WRITE(*,*) 'matrix sent to trace function is not square'
            STOP 1
        END IF

        DO i = 1 , SIZE( matrix(:,1) )
            res = res + matrix(i,i)
        END DO

    END FUNCTION trace_complex_matrix

END MODULE linear_algebra