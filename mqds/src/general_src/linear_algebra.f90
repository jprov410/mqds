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

    SUBROUTINE invert_real_symmetric
        USE kinds
        IMPLICIT NONE



    END SUBROUTINE invert_real_symmetric


END MODULE linear_algebra