!> This module containes subroutines and functions
!! necessary for imaginary time path integral calculations
MODULE imaginary_time
    USE kinds
    IMPLICIT NONE

CONTAINS

    SUBROUTINE get_covariance_matrix( covariance_matrix, beta )
    USE parameters
    USE input_output
    IMPLICIT NONE
    INTEGER :: npts
    REAL(dp), INTENT(out) :: covariance_matrix(:,:,:)
    REAL(dp), INTENT(in) :: beta
    REAL(dp), ALLOCATABLE :: sd(:,:), omega(:)
    CHARACTER :: junk

    npts = 0
    OPEN(UNIT=10, FILE=SPECTRALIN)
    DO
        READ(10, *, ERR=11, END=11) junk
        npts=npts+1
    END DO
    11  CONTINUE ; CLOSE(10)

    ! Warn if there was no continuum spectral density provided
    IF (npts==0) THEN
        OPEN(unit=20,file=ERRORLOG)
        WRITE(20,*) 'No continuum spectral density file found!'
        CLOSE(20)
    END IF




    END SUBROUTINE get_covariance_matrix

END MODULE imaginary_time