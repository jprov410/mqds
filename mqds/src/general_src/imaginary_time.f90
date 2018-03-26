!> This module containes subroutines and functions
!! necessary for imaginary time path integral calculations
MODULE imaginary_time
    USE kinds
    IMPLICIT NONE

CONTAINS

    !> This subroutine calculated the cavariance matrix for imaginary time
    !! path integral simulation for equilibrium reduced density matrix calculations.
    !! It is the imaginary time version of the force autocorrelation function of
    !! the harmonic bath. Each state has its own bath.
    SUBROUTINE get_covariance_matrix( covariance_matrix, beta )
        USE parameters
        USE input_output
        USE unit_conversions
        IMPLICIT NONE
        INTEGER :: npts, istate, ipt, istep, istept
        REAL(dp), INTENT(out) :: covariance_matrix( 0 : nbstep, 0 : nbstep, nstate )
        REAL(dp), INTENT(in) :: beta
        REAL(dp), ALLOCATABLE :: sd(:,:), freq(:)
        REAL(dp) :: dw, dt, t1, t2
        CHARACTER :: junk
        npts = 0 ; dw = 0.0_dp ; dt = beta / nbstep

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

        ALLOCATE( sd( nstate, npts ), freq( npts ) )

        OPEN(UNIT=10, FILE=SPECTRALIN)
        DO ipt=1, npts
            READ(10, *) freq(ipt), ( sd(istate,ipt), istate=1, nstate )
        END DO
        CLOSE(10)

        freq = freq * convert('wvnbr','au_ang_freq')
        sd = sd * convert('wvnbr','au_energy')

        dw = freq(2) - freq(1)

        DO istep = 0, nbstep
            t1 = istep * dw
            DO istept = 0, istep
                t2 = istept * dw
                DO ipt = 1 , npts
                    covariance_matrix( istep, istept, : ) = ( dw / pi ) &
                            * sd( :, ipt ) &
                            * COSH( ( beta * freq(ipt) ) / 2.0_dp - freq(ipt) * ABS(t1-t2) ) &
                            / SINH( beta * freq(ipt) / 2.0_dp )
                END DO
                covariance_matrix( istept, istep, : ) = covariance_matrix( istep, istept, : )
            END DO
        END DO


        DEALLOCATE( sd, freq )

    END SUBROUTINE get_covariance_matrix

END MODULE imaginary_time