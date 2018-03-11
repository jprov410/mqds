PROGRAM fourier_transform
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: fs2wvnbr = 33356.41
    DOUBLE PRECISION, PARAMETER :: pi = 3.1415926536
    DOUBLE PRECISION :: dt1,dt2,dt3,dummy1,dummy2
    DOUBLE PRECISION :: w1_min, w2_min, w3_min
    DOUBLE PRECISION :: w1_max, w2_max, w3_max
    DOUBLE PRECISION :: dw1, dw2, dw3, w1, w2, w3
    INTEGER :: calc_type
    INTEGER :: itime1,itime2,itime3
    INTEGER :: ntime1,ntime2,ntime3
    INTEGER :: iw1,iw2,iw3
    INTEGER :: nw1,nw2,nw3
    CHARACTER(30) :: filename,junk
    CHARACTER(33) :: outfile
    COMPLEX(KIND=16), ALLOCATABLE :: lin_resp(:)
    COMPLEX(KIND=16), ALLOCATABLE :: ft_lin_resp(:)
    COMPLEX(KIND=16), ALLOCATABLE :: nonlin_resp(:,:,:)
    COMPLEX(KIND=16), ALLOCATABLE :: ft_nonlin_resp(:)
    ntime1 = 0 ; ntime2 = 0 ; ntime3 = 0


    ! FIGURE OUT WHICH KIND OF RESPONSE TO ANALYZE
    WRITE(*,*) 'WELCOME, PLEASE CHOOSE ONE OF THE FOLLOWING OPTIONS:'
    WRITE(*,*) 'LINEAR ABSORPTION SPECTRUM:                       1 '
    WRITE(*,*) '3rd ORDER NONLINEAR SPECTRUM:                     2 '
    READ(*,*) calc_type

    ! COMPUTE 1D FOURIER TRANSFORM
    IF ( calc_type == 1 ) THEN
        WRITE(*,*) 'INSERT LINEAR RESPONSE FILE NAME'
        WRITE(*,*) 'EXPECTING FORMAT time(fs), Re(resp), Im(resp)'
        READ(*,*) filename
        ! CHOOSE OUTPUT FILE NAME
        WRITE( outfile, '(A,A)') 'FT_',filename

        ! FIGURE OUT NUMBER OF TIME POINTS
        OPEN(UNIT=10, FILE=filename)
        88  CONTINUE
        READ(10,*,END=99,ERR=99) junk
        ntime1 = ntime1 + 1 ; GO TO 88
        99 CONTINUE ; REWIND(10)
        ALLOCATE( lin_resp( 0 : ntime1-1 ) )

        ! READ RESPONSE FUNCTION
        DO itime1 = 0, ntime1-1
            READ(10,*) dt1, dummy1, dummy2
            lin_resp(itime1) = CMPLX( dummy1, dummy2 )
        END DO
        CLOSE(10)
        dt1 = dt1 / (ntime1-1)
        WRITE(*,*) 'FOUND', ntime1, 'TIME POINTS, WITH dt=', dt1

        ! APPLY THE WINDOW FUNCTION
        WRITE(*,*) 'INTRODUCING COSINE WINDOW FUNCTION TO RESPONSE FUNCTION'
        DO itime1 = 0, ntime1-1
            lin_resp(itime1) = lin_resp(itime1) &
                    * COS( pi * itime1/ (2 * (ntime1-1)))
        END DO

        ! CHOOSE FREQUENCY RANGE AND NUMBER OF POINTS
        WRITE(*,*) 'INSERT FREQUENCY RANGE IN WAVENUMBERS (w_min, w_max, nw)'
        READ(*,*) w1_min, w1_max, nw1
        w1_min = w1_min * 2.d0 * pi / fs2wvnbr ; w1_max = w1_max * 2.d0 * pi / fs2wvnbr
        dw1 = (w1_max - w1_min) / DBLE(nw1)
        ALLOCATE( ft_lin_resp( 0 : nw1 ) )
        ft_lin_resp = 0.d0

        ! DO FOURIER TRANSFORM AT EACH FREQUENCY
        DO iw1=0, nw1
            w1 = w1_min + DBLE(iw1) * dw1
            ft_lin_resp( iw1 ) = ft( lin_resp, w1, dt1 )
        END DO

        ! WRITE THE FOURIER TRANSFORMED FUNCTION
        OPEN(UNIT=20, FILE=outfile)
        DO iw1 = 0, nw1
            WRITE(20,*) (w1_min + iw1*dw1) * fs2wvnbr / (2.d0 * pi),&
                    REAL(ft_lin_resp(iw1)), AIMAG(ft_lin_resp(iw1))
        END DO
        CLOSE(20)
    END IF


    IF ( calc_type == 2 ) THEN
        WRITE(*,*) 'APOLOGIES, THE NONLINEAR PORTION OF THE CODE HAS &
                BEEN IMPLEMENTED YET :( ...'
        STOP
    END IF


    IF ( ALLOCATED(lin_resp) .EQV. .TRUE. ) THEN
        DEALLOCATE(lin_resp)
    END IF
    IF ( ALLOCATED(nonlin_resp) .EQV. .TRUE. ) THEN
        DEALLOCATE(nonlin_resp)
    END IF

CONTAINS
    FUNCTION ft( response, w, dt ) RESULT( res )
        IMPLICIT NONE
        INTEGER :: t
        DOUBLE PRECISION, INTENT(IN) :: w, dt
        COMPLEX(KIND=16) :: res
        COMPLEX(KIND=16), PARAMETER :: eye = (0.d0, 1.d0)
        COMPLEX(KIND=16), INTENT(IN) :: response(:)
        res = 0.d0

        DO t = 1, SIZE(response)
            res = res + dt * exp( eye * w * (t-1) * dt ) * response(t)
        END DO

    END FUNCTION ft

END PROGRAM fourier_transform