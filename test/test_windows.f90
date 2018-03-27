!> Tests the windows module for SQC calculations.
PROGRAM test_windows
    USE windows
    USE random_numbers
    USE kinds
    USE input_output
    IMPLICIT NONE
    REAL(dp) :: x(2), p(2)
    REAL(dp), PARAMETER :: pi = 3.1415926536_dp
    COMPLEX(dp) :: test(2,2,0:1)
    nbstep = 1 ; dump = 1
    window = 0.5_dp ; nstate = 2
    initstate = 1 ; initstatet = 1
    zpe = 0.5_dp ; test = 0.0_dp
    windowshape = square
    CALL initialize_windows

    CALL sample_sqc_map(x,p)
    IF ( 0.5_dp * (x(1)**2 + p(1)**2) - zpe >= 1.0_dp + window &
            .OR. 0.5_dp * (x(1)**2 + p(1)**2) - zpe < 1.0_dp - window ) THEN
        WRITE(*,*) 'Initialization of windowed action not working correctly'
        STOP 1
    END IF
    IF ( -ATAN2( p(1), x(1) ) >= pi &
            .OR. -ATAN2( p(1), x(1) ) < -pi ) THEN
        WRITE(*,*) 'Initialization of angle variable not working correctly'
        STOP 1
    END IF

    test(:,:,0) = sqc_redmat(x,p,0) + sqc_redmat(x,p,0)
    test(:,:,1) = sqc_redmat(x,p,1)

    CALL normalize_sqc_redmat( test )
    IF ( test(1,1,0) /= 1.0_dp ) THEN
        WRITE(*,*) 'SQC redmat or normalization is not working correctly'
        STOP 1
    END IF


    CALL finalize_windows



END PROGRAM test_windows