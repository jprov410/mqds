PROGRAM test_spectroscopy
    USE kinds
    USE spectroscopy
    USE input_output
    IMPLICIT NONE
    COMPLEX(dp) :: test(2,2), check(2,2)
    COMPLEX(dp) :: trace_test(2,2,0:1), trace_check(0:1)
    nstate = 2 ; nbstep = 1 ; dump = 1
    test = 0.0_dp ; test(1,1) = 1.0_dp
    check = 0.0_dp ; check(2,1) = 1.0_dp
    trace_check = 0.0_dp ; trace_test = 0.0_dp

    CALL initialize_spectroscopy

    test = dipole_operator( test )
    IF ( test(1,1) /= check(1,1) .OR. &
            test(1,2) /= check(1,2) .OR. &
            test(2,1) /= check(2,1) .OR. &
            test(2,2) /= check(2,2) ) THEN
        WRITE(*,*) 'Dipole operator not working correctly!'
        STOP 1
    END IF

    test = 0.0_dp ; test(1,1) = 1.0_dp
    check = 0.0_dp ; check(2,1) = 1.0_dp ; check(1,2) = -1.0_dp

    test = dipole_commutator( test )
    IF ( test(1,1) /= check(1,1) .OR. &
            test(1,2) /= check(1,2) .OR. &
            test(2,1) /= check(2,1) .OR. &
            test(2,2) /= check(2,2) ) THEN
        WRITE(*,*) 'Dipole commutator not working correctly!'
        STOP 1
    END IF

    trace_test(1,1,0) = 1.0_dp ; trace_test(2,2,0) = 1.0_dp
    trace_check = system_trace( trace_test )
    IF ( trace_check(0) /= 2.0_dp ) THEN
        WRITE(*,*) 'System trace is not working correctly!'
        STOP 1
    END IF

    CALL finalize_spectroscopy

END PROGRAM test_spectroscopy