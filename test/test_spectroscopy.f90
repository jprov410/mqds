!> Tests the spectroscopy module.
PROGRAM test_spectroscopy
    USE kinds
    USE spectroscopy
    USE input_output
    USE mapping_variables
    IMPLICIT NONE
    INTEGER :: k1, k2, k3
    REAL(dp) :: x(2), p(2), xt(2), pt(2)
    COMPLEX(dp) :: test(2,2), check(2,2), w, resp(8) = 0.0_dp
    COMPLEX(dp) :: trace_test(2,2,0:1), trace_check(0:1)
    nstate = 2 ; nbstep = 1 ; dump = 1 ; prod = 1.0_dp
    test = 0.0_dp ; test(1,1) = 1.0_dp
    check = 0.0_dp ; check(2,1) = 1.0_dp
    trace_check = 0.0_dp ; trace_test = 0.0_dp
    x = 0.0_dp ; p = 0.0_dp ; xt = 0.0_dp ; pt = 0.0_dp
    w = 1.0_dp ; k1 = 1 ; k2 = 1 ; k3 = 1

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

    initstate = 1 ; initstatet = 1
    x(2) = 1.0_dp ; p(2) = 1.0_dp
    xt(2) = 1.0_dp ; pt(2) = 1.0_dp
    CALL dipole_focus(x, p, xt, pt, w, k1)
    IF (initstate == 1 .AND. &
            initstatet == 1) THEN
        WRITE(*,*) 'Dipole focus is not working correctly!'
        STOP 1
    END IF

    x = 1.0_dp ; p = 1.0_dp
    xt = 1.0_dp ; pt = 1.0_dp
    k1=-1 ; k2=+1 ; k3=+1 ; w = 1.0_dp
    resp = nonlinear_response(x,p,xt,pt,w,k1,k2,k3)
    IF ( resp( 1 ) == (0.0_dp,0.0_dp) ) THEN
        WRITE(*,*) 'Nonlinear response is not working correctly!'
        STOP 1
    END IF

    k1=+1 ; k2=-1 ; k3=-1 ; w = 1.0_dp
    resp = nonlinear_response(x,p,xt,pt,w,k1,k2,k3)
    IF ( resp( 2 ) == (0.0_dp,0.0_dp) ) THEN
        WRITE(*,*) 'Nonlinear response is not working correctly!'
        STOP 1
    END IF

    k1=+1 ; k2=-1 ; k3=+1 ; w = 1.0_dp
    resp = nonlinear_response(x,p,xt,pt,w,k1,k2,k3)
    IF ( resp( 3 ) == (0.0_dp,0.0_dp) ) THEN
        WRITE(*,*) 'Nonlinear response is not working correctly!'
        STOP 1
    END IF

    k1=-1 ; k2=+1 ; k3=-1 ; w = 1.0_dp
    resp = nonlinear_response(x,p,xt,pt,w,k1,k2,k3)
    IF ( resp( 4 ) == (0.0_dp,0.0_dp) ) THEN
        WRITE(*,*) 'Nonlinear response is not working correctly!'
        STOP 1
    END IF

    k1=+1 ; k2=+1 ; k3=-1 ; w = 1.0_dp
    resp = nonlinear_response(x,p,xt,pt,w,k1,k2,k3)
    IF ( resp( 5 ) == (0.0_dp,0.0_dp) ) THEN
        WRITE(*,*) 'Nonlinear response is not working correctly!'
        STOP 1
    END IF

    k1=-1 ; k2=-1 ; k3=+1 ; w = 1.0_dp
    resp = nonlinear_response(x,p,xt,pt,w,k1,k2,k3)
    IF ( resp( 6 ) == (0.0_dp,0.0_dp) ) THEN
        WRITE(*,*) 'Nonlinear response is not working correctly!'
        STOP 1
    END IF

    k1=+1 ; k2=+1 ; k3=+1 ; w = 1.0_dp
    resp = nonlinear_response(x,p,xt,pt,w,k1,k2,k3)
    IF ( resp( 7 ) == (0.0_dp,0.0_dp) ) THEN
        WRITE(*,*) 'Nonlinear response is not working correctly!'
        STOP 1
    END IF

    k1=-1 ; k2=-1 ; k3=-1 ; w = 1.0_dp
    resp = nonlinear_response(x,p,xt,pt,w,k1,k2,k3)
    IF ( resp( 8 ) == (0.0_dp,0.0_dp) ) THEN
        WRITE(*,*) 'Nonlinear response is not working correctly!'
        STOP 1
    END IF

    CALL finalize_spectroscopy

END PROGRAM test_spectroscopy