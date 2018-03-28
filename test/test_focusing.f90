!> Tests the focusing procedure.
PROGRAM test_focusing
    USE kinds
    USE focusing
    USE input_output
    IMPLICIT NONE
    REAL(dp) :: x(2), p(2), xt(2), pt(2)
    COMPLEX(dp) :: test(2,2), test_weight
    nstate = 2 ; test(:,:) = 0.0_dp ; test(1,1) = 1.0_dp
    test_weight = 1.0_dp

    ! Try to focus
    initstate = 1 ; initstatet = 1
    CALL focus_pldm_map(x, p, xt, pt)
    ! Check if the focusing worked properly
    IF ( x(1) /= 1.0_dp .OR. x(2) /= 0.0_dp) THEN
        WRITE(*,*) 'focus_pldm_map subroutine did not work properly!'
        STOP 1
    END IF

    initstate = 2 ; initstatet=2
    CALL focus_redmat(test, test_weight)
    IF ( initstate /= 1 .OR. initstatet /= 1 .OR. test_weight /= 1.0_dp ) THEN
        WRITE(*,*) 'focus_redmat subroutine did not work properly!'
        STOP 1
    END IF

END PROGRAM test_focusing