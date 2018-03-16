!> Tests the focusing procedure.
PROGRAM test_focusing
    USE kinds
    USE focusing
    USE input_output
    IMPLICIT NONE
    REAL(dp) :: x(2), p(2), xt(2), pt(2)
    nstate = 2

    ! Try to focus
    initstate = 1 ; initstatet = 1
    CALL focus_pldm_map(x, p, xt, pt)
    ! Check if the focusing worked properly
    IF ( x(1) /= 1.0_dp .OR. x(2) /= 0.0_dp) THEN
        WRITE(*,*) 'focus_pldm_map subroutine did not work properly!'
        STOP 1
    END IF

END PROGRAM test_focusing