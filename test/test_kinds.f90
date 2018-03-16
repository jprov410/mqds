!> Tests the defined single- and double-precision kinds.
PROGRAM test_kinds
    USE kinds
    IMPLICIT NONE
    REAL(dp) :: testdp = 1.0_dp
    REAL(sp) :: testsp = 1.0_sp

    IF( KIND(testdp) /= 8 .OR. KIND(testsp) /= 4 ) THEN
        WRITE(*,*) 'The defined kinds are not working correctly!'
        STOP 1
    END IF

END PROGRAM test_kinds