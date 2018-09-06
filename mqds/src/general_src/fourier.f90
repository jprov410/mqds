MODULE fourier
    USE kinds
    USE parameters
CONTAINS
    FUNCTION ft_at_omega(f,w,dt) RESULT(res)
        IMPLICIT NONE
        INTEGER :: itime
        REAL(dp), INTENT(in) :: f(:)
        REAL(dp), INTENT(in) :: w, dt
        REAL(dp) :: res
        res = 0.0_dp

        DO itime = 1, SIZE(f)
            res = res + dt * f(itime) * exp(- eye * w * (itime-1) * dt)
        END DO
    END FUNCTION ft_at_omega

END MODULE fourier