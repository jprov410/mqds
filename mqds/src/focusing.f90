!
! Module for the monte carlo "focusing" procedure defined in
!
MODULE focusing
    USE kinds
    USE input_output
    USE parameters
    USE mapping_variables
CONTAINS

    SUBROUTINE focus_pldm_map(x_init, p_init, xt_init, pt_init)
        IMPLICIT NONE
        REAL(dp), INTENT(out) :: x_init(nstate), p_init(nstate)
        REAL(dp), INTENT(out) :: xt_init(nstate), pt_init(nstate)
        x_init = 0.0_dp ; p_init = 0.0_dp
        xt_init = 0.0_dp ; pt_init = 0.0_dp

        prod = (0.0_dp, 0.0_dp)

        x_init(initstate) = 1.0_dp
        p_init(initstate) = 1.0_dp
        xt_init(initstatet) = 1.0_dp
        pt_init(initstatet) = 1.0_dp

        prod = 0.5_dp * ( x_init(initstate) - eye * p_init(initstate) ) &
                * ( xt_init(initstatet) + eye * pt_init(initstatet) )
    END SUBROUTINE focus_pldm_map

END MODULE focusing