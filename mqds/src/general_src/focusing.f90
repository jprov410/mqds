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

    SUBROUTINE pldm_focus_manifold_map(x_init, p_init, xt_init, pt_init)
        USE random_numbers
        USE parameters
        USE input_output
        USE kinds
        IMPLICIT NONE
        REAL(dp), INTENT(out) :: x_init(nstate), p_init(nstate)
        REAL(dp), INTENT(out) :: xt_init(nstate), pt_init(nstate)

        prod = (0.0_dp, 0.0_dp)

        IF ( calculation == 'absorption' ) THEN
            IF ( initstate == 1 ) THEN
                x_init(:) = 0.0_dp ; p_init(:) = 0.0_dp
                x_init(1) = 1.0_dp ; p_init(1) = 1.0_dp
                xt_init(:) = gaussian_rn(xt_init)
                pt_init(:) = gaussian_rn(pt_init)
            ELSE IF ( initstatet == 1 ) THEN
                xt_init(:) = 0.0_dp ; pt_init(:) = 0.0_dp
                xt_init(1) = 1.0_dp ; pt_init(1) = 1.0_dp
                x_init(:) = gaussian_rn(x_init)
                p_init(:) = gaussian_rn(p_init)
            ELSE
                OPEN( UNIT=10, FILE=ERRORLOG )
                WRITE(10, *) 'There has been an error in the initialization of mapping &
                        variables for absorption spectrum calculation, neither state is the ground state'
                CLOSE(10)
            END IF
        END IF

        prod = 0.5_dp * ( x_init(initstate) - eye * p_init(initstate) ) &
                * ( xt_init(initstatet) + eye * pt_init(initstatet) )

    END SUBROUTINE pldm_focus_manifold_map

END MODULE focusing