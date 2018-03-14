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

    SUBROUTINE focus_redmat( input_redmat, weight )
        USE random_numbers
        IMPLICIT NONE
        INTEGER :: istate, istatet
        COMPLEX(dp), INTENT(in) :: input_redmat(nstate,nstate)
        COMPLEX(dp), INTENT(inout) :: weight
        COMPLEX(dp) :: phase
        REAL(dp) :: unif(1), norm, cdf

        cdf = 0.0_dp ; unif = 0.0_dp ; norm = 0.0_dp

        ! Pull uniform random number to sample CDF
        unif = uniform_rn(unif)

        ! Build normalization factor
        DO istate = 1, nstate
            DO istatet = 1, nstate
                norm = norm + SQRT(input_redmat(istate,istatet) &
                        * CONJG(input_redmat(istate,istatet)))
            END DO
        END DO


        ! Identify the state redmat element to focus on by randomly
        ! sampling the integrated pdf
        DO istate = 1, nstate
            DO istatet = 1, nstate
                cdf = cdf + SQRT(input_redmat(istate,istatet) &
                        * CONJG(input_redmat(istate,istatet))) &
                        / norm
                IF ( cdf >= unif(1) ) THEN
                    initstate = istate
                    initstatet = istatet
                    GO TO 55
                END IF
            END DO
        END DO

        55 CONTINUE

        phase = input_redmat(initstate,initstatet) &
                / SQRT(input_redmat(initstate,initstatet) &
                        * CONJG(input_redmat(initstate,initstatet)))

        weight = weight * norm * phase

    END SUBROUTINE focus_redmat

END MODULE focusing