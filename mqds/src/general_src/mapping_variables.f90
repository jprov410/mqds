!> This module contains the necessary functions and
!! subroutines to use the Meyer-Miller mapping model
!! which maps \f$|\alpha\rangle \rightarrow \hat{a}_\alpha^\dagger =
!! \frac{1}{\sqrt{2}}(\hat{x}_\alpha - i\hat{p}_\alpha) \f$
MODULE mapping_variables
  USE kinds
  IMPLICIT NONE
  REAL(dp), ALLOCATABLE :: x_map(:), p_map(:)
  REAL(dp), ALLOCATABLE :: xt_map(:), pt_map(:)
  COMPLEX(dp) :: prod, weight_f, weight_b

CONTAINS

    !----------------allocate and deallocate mapping variables----------------!

    !> Allocate the mapping variable arrays for PLDM calculation
    SUBROUTINE initialize_pldm_map
        USE input_output
        IMPLICIT NONE

        ALLOCATE(x_map(nstate),xt_map(nstate),p_map(nstate),pt_map(nstate))

    END SUBROUTINE initialize_pldm_map
  
    !> Deallocate the mapping variable arrays for PLDM calculation
    SUBROUTINE finalize_pldm_map
    IMPLICIT NONE
    
    DEALLOCATE( x_map, xt_map, p_map, pt_map )
    
  END SUBROUTINE finalize_pldm_map

    !> Allocate the mapping variable arrays for TWA or SQC calculation
    SUBROUTINE initialize_twa_map
    USE input_output
    IMPLICIT NONE
    
    ALLOCATE( x_map(nstate), p_map(nstate) )
    
  END SUBROUTINE initialize_twa_map
  
    !> Deallocate the mapping variable arrays for TWA or SQC calculation
    SUBROUTINE finalize_twa_map
    IMPLICIT NONE

    DEALLOCATE( x_map, p_map )
    
  END SUBROUTINE finalize_twa_map

    !----------sample the initial distributions of mapping variables----------!

    !> Sample the initial distribution of mapping variables for PLDM from
    !! gaussian distributions
    SUBROUTINE sample_pldm_map(x_init, p_init, xt_init, pt_init)
        USE random_numbers
        USE parameters
        USE input_output
        USE kinds
        IMPLICIT NONE
        REAL(dp), INTENT(out) :: x_init(nstate), p_init(nstate)
        REAL(dp), INTENT(out) :: xt_init(nstate), pt_init(nstate)

        prod = (0.0_dp, 0.0_dp)

        x_init(:) = gaussian_rn(x_init)
        p_init(:) = gaussian_rn(p_init)
        xt_init(:) = gaussian_rn(xt_init)
        pt_init(:) = gaussian_rn(pt_init)

        prod = 0.5_dp * ( x_init(initstate) - eye * p_init(initstate) ) &
                * ( xt_init(initstatet) + eye * pt_init(initstatet) )
    END SUBROUTINE sample_pldm_map

    !> Subroutine that takes the current density matrix and uses it
    !! for coefficients in a linear combination of new initial Hermite
    !! polynomials

    SUBROUTINE pldm_map_hop( c_f, c_b, x, p, xt, pt )
    !SUBROUTINE pldm_map_hop( c, x, p, xt, pt, islice )
        USE input_output
        USE parameters
        USE random_numbers
        IMPLICIT NONE
        INTEGER :: i, j, islice
        REAL(dp) :: r( nstate ), theta( nstate )
        REAL(dp) :: rt( nstate ), thetat( nstate )
        !REAL(dp) :: norm_f, norm_b
        REAL(dp) :: norm
        REAL(dp), INTENT(out) :: x( nstate ), p( nstate )
        REAL(dp), INTENT(out) :: xt( nstate ), pt( nstate )
        !COMPLEX(dp), INTENT(in) :: c( nstate, nstate )
        COMPLEX(dp),INTENT(in) :: c_f( nstate ), c_b( nstate )
        !COMPLEX(dp) :: xi, A
        COMPLEX(dp) :: xi_f, xi_b ! fwd and bkwd lin comb of phases
        COMPLEX(dp) :: P_f, P_b ! sampled probability distributions
        REAL(dp) :: N_f, N_b ! normalizations

        xi_f = 0.0_dp ; xi_b = 0.0_dp ; N_f = 0.0_dp
        r = 0.0_dp ; rt = 0.0_dp ; N_b = 0.0_dp
        theta = 0.0_dp ; thetat = 0.0_dp !; xi = 0.0_dp ; A = 0.0_dp
        N_f = 0.0_dp ; N_b = 0.0_dp
        P_f = 0.0_dp ; P_b = 0.0_dp

        !IF ( islice == 1 ) THEN
        !    norm = 1.0_dp
        !    A = 1.0_dp
        !END IF

        ! Sample angle variables
        theta = 2.0_dp * pi * uniform_rn( theta )
        thetat = 2.0_dp * pi * uniform_rn( thetat )

        CALL sample_current_cdfs( r, rt )

        DO i = 1, nstate
            !DO j = 1, nstate
            !    xi = xi + r(i) * rt(j) * c(i,j) * EXP( -eye * (theta(i) - thetat(j) ) )
            !    !A = A + r(i) * rt(j) * SQRT( c(i,j) * CONJG( c(i,j) ) )
                !norm = norm + SQRT( c(i,j) * CONJG( c(i,j) ) )
            !END DO
            !j=i
            !A = A + r(i) * rt(j) * SQRT( c(i,j) * CONJG( c(i,j) ) )
            !norm = norm + SQRT( c(i,j) * CONJG( c(i,j) ) )

            ! Inintial forward boundary terms (xi_f/A_f)
            xi_f = xi_f + r(i) * c_f(i) * EXP( -eye * theta(i) )
            P_f = P_f + SQRT( c_f(i) * CONJG( c_f(i) ) ) * r(i)
            N_f = N_f + SQRT( c_f(i) * CONJG( c_f(i) ) )

            !A_f = A_f + SQRT(SQRT( c(i,i) * CONJG( c(i,i) ) ) ) * r(i)
            !norm_f = norm_f + SQRT(SQRT( c(i,i) * CONJG( c(i,i) ) ))

            !A_b = A_b + SQRT(SQRT( c(i,i) * CONJG( c(i,i) ) ) ) * rt(i)
            !norm_b = norm_b + SQRT(SQRT( c(i,i) * CONJG( c(i,i) ) ))

            ! Inintial backward boundary term (xi_b/A_b)
            xi_b = xi_b + rt(i) * c_b(i) * EXP( eye * thetat(i) )
            P_b = P_b + SQRT( c_b(i) * CONJG( c_b(i) ) ) * rt(i)
            N_b = N_b + SQRT( c_b(i) * CONJG( c_b(i) ) )
        END DO
        N_f = N_f * DSQRT(0.5_dp * pi)
        N_b = N_b * DSQRT(0.5_dp * pi)

        weight_f = DSQRT(0.5_dp) * xi_f / ( P_f / N_f )
        weight_b = DSQRT(0.5_dp) * xi_b / ( P_b / N_b )

        prod = weight_f * weight_b
        !print*, norm_f, norm_b

        !norm = norm * norm_f * norm_b
        !norm = norm * pi / 2.0_dp
        !print*, 'norm =', norm

        !A = A_f * A_b
        !A = A * A_f * A_b
        !prod = 0.5_dp * norm * xi / A
        !prod = 0.5_dp *  xi * (norm / A)!(A_f * A_b)



        x( : ) = r( : ) * DCOS( theta )
        xt( : ) = rt( : ) * DCOS( thetat )
        p( : ) = r( : ) * DSIN( theta )
        pt( : ) = rt( : ) * DSIN( thetat )

    END SUBROUTINE pldm_map_hop

    !> Sample the initial distribution of mapping variables for
    !! truncated wigner approximation in action-angle variables
    SUBROUTINE sample_twa_map(x_init, p_init)
        USE random_numbers
        USE parameters
        USE input_output
        USE kinds
        IMPLICIT NONE
        REAL(dp), INTENT(out) :: x_init(nstate), p_init(nstate)
        REAL(dp) :: n(nstate), q(nstate)
        prod = (0.0_dp, 0.0_dp)

        p_init = gaussian_rn(p_init) / 2.0_dp
        x_init = gaussian_rn(x_init) / 2.0_dp

        prod = (2.0_dp)**(nstate + 1) &
                * ( x_init(initstate) - eye * p_init(initstate) ) &
                * ( x_init(initstatet) + eye * p_init(initstatet) )

        IF (initstate == initstatet) prod = prod - 0.5_dp * (2.0_dp)**(nstate + 1)

    END SUBROUTINE sample_twa_map

    !-----------------calculate the reduced density matrix------------------!

    !> Calculate the reduced density matrix using PLDM with the mapping variables
    !! from the product of "Hermite polynomials"
    FUNCTION pldm_redmat(x, p, xt, pt) RESULT( res )
        USE kinds
        USE parameters
        USE input_output
        IMPLICIT NONE
        INTEGER :: i, j
        REAL(dp), INTENT(in) :: x(nstate), p(nstate), xt(nstate), pt(nstate)
        COMPLEX(dp) :: res( nstate, nstate )

        DO i=1, nstate
            DO j=1, nstate
                res(j,i) = 0.5_dp * ( x(j) + eye * p(j) ) * ( xt(i) - eye * pt(i) ) * prod
            END DO
        END DO

    END FUNCTION pldm_redmat

    !> Calculate the reduced density matrix using TWA with the mapping variables
    FUNCTION twa_redmat(x, p) RESULT( res )
        USE kinds
        USE parameters
        USE input_output
        IMPLICIT NONE
        INTEGER :: i, j
        REAL(dp), INTENT(in) :: x(nstate), p(nstate)
        COMPLEX(dp) :: res( nstate, nstate )
        res = 1.0_dp

        !DO i=1, nstate
        !    DO j=1, nstate
        !        res(j,i) = 0.5_dp * ( x(j) + eye * p(j) ) * ( x(i) - eye * p(i) ) * prod
        !        IF (i == j) res(j,i) = res(j,i) - 0.5_dp * prod
        !    END DO
        !END DO

    END FUNCTION twa_redmat

    !-------------------propagators for mapping variables-------------------!

    !> propagates x and p mapping variables using the Hamiltonian by
    !! integrating Hamilton's equations with the Verlet integrator.
    SUBROUTINE verlet_mapping_variables(x, p, H, dt)
        USE kinds
        USE input_output
        IMPLICIT NONE
        INTEGER :: istep, i, j
        REAL(dp), INTENT(inout) :: x(nstate), p(nstate)
        REAL(dp), INTENT(in) :: H(nstate, nstate)
        REAL(dp), INTENT(in) :: dt
        REAL(dp) :: dxdt(nstate), d2xdt2(nstate)
        REAL(dp) :: dpdt(nstate)

        DO istep=1, nlit

            !Calculate initial dpdt, dxdt, d2xdt2
            dpdt = 0.0_dp
            DO i=1, nstate
                DO j=1, nstate
                    dpdt(i) = dpdt(i) - x(j) * H(j,i)
                END DO
            END DO

            dxdt = 0.0_dp ; d2xdt2 = 0.0_dp
            DO i=1, nstate
                DO j=1, nstate
                    dxdt(i) = dxdt(i) + p(j) * H(j,i)
                    d2xdt2(i) = d2xdt2(i) + dpdt(j) * H(j,i)
                END DO
            END DO


            ! First half of verlet
            p = p + 0.5_dp * dpdt * dt
            x = x + dxdt * dt + d2xdt2 * 0.5_dp * dt ** 2

            ! Calculate updated dpdt
            dpdt = 0.0_dp
            DO i=1, nstate
                DO j=1, nstate
                    dpdt(i) = dpdt(i) - x(j) * H(j,i)
                END DO
            END DO

            ! Second half of verlet
            p = p + 0.5_dp * dpdt * dt

        END DO

    END SUBROUTINE verlet_mapping_variables

END MODULE mapping_variables
