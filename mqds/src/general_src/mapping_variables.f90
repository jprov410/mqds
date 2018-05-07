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
        USE input_output
        USE parameters
        USE random_numbers
        IMPLICIT NONE
        INTEGER :: i, j
        REAL(dp) :: r( nstate ), theta( nstate ), A_f, A_b
        REAL(dp) :: rt( nstate ), thetat( nstate )
        REAL(dp), INTENT(out) :: x( nstate ), p( nstate )
        REAL(dp), INTENT(out) :: xt( nstate ), pt( nstate )
        COMPLEX(dp) :: c_f( nstate ), c_b( nstate )
        !COMPLEX(dp) :: xi_f( nstate ), xi_b( nstate )
        COMPLEX(dp) :: prod_f, prod_b
        COMPLEX(dp) :: xi_f, xi_b

        xi_f = 0.0_dp ; xi_b = 0.0_dp ; A_f = 0.0_dp
        r = 0.0_dp ; rt = 0.0_dp ; A_b = 0.0_dp
        theta = 0.0_dp ; thetat = 0.0_dp
        prod_f = 1.0_dp ; prod_b = 1.0_dp


        ! Sample angle variables
        theta = 2.0_dp * pi * uniform_rn( theta )
        thetat = 2.0_dp * pi * uniform_rn( thetat )
        !print*, theta, thetat

        ! change this back it didnt work!
        !!!!below

        A_b = 1.0_dp
        A_f = 1.0_dp
        !A_f = SQRT( DOT_PRODUCT( c_f, CONJG(c_f) ) )
        !A_b = SQRT( DOT_PRODUCT( c_b, CONJG(c_b) ) )

        !!! Maybe devided by radius???
        DO i = 1, nstate
        !    A_f = A_f + SQRT( c_f(i) * CONJG( c_f(i) ) )
        !    A_b = A_f + SQRT( c_b(i) * CONJG( c_b(i) ) )
        !!
        !
            xi_f = xi_f + c_f(i) * EXP( -eye * theta(i) )
            xi_b = xi_b + c_b(i) * EXP( eye * thetat(i) )
            prod_f = prod_f * radial_f(i, SIZE(radial_f(i,:)) - 1)
            prod_b = prod_b * radial_b(i, SIZE(radial_b(i,:)) - 1)
        END DO

        !weight_f =  0.5_dp * SQRT(pi) * xi_f
        !weight_b =  0.5_dp * SQRT(pi) * xi_b

        weight_f = A_f * prod_f * SQRT( pi ) * 0.5_dp * xi_f
        weight_b = A_b * prod_b * SQRT( pi ) * 0.5_dp * xi_b
        !print*, A_f, prod_f, xi_f
        !print*, A_b, prod_b, xi_b
        !print*, weight_f, weight_b

        prod = weight_f * weight_b * SQRT( 0.5_dp * pi )

        CALL sample_current_cdfs( r, rt )

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
