!
! This module contains the necessary functions and
! subroutines to use the meyer-miller mapping model
! within PLDM and TWA
!
MODULE mapping_variables
  USE kinds
  IMPLICIT NONE
  REAL(dp), ALLOCATABLE :: x_map(:), p_map(:)
  REAL(dp), ALLOCATABLE :: xt_map(:), pt_map(:)
  COMPLEX(dp) :: prod
  
CONTAINS
  
  ! Allocate the mapping variable arrays for PLDM calculation
  SUBROUTINE initialize_pldm_map    
    USE input_output
    IMPLICIT NONE
    
    ALLOCATE(x_map(nstate),xt_map(nstate),p_map(nstate),pt_map(nstate))
    
  END SUBROUTINE initialize_pldm_map
  
  ! Deallocate the mapping variable arrays for PLDM calculation
  SUBROUTINE finalize_pldm_map
    IMPLICIT NONE
    
    DEALLOCATE(x_map,xt_map,p_map,pt_map)
    
  END SUBROUTINE finalize_pldm_map

  ! Allocate the mapping variable arrays for TWA calculation
  SUBROUTINE initialize_twa_map    
    USE input_output
    IMPLICIT NONE
    
    ALLOCATE( x_map(nstate) ,p_map(nstate) )
    
  END SUBROUTINE initialize_twa_map
  
  ! Deallocate the mapping variable arrays for TWA calculation
  SUBROUTINE finalize_twa_map
    IMPLICIT NONE
    
    DEALLOCATE( x_map, p_map )
    
  END SUBROUTINE finalize_twa_map
  
  ! Sample the initial distribution of mapping variables for PLDM
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

  ! Sample the initial distribution of mapping variables for 
  ! truncated wigner approximation in action-angle variables
  SUBROUTINE sample_twa_map(x_init, p_init)
    USE random_numbers
    USE parameters
    USE input_output
    USE kinds
    IMPLICIT NONE
    REAL(dp), INTENT(out) :: x_init(nstate), p_init(nstate)
    REAL(dp) :: n(nstate), q(nstate)
    n = 0.0_dp ; q = 0.0_dp
    prod = (0.0_dp, 0.0_dp)
    
    n(initstate) = n(initstate) + 0.5_dp
    n(initstatet) = n(initstatet) + 0.5_dp

    q(:) = 2.0_dp * pi * uniform_rn(q)
    
    prod = & !DSQRT( n( initstate ) + zpe ) * DSQRT( n( initstatet ) + zpe ) &
         EXP( eye * ( q(initstate) - q(initstatet) ) )

    p_init = -DSQRT( 2.0_dp * ( n + zpe ) ) * DSIN( q ) 
    x_init = DSQRT( 2.0_dp * ( n + zpe ) ) * DCOS( q ) 

  END SUBROUTINE sample_twa_map

    ! propagates x and p mapping variables using the Hamiltonian (H)
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

  ! Calculate the reduced density matrix using PLDM with the mapping variables
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


  ! Calculate the reduced density matrix using PLDM with the mapping variables
  FUNCTION twa_redmat(x, p) RESULT( res )
    USE kinds
    USE parameters
    USE input_output
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(dp), INTENT(in) :: x(nstate), p(nstate)
    REAL(dp) :: n(nstate), q(nstate)
    COMPLEX(dp) :: res( nstate, nstate )
    n = 0.5_dp * ( p**2 + x**2 ) - zpe
    q = -DATAN2( p, x )

    DO i=1, nstate
       DO j=1, nstate
          res(j,i) = &  !DSQRT( n(i) + zpe ) * DSQRT( n(j) + zpe ) &
               EXP( eye * ( q(j) - q(i) ) ) * prod
!          res(j,i) = 0.5_dp * ( x(j) + eye * p(j) ) * ( x(i) - eye * p(i) ) * prod
       END DO
    END DO

  END FUNCTION twa_redmat

END MODULE mapping_variables
