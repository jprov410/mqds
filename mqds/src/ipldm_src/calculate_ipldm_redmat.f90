!
! This subroutine calculates the reduced density matrix
! using ITERATIVE PLDM with a harmonic bath model in
! the diabatic basis with bilinear coupling
!
SUBROUTINE calculate_ipldm_redmat
  USE kinds
  USE unit_conversions
  USE hamiltonians
  USE mapping_variables
  USE harmonic_bath
  USE input_output
  USE focusing
  IMPLICIT NONE
  INTEGER :: i, j, k
  INTEGER :: istep, itraj, itime, islice
  INTEGER :: initstate_save, initstatet_save
  REAL(dp) :: beta
  REAL(dp) :: ham( nstate, nstate )
  REAL(dp) :: bath_force( nosc * nbath )
  REAL(dp) :: dt_bath, dt_map, printstep
  COMPLEX(dp) :: redmat( nstate, nstate, 0 : nbstep / dump )
  COMPLEX(dp) :: focusing_weight, focusing_redmat(nstate,nstate)
  ham = 0.0_dp ; bath_force = 0.0_dp ; dt_bath = 0.0_dp ; dt_map = 0.0_dp
  redmat = (0.0_dp, 0.0_dp) ; beta = 0.0_dp ; focusing_redmat = (0.0_dp, 0.0_dp)
  initstate_save = initstate ; initstatet_save = initstatet

  ! Calculate Beta for thermal sampling
  beta = 1.0_dp / ( temperature * convert('kelvin','au_energy') )

  ! Setup the necessary timings
  printstep = runtime * REAL( dump ) / REAL( nbstep )
  runtime = runtime * convert('fs','au_time')  
  dt_bath = runtime / REAL(nbstep)
  dt_map = dt_bath / REAL(nlit)

  ! Read the input electonic Hamiltonian
  CALL read_hel
  hel = hel * convert('wvnbr','au_energy')

  ! Make sure nbstep is integer multiple of nslice
  IF ( MOD( nbstep, nslice ) /= 0 ) THEN
      OPEN(UNIT=10, FILE=ERRORLOG)
      WRITE(10,*) 'nbstep must be an integer multiple of nslice for ipldm calculations'
      STOP 1
  END IF

  DO itraj=1, ntraj
      initstate = initstate_save ; initstatet = initstatet_save
      focusing_weight = 1.0_dp
      itime=0
     
     ! Sample the initial conditions for system mapping and bath DOFs
     CALL sample_thermal_wigner(x_bath, p_bath, beta)
     CALL sample_pldm_map(x_map, p_map, xt_map, pt_map)

     ! Calculate the t=0 redmat
     redmat(:, :, itime) = redmat(:, :, itime) + &
             focusing_weight * pldm_redmat(x_map, p_map, xt_map, pt_map)

     ! Find the initial bath force
     bath_force = bilinear_harmonic_force_pldm(x_bath, x_map, p_map, xt_map, pt_map)
     DO islice = 1, nslice
         DO istep=1, INT(nbstep/nslice)

             ! First half of the verlet
             x_bath = x_bath + p_bath * dt_bath + bath_force * 0.5_dp * dt_bath ** 2
             p_bath = p_bath + bath_force * 0.5_dp * dt_bath


             ! Update the full hamiltonian
             ham = diabatic_bilinear_coupling_hamiltonian(x_bath, c)

             ! Propagate the mapping variable
             CALL verlet_mapping_variables(x_map, p_map, ham, dt_map)
             CALL verlet_mapping_variables(xt_map, pt_map, ham, dt_map)

             ! If the step is divisible by dump, compute the redmat
             IF ( MOD( istep, dump) == 0 ) THEN
                 itime = itime + 1
                 redmat(:, :, itime) = redmat(:, :, itime) + &
                         focusing_weight * pldm_redmat(x_map, p_map, xt_map, pt_map)
             END IF

             ! Update the force and finish the verlet
             bath_force = bilinear_harmonic_force_pldm(x_bath, x_map, p_map, xt_map, pt_map)
             p_bath = p_bath + bath_force * 0.5_dp * dt_bath
         END DO
         ! Importance sampling of redmat
         ! sets new initial states for subsequent focus_pldm_map call
         focusing_redmat = pldm_redmat(x_map, p_map, xt_map, pt_map)
         CALL focus_redmat( focusing_redmat, focusing_weight )
         ! "Steepest descent" of intermediate integrals
         CALL focus_pldm_map( x_map, p_map, xt_map, pt_map )

     END DO

  END DO

  redmat = redmat / REAL( ntraj )
  
  CALL write_redmat(method, redmat, printstep)


END SUBROUTINE calculate_ipldm_redmat
