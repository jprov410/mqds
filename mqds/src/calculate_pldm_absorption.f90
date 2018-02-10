!
! This subroutine calculates the reduced density matrix
! using PLDM with a harmonic bath model in the diabatic basis
! with bilinear coupling
!
SUBROUTINE calculate_pldm_absorption
  USE kinds
  USE spectroscopy
  USE unit_conversions
  USE hamiltonians
  USE mapping_variables
  USE harmonic_bath
  USE input_output
  IMPLICIT NONE
  INTEGER :: i, j, k
  INTEGER :: istep, itraj, itime
  INTEGER :: istate, istatet
  REAL(dp) :: beta
  REAL(dp) :: weight
  REAL(dp) :: ham(nstate, nstate)
  REAL(dp) :: bath_force( nosc * nbath )
  REAL(dp) :: dt_bath, dt_map, printstep
  COMPLEX(dp) :: redmat(nstate, nstate, 0 : nbstep / dump)
  COMPLEX(dp) :: resp_func(0 : nbstep / dump)
  COMPLEX(dp) :: dipcom_gstate(nstate, nstate) ! dictates initial states for prop
  ham = 0.0_dp ; bath_force = 0.0_dp ; dt_bath = 0.0_dp ; dt_map = 0.0_dp
  redmat = ( 0.0_dp, 0.0_dp ) ; beta = 0.0_dp ; resp_func = ( 0.0_dp, 0.0_dp )
  dipcom_gstate = (0.0_dp, 0.0_dp) ; dipcom_gstate(1,1) = (1.0_dp, 0.0_dp)
  weight = 0.0_dp

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

  ! Allocate arrays for spectroscopy calculation
  CALL initialize_spectroscopy

  ! Get initially occupied coherences for calculation
  ! of linear spectroscopy assuming originally in
  ! system state (1,1) (initialized above)
  dipcom_gstate = dipole_commutator( dipcom_gstate )

  DO istate=1, nstate
      DO istatet=1, nstate
          ! check if dynamics from this state are necessary
          IF ( dipcom_gstate(istate, istatet) /= 0.0_dp ) THEN

              ! Set initial staties
              initstate = istate ; initstatet = istatet

              ! Set weight from transition dipole
              weight = REAL(dipcom_gstate(initstate, initstatet))

              ! Do the Dynamics
              DO itraj=1, ntraj
                  itime=0

                  ! Sample the initial conditions for system mapping and bath DOFs
                  CALL sample_thermal_wigner(x_bath, p_bath, beta)
                  CALL sample_pldm_map(x_map, p_map, xt_map, pt_map)

                  ! Calculate the t=0 redmat weighted by initial dipole
                  redmat(:, :, itime) = redmat(:, :, itime) &
                          + weight * pldm_redmat(x_map, p_map, xt_map, pt_map)

                  ! Find the initial bath force
                  bath_force = bilinear_harmonic_force_pldm(x_bath, x_map, p_map, xt_map, pt_map)

                  ! Propagate
                  DO istep=1, nbstep

                      ! First half of the verlet
                      x_bath = x_bath + p_bath * dt_bath + bath_force * 0.5_dp * dt_bath ** 2
                      p_bath = p_bath + bath_force * 0.5_dp * dt_bath


                      ! Update the full hamiltonian
                      ham = diabatic_bilinear_coupling_hamiltonian(x_bath, c)

                      ! Propagate the mapping variable
                      CALL verlet_mapping_variables(x_map, p_map, ham, dt_map)
                      CALL verlet_mapping_variables(xt_map, pt_map, ham, dt_map)

                      ! If the step is divisible by dump, compute the redmat weighted by initial dipolt
                      IF ( MOD( istep, dump) == 0 ) THEN
                          itime = itime + 1
                          redmat(:, :, itime) = redmat(:, :, itime) &
                                  + weight * pldm_redmat(x_map, p_map, xt_map, pt_map)
                      END IF

                      ! Update the force and finish the verlet
                      bath_force = bilinear_harmonic_force_pldm(x_bath, x_map, p_map, xt_map, pt_map)
                      p_bath = p_bath + bath_force * 0.5_dp * dt_bath

                  END DO

              END DO

          END IF

      END DO

  END DO


  ! Normalize the perturbed density matrix
  redmat = redmat / REAL( ntraj )

  ! Expectation value of the dipole operator with perturbed density matrix
  DO itime=0, ( nbstep / dump )
      redmat(:, :, itime) = dipole_operator( redmat(:, :, itime) )
  END DO
  ! calculate trace at all times
  resp_func = linear_response( redmat )

  ! Write the response function
  CALL write_linear_response(method, resp_func, printstep)

  ! Deallocate spectroscopy arrays
  CALL finalize_spectroscopy

END SUBROUTINE calculate_pldm_absorption
