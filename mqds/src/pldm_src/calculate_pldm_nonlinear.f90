!> This subroutine calculates the reduced density matrix
!! using PLDM with a harmonic bath model in the diabatic basis
!! with bilinear coupling.
SUBROUTINE calculate_pldm_nonlinear
  USE kinds
  USE spectroscopy
  USE unit_conversions
  USE hamiltonians
  USE mapping_variables
  USE focusing
  USE parameters
  USE harmonic_bath
  USE input_output
  IMPLICIT NONE
  INTEGER :: i, j, k
  INTEGER :: istep, itraj
  INTEGER :: istate, istatet
  INTEGER :: itime1, itime2, itime3
  INTEGER :: tdim1, tdim2, tdim3
  INTEGER :: k1_sign
  INTEGER, ALLOCATABLE :: k2_sign( : ), k3_sign( :, : )
  INTEGER, ALLOCATABLE :: istate_save( :, : ), istatet_save( :, : )
  REAL(dp) :: beta
  REAL(dp) :: weight
  REAL(dp) :: ham(nstate, nstate)
  REAL(dp) :: bath_force( nosc * nbath )
  REAL(dp) :: dt_bath_1, dt_map_1
  REAL(dp) :: dt_bath_2, dt_map_2
  REAL(dp) :: dt_bath_3, dt_map_3
  REAL(dp) :: dt_bath, dt_map
  REAL(dp), ALLOCATABLE :: x_bath_save( :,:,: ), p_bath_save( :,:,: )
  COMPLEX(dp) :: dipcom_gstate(nstate, nstate) ! dictates initial states for prop
  COMPLEX(dp), ALLOCATABLE :: monte_carlo_weight( :, : )
  COMPLEX(dp), ALLOCATABLE :: resp_func( :, :, :, : )
  ham = 0.0_dp ; bath_force = 0.0_dp
  beta = 0.0_dp ; resp_func = ( 0.0_dp, 0.0_dp )
  dipcom_gstate = (0.0_dp, 0.0_dp) ; dipcom_gstate(1,1) = (1.0_dp, 0.0_dp)
  weight = 0.0_dp ; tdim1 = INT( nbstep1 / branch1 ) ; tdim2 = INT( nbstep2 / branch2 )
  tdim3 = INT( nbstep3 / branch3 )

  ! Setup timings
  tdelay1 = tdelay1 * convert('fs','au_time')
  tdelay2 = tdelay2 * convert('fs','au_time')
  tdelay3 = tdelay3 * convert('fs','au_time')
  dt_bath_1 = tdelay1 / DBLE( nbstep1 ) ; dt_map_1 = dt_bath_1 / DBLE( nlit )
  dt_bath_2 = tdelay2 / DBLE( nbstep2 ) ; dt_map_2 = dt_bath_2 / DBLE( nlit )
  dt_bath_3 = tdelay3 / DBLE( nbstep3 ) ; dt_map_3 = dt_bath_3 / DBLE( nlit )

  ! Calculate Beta for thermal sampling
  beta = 1.0_dp / ( temperature * convert('kelvin','au_energy') )

  ! Read the input electonic Hamiltonian
  CALL read_hel
  hel = hel * convert('wvnbr','au_energy')

  ! If the number of steps is 0, set time interval dimension bound to 0
  IF ( nbstep1 == 0 ) THEN
      tdim1 = 0
  END IF
  IF ( nbstep2 == 0 ) THEN
      tdim2 = 0
  END IF
  IF ( nbstep3 == 0 ) THEN
      tdim3 = 0
  END IF

  ALLOCATE( monte_carlo_weight( 0 : tdim1, 0 : tdim2 ), &
          resp_func( 8, 0 : tdim1, 0 : tdim2, 0 : tdim3 ), &
          x_bath_save( 0 : tdim1, 0 : tdim2, nosc * nbath ), &
          p_bath_save( 0 : tdim1, 0 : tdim2, nosc * nbath ), &
          istate_save( 0 : tdim1, 0 : tdim2 ), &
          istatet_save( 0 : tdim1, 0 : tdim2 ), &
          k2_sign( 0 : tdim1 ), &
          k3_sign( 0 : tdim1, 0 : tdim2 ) )

  monte_carlo_weight = 0.0_dp ; resp_func = 0.0_dp
  x_bath_save = 0.0_dp ; p_bath_save = 0.0_dp
  istate_save = 0 ; istatet_save = 0
  k1_sign = 0 ; k2_sign = 0 ; k3_sign = 0


  ! Allocate arrays for spectroscopy calculation
  CALL initialize_spectroscopy

  ! Get initially occupied coherences for calculation
  ! of linear spectroscopy assuming originally in
  ! system state (1,1)
  dipcom_gstate = dipole_commutator( dipcom_gstate )

  ! Setup the necessary timings
!  printstep = runtime * REAL( dump ) / REAL( nbstep )
!  runtime = runtime * convert('fs','au_time')
!  dt_bath = runtime / REAL(nbstep)
!  dt_map = dt_bath / REAL(nlit)

  DO istate=1, nstate
      DO istatet=1, nstate
          ! check if dynamics from this state are necessary
          IF ( dipcom_gstate(istate, istatet) /= 0.0_dp ) THEN
              ! Set initial states
              initstate = istate ; initstatet = istatet
              ! Set weight from transition dipole
              weight = REAL(dipcom_gstate(initstate, initstatet))

              ! Assign plus or minus k1 sign based on initial states
              IF ( istate > istatet ) THEN
                  k1_sign = 1
              ELSE
                  k1_sign = -1
              END IF

              ! Do the Dynamics
              DO itraj=1, ntraj

                  ! Set the monte carlo weight to be the dipole transition weight
                  monte_carlo_weight = weight

                  ! Sample the initial conditions for system mapping and bath DOFs
                  CALL sample_thermal_wigner(x_bath, p_bath, beta)
                  CALL focus_pldm_map(x_map, p_map, xt_map, pt_map)

                  ! Set initial time index to 0 and assign values to dt
                  itime1= 0 ; dt_map = dt_map_1 ; dt_bath = dt_bath_1

                  !!!!!!!!!DO FOCUSING AT ITIME = 0!!!!!!!!!!!
                  !! focusing part assigns new initstate(t) values
                  CALL dipole_focus( x_map, p_map, xt_map, pt_map,&
                          monte_carlo_weight( itime1, 0 ) , k2_sign( itime1 ) )
                  monte_carlo_weight( itime1, : ) = monte_carlo_weight( itime1, 0 )
                  istate_save( itime1, : ) = initstate
                  istatet_save( itime1, : ) = initstatet
                  x_bath_save( itime1, tdim2, : ) = x_bath( : )
                  p_bath_save( itime1, tdim2, : ) = p_bath( : )


                  ! Find the initial bath force
                  bath_force = bilinear_harmonic_force_pldm(x_bath, x_map, p_map, xt_map, pt_map)

                  ! Propagate
                  DO istep=1, nbstep1

                      ! First half of the verlet
                      x_bath = x_bath + p_bath * dt_bath + bath_force * 0.5_dp * dt_bath ** 2
                      p_bath = p_bath + bath_force * 0.5_dp * dt_bath

                      ! Update the full hamiltonian
                      ham = diabatic_bilinear_coupling_hamiltonian(x_bath, c)

                      ! Propagate the mapping variable
                      CALL verlet_mapping_variables(x_map, p_map, ham, dt_map)
                      CALL verlet_mapping_variables(xt_map, pt_map, ham, dt_map)

                      ! If the step is divisible by dump, compute the redmat weighted by initial dipolt
                      IF ( MOD( istep , branch1 ) == 0 ) THEN
                          itime1 = itime1 + 1
                          !!!!!!!DO FOCUSING AT ITIME > 0!!!!!!!!!!!!!
                          CALL dipole_focus( x_map, p_map, xt_map, pt_map,&
                                  monte_carlo_weight( itime1, 0 ) , k2_sign( itime1 ) )
                          monte_carlo_weight( itime1, : ) = monte_carlo_weight( itime1, 0 )
                          istate_save( itime1, : ) = initstate
                          istatet_save( itime1, : ) = initstatet
                          x_bath_save( itime1, tdim2, : ) = x_bath( : )
                          p_bath_save( itime1, tdim2, : ) = p_bath( : )
                      END IF

                      ! Update the force and finish the verlet
                      bath_force = bilinear_harmonic_force_pldm(x_bath, x_map, p_map, xt_map, pt_map)
                      p_bath = p_bath + bath_force * 0.5_dp * dt_bath

                  END DO

                  ! second time interval
                  dt_bath = dt_bath_2 ; dt_map = dt_map_2

                  DO itime1 = 0, tdim1
                      itime2 = 0

                      ! Reinitialize from the distribution at the end of t1 propagation
                      x_bath(:) = x_bath_save( itime1, tdim2, : )
                      p_bath(:) = p_bath_save( itime1, tdim2, : )
                      initstate = istate_save( itime1, tdim2 )
                      initstatet = istatet_save( itime1, tdim2 )
                      CALL focus_pldm_map( x_map, p_map, xt_map, pt_map )

                      !!!!!!!!!DO FOCUSING AT ITIME = 0!!!!!!!!!!!
                      !! focusing part assigns new initstate(t) values
                      CALL dipole_focus( x_map, p_map, xt_map, pt_map,&
                              monte_carlo_weight( itime1, itime2 ) , k3_sign( itime1, itime2 ) )
                      istate_save( itime1, itime2 ) = initstate
                      istatet_save( itime1, itime2 ) = initstatet
                      x_bath_save( itime1, itime2, : ) = x_bath( : )
                      p_bath_save( itime1, itime2, : ) = p_bath( : )

                      ! Find the initial bath force
                      bath_force = bilinear_harmonic_force_pldm(x_bath, x_map, p_map, xt_map, pt_map)

                      ! Propagate
                      DO istep=1, nbstep2

                          ! First half of the verlet
                          x_bath = x_bath + p_bath * dt_bath + bath_force * 0.5_dp * dt_bath ** 2
                          p_bath = p_bath + bath_force * 0.5_dp * dt_bath

                          ! Update the full hamiltonian
                          ham = diabatic_bilinear_coupling_hamiltonian(x_bath, c)

                          ! Propagate the mapping variable
                          CALL verlet_mapping_variables(x_map, p_map, ham, dt_map)
                          CALL verlet_mapping_variables(xt_map, pt_map, ham, dt_map)

                          ! If the step is divisible by dump, compute the redmat weighted by initial dipolt
                          IF ( MOD( istep , branch2 ) == 0 ) THEN
                              itime2 = itime2 + 1

                              !!!!!!!DO FOCUSING AT ITIME > 0!!!!!!!!!!!!!
                              !! focusing part assigns new initstate(t) values
                              CALL dipole_focus( x_map, p_map, xt_map, pt_map,&
                                      monte_carlo_weight( itime1, itime2 ), k3_sign( itime1, itime2 ) )
                              istate_save( itime1, itime2 ) = initstate
                              istatet_save( itime1, itime2 ) = initstatet
                              x_bath_save( itime1, itime2, : ) = x_bath( : )
                              p_bath_save( itime1, itime2, : ) = p_bath( : )

                          END IF

                          ! Update the force and finish the verlet
                          bath_force = bilinear_harmonic_force_pldm(x_bath, x_map, p_map, xt_map, pt_map)
                          p_bath = p_bath + bath_force * 0.5_dp * dt_bath

                      END DO

                  END DO

                  ! Third time interval
                  dt_bath = dt_bath_3 ; dt_map = dt_map_3

                  DO itime1 = 0, tdim1
                      DO itime2 = 0, tdim2
                          itime3 = 0

                          ! Reinitialize from the distribution at the end of t2 propagation
                          x_bath(:) = x_bath_save( itime1, itime2, : )
                          p_bath(:) = p_bath_save( itime1, itime2, : )
                          initstate = istate_save( itime1, itime2 )
                          initstatet = istatet_save( itime1, itime2 )
                          CALL focus_pldm_map( x_map, p_map, xt_map, pt_map )

                          ! GET FINAL RESPONSE AT T3 = 0
                          resp_func( :, itime1, itime2, itime3 ) = resp_func( :, itime1, itime2, itime3 ) &
                                  + nonlinear_response( x_map, p_map, xt_map, pt_map, monte_carlo_weight( itime1, itime2 ), &
                                          k1_sign, k2_sign( itime1 ), k3_sign( itime1, itime2 ) )

                          ! Find the initial bath force
                          bath_force = bilinear_harmonic_force_pldm(x_bath, x_map, p_map, xt_map, pt_map)

                          ! Propagate
                          DO istep=1, nbstep3

                              ! First half of the verlet
                              x_bath = x_bath + p_bath * dt_bath + bath_force * 0.5_dp * dt_bath ** 2
                              p_bath = p_bath + bath_force * 0.5_dp * dt_bath

                              ! Update the full hamiltonian
                              ham = diabatic_bilinear_coupling_hamiltonian(x_bath, c)

                              ! Propagate the mapping variable
                              CALL verlet_mapping_variables(x_map, p_map, ham, dt_map)
                              CALL verlet_mapping_variables(xt_map, pt_map, ham, dt_map)

                              ! If the step is divisible by dump, compute the redmat weighted by initial dipolt
                              IF ( MOD( istep , branch3 ) == 0 ) THEN
                                  itime3 = itime3 + 1

                                  ! GET FINAL RESPONSE AT T3 > 0
                                  resp_func( :, itime1, itime2, itime3 ) = resp_func( :, itime1, itime2, itime3 ) &
                                          + nonlinear_response( x_map, p_map, xt_map, pt_map, monte_carlo_weight( itime1, itime2 ),&
                                                   k1_sign, k2_sign( itime1 ), k3_sign( itime1, itime2 ) )

                              END IF

                              ! Update the force and finish the verlet
                              bath_force = bilinear_harmonic_force_pldm(x_bath, x_map, p_map, xt_map, pt_map)
                              p_bath = p_bath + bath_force * 0.5_dp * dt_bath

                          END DO

                      END DO

                  END DO

              END DO

          END IF

      END DO

  END DO


  ! Write the response function give time increments in fs
  dt_bath_1 = dt_bath_1 * convert('au_time','fs')
  dt_bath_2 = dt_bath_2 * convert('au_time','fs')
  dt_bath_3 = dt_bath_3 * convert('au_time','fs')
  CALL write_nonlinear_response( method, resp_func, dt_bath_1 * branch1, &
          dt_bath_2 * branch2, dt_bath_3 * branch3 )


  ! Deallocate spectroscopy arrays
  CALL finalize_spectroscopy
  DEALLOCATE( monte_carlo_weight, resp_func, x_bath_save, p_bath_save, &
          istate_save, istatet_save, k2_sign, k3_sign )

END SUBROUTINE calculate_pldm_nonlinear
