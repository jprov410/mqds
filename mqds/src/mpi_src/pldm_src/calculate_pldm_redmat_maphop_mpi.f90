!> This subroutine calculates the reduced density matrix
!! using PLDM with a harmonic bath model in the diabatic basis
!! with bilinear coupling.
SUBROUTINE calculate_pldm_redmat_maphop_mpi
  USE kinds
  USE unit_conversions
  USE hamiltonians
  USE mapping_variables
  USE harmonic_bath
  USE input_output
  USE mpi_variables
  USE random_numbers
  USE parameters
  IMPLICIT NONE
  INTEGER :: i, j, k
  INTEGER :: istep, itraj, itime, islice
  REAL(dp) :: beta
  REAL(dp) :: ham( nstate, nstate )
  REAL(dp) :: bath_force( nosc * nbath )
  REAL(dp) :: dt_bath, dt_map, printstep
  REAL(dp) :: x_bath_save( nosc * nbath, INT( ntraj / npes ) )
  REAL(dp) :: p_bath_save( nosc * nbath, INT( ntraj / npes ) )
  REAL(dp) :: bath_force_save( nosc * nbath, INT( ntraj / npes ) )
  COMPLEX(dp) :: redmat( nstate, nstate, 0 : nbstep / dump )
  COMPLEX(dp) :: redmat_save( nstate, nstate )
  COMPLEX(dp) :: coeff_fwd_save( nstate ), coeff_bkwd_save( nstate )
  COMPLEX(dp) :: coeff_fwd( nstate ), coeff_bkwd( nstate )
  COMPLEX(dp) :: sumredmat( nstate, nstate, 0 : nbstep / dump )
  ham = 0.0_dp ; bath_force = 0.0_dp ; dt_bath = 0.0_dp ; dt_map = 0.0_dp
  redmat = ( 0.0_dp , 0.0_dp ) ; beta = 0.0_dp ; x_bath_save = 0.0_dp
  p_bath_save = 0.0_dp ; coeff_fwd = 0.0_dp ; coeff_bkwd = 0.0_dp
  coeff_fwd_save = 0.0_dp ; coeff_bkwd_save = 0.0_dp
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


  DO islice = 1, nslice
      IF (islice == 1 ) THEN
          coeff_fwd( : ) = (0.0_dp,0.0_dp)
          coeff_fwd( 1 ) = (0.7071067812,0.0_dp)
          coeff_fwd( 2 ) = (0.7071067812,0.0_dp)
          coeff_bkwd( : ) = (0.0_dp,0.0_dp)
          coeff_bkwd( 1 ) = (0.7071067812,0.0_dp)
          coeff_bkwd( 2 ) = (0.7071067812,0.0_dp)
          !coeff_fwd( : ) = (0.0_dp,0.0_dp) ; coeff_fwd( initstate ) = (1.0_dp,0.0_dp)
          !coeff_bkwd( : ) = (0.0_dp,0.0_dp) ; coeff_bkwd( initstatet ) = (1.0_dp,0.0_dp)

      ELSE
          coeff_fwd( : ) = coeff_fwd_save(:)
          coeff_bkwd( : ) = coeff_bkwd_save(:)
      END IF

      CALL build_current_cdfs(coeff_fwd, coeff_bkwd)

      DO itraj=1, INT( ntraj / npes )
          itime = (islice - 1) * ( nbstep / nslice ) / dump
          ! Two different initializations, if on second slice initialize
          ! from old bath values, resample mapping and make new initial produce
          IF ( islice == 1 ) THEN

              CALL pldm_map_hop( coeff_fwd, coeff_bkwd , &
                      x_map, p_map, xt_map, pt_map)

              CALL sample_thermal_wigner(x_bath, p_bath, beta)
              ! Calculate the t=0 redmat
              redmat(:, :, itime) = redmat(:, :, itime) + pldm_redmat(x_map, p_map, xt_map, pt_map)
              ! Find the initial bath force
              bath_force = bilinear_harmonic_force_pldm(x_bath, x_map, p_map, xt_map, pt_map)
          ELSE
              CALL pldm_map_hop( coeff_fwd, coeff_bkwd , &
                      x_map, p_map, xt_map, pt_map)

              x_bath( : ) = x_bath_save( :, itraj )
              p_bath( : ) = p_bath_save( :, itraj )
              bath_force( : ) = bath_force_save( :, itraj )
              !print*, p_bath(:)
              !bath_force = bilinear_harmonic_force_pldm( x_bath, x_map, p_map, xt_map, pt_map )
          END IF

          DO istep = 1, INT( nbstep / nslice )

              ! First half of the verlet
              x_bath = x_bath + p_bath * dt_bath + bath_force * 0.5_dp * dt_bath ** 2
              p_bath = p_bath + bath_force * 0.5_dp * dt_bath

              ! Update the full hamiltonian
              ham = diabatic_bilinear_coupling_hamiltonian(x_bath, c)

              ! Propagate the mapping variables
              CALL verlet_mapping_variables(x_map, p_map, ham, dt_map)
              CALL verlet_mapping_variables(xt_map, pt_map, ham, dt_map)

              ! If the step is divisible by dump, compute the redmat
              IF ( MOD( istep, dump) == 0 ) THEN
                  itime = itime + 1
                  redmat(:, :, itime) = redmat(:, :, itime) + pldm_redmat(x_map, p_map, xt_map, pt_map)
              END IF

              ! Update the force and finish the verlet
              bath_force = bilinear_harmonic_force_pldm(x_bath, x_map, p_map, xt_map, pt_map)
              p_bath = p_bath + bath_force * 0.5_dp * dt_bath

          END DO
          !! need to figure out most efficient way to do this part
          coeff_fwd_save( : ) = coeff_fwd_save( : ) + 1.0_dp / DSQRT( 2.0_dp )  &
                  * ( x_map( : ) + eye * p_map( : ) ) * weight_f
          coeff_bkwd_save( : ) = coeff_bkwd_save( : ) + 1.0_dp / DSQRT( 2.0_dp )  &
                  * ( xt_map( : ) - eye * pt_map( : ) ) * weight_b
          x_bath_save( :, itraj ) = x_bath( : )
          p_bath_save( :, itraj ) = p_bath( : )
          bath_force_save( :, itraj ) = bath_force( : )
      END DO
      redmat( :, :, (islice - 1) * nbstep / nslice + 1 : islice * nbstep / nslice ) =  &
              redmat( :, :, (islice - 1) * nbstep / nslice + 1 : islice * nbstep / nslice ) / INT( ntraj / npes )
      !redmat_save( :, : ) = redmat( :, :, itime )
      coeff_fwd_save = coeff_fwd_save /  INT( ntraj / npes )
      coeff_bkwd_save = coeff_bkwd_save /  INT( ntraj / npes )
  END DO

  redmat( :, :, 0 ) = redmat( :, :, 0 ) / INT( ntraj / npes )

  result_size = 2 * SIZE(redmat)

  CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
  IF ( mype /= 0 ) THEN
      CALL MPI_Send(redmat, result_size, MPI_Complex, &
              0, tag, MPI_COMM_WORLD, ierr)
  ELSE
      sumredmat = redmat
      DO ipe=1, npes-1
          CALL MPI_recv(redmat, result_size, MPI_Complex,&
                  ipe, tag, MPI_COMM_WORLD, status,ierr)
          sumredmat = sumredmat + redmat
      END DO
      sumredmat = sumredmat / DBLE(npes)
      CALL write_redmat(method, sumredmat, printstep)
  END IF

  DEALLOCATE( radial_f, radial_b)

  CALL write_redmat(method, redmat, printstep)


END SUBROUTINE calculate_pldm_redmat_maphop_mpi
