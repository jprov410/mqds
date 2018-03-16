!> This subroutine calculates the reduced density matrix
!! using truncated wigner approximation with a harmonic bath model
!! in the diabatic basis with bilinear coupling.
SUBROUTINE calculate_sqc_redmat_mpi
  USE kinds
  USE unit_conversions
  USE hamiltonians
  USE mapping_variables
  USE windows
  USE harmonic_bath
  USE input_output
  USE mpi_variables
  IMPLICIT NONE
  INTEGER :: i, j, k
  INTEGER :: istep, itraj, itime
  REAL(dp) :: beta
  REAL(dp) :: ham( nstate, nstate )
  REAL(dp) :: bath_force( nosc * nbath )
  REAL(dp) :: dt_bath, dt_map, printstep
  COMPLEX(dp) :: redmat( nstate, nstate, 0 : nbstep / dump )
  COMPLEX(dp) :: sumredmat( nstate, nstate, 0 : nbstep / dump )
  ham = 0.0_dp ; bath_force = 0.0_dp ; dt_bath = 0.0_dp ; dt_map = 0.0_dp
  sumredmat = ( 0.0_dp , 0.0_dp ) ; redmat = ( 0.0_dp , 0.0_dp ) ; beta = 0.0_dp
  
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

  DO itraj=1, INT( ntraj/npes )
     itime=0
     
     ! Sample the initial conditions for system mapping and bath DOFs
     CALL sample_thermal_wigner(x_bath, p_bath, beta)
     CALL sample_sqc_map(x_map, p_map)

     ! Calculate the t=0 redmat
     redmat(:, :, itime) = redmat(:, :, itime) + sqc_redmat(x_map, p_map, itime)

     ! Find the initial bath force
     bath_force = bilinear_harmonic_force_twa(x_bath, x_map, p_map)

     DO istep=1, nbstep

        ! First half of the verlet
        x_bath = x_bath + p_bath * dt_bath + bath_force * 0.5_dp * dt_bath ** 2
        p_bath = p_bath + bath_force * 0.5_dp * dt_bath

        ! Update the full hamiltonian
        ham = diabatic_bilinear_coupling_hamiltonian(x_bath, c)
        
        ! Propagate the mapping variable
        CALL verlet_mapping_variables(x_map, p_map, ham, dt_map)
   
        ! If the step is divisible by dump, compute the redmat
        IF ( MOD( istep, dump) == 0 ) THEN
           itime = itime + 1           
           redmat(:, :, itime) = redmat(:, :, itime) + sqc_redmat(x_map, p_map, itime)
        END IF

        ! Update the force and finish the verlet
        bath_force = bilinear_harmonic_force_twa(x_bath, x_map, p_map)
        p_bath = p_bath + bath_force * 0.5_dp * dt_bath
        
     END DO
     
  END DO

  CALL normalize_sqc_redmat(redmat)
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



END SUBROUTINE calculate_sqc_redmat_mpi
