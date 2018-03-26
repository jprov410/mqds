!> Module to handle a harmonic-oscillator model for the
!! bath in quantum dynamics simulations
MODULE harmonic_bath
  USE kinds
  IMPLICIT NONE
  REAL(dp), ALLOCATABLE :: x_bath(:), p_bath(:) ! bath phase space DOFs
  REAL(dp), ALLOCATABLE :: omega(:), c(:,:) ! bath frequencies and coupling matrix elements

CONTAINS

    !> Allocates the necessary arrays to describe the bath and system-bath interactions. Also, the
    !! spectral density is uniformly sampled according to frequency-dependent reorganization energy.
  SUBROUTINE initialize_bath
    USE kinds
    USE input_output
    USE parameters
    USE unit_conversions
    IMPLICIT NONE
    INTEGER :: i,j,k,l,imode
    INTEGER :: npts ! Number of points in spectral density file
    INTEGER :: smode ! Sampled mode index
    CHARACTER :: junk ! Dummy character used to navigate input files
    REAL(dp) :: dw ! Discrete frequency element i.e.) freq(2) - freq(1)
    REAL(dp), ALLOCATABLE :: freq(:) ! Frequencies from continuumsd.in file
    REAL(dp), ALLOCATABLE :: sd(:,:), lambda(:,:) ! Spectral density from continuumsd.in, reorg engy
    REAL(dp), ALLOCATABLE :: norm(:) ! Normalization factor for each bath
    REAL(dp), ALLOCATABLE :: cfactor(:,:) ! Factor to turn on/off site<->bath coupling
    REAL(dp), ALLOCATABLE :: sampledw(:,:), sampledc(:,:), sampledj(:,:) ! Samplings from continuumsd.in
    dw = 0.0_dp
    !> Note: This is restricted to a single range and discretization
    !!       of frequencies (defined by a single dw)

    ! Figure out how many lines in spectral density file
    npts = 0
    OPEN(UNIT=10, FILE=SPECTRALIN)
    DO
        READ(10, *, ERR=11, END=11) junk
        npts=npts+1
    END DO
    11  CONTINUE ; CLOSE(10)

    ! Warn if there was no continuum spectral density provided
    IF (npts==0) THEN
        OPEN(unit=20,file=ERRORLOG)
        WRITE(20,*) 'No continuum spectral density file found!'
        CLOSE(20)
    END IF

    ! Figure out how many baths and how many oscillators per bath
    OPEN(UNIT=30,FILE=COUPLINGS)
    READ(30,*) nbath, junk, nosc
    READ(30,*) junk

    ! Allocate all arrays needed to describe the bath
    ALLOCATE(freq(npts), sd(nbath, npts), lambda(nbath, npts), norm(nbath), &
         cfactor(nbath, nstate), sampledw(nbath, nosc), sampledc(nbath, nosc), &
         sampledj(nbath, nosc))

    freq = 0.0_dp ; sd = 0.0_dp ; lambda = 0.0_dp ; cfactor = 0.0_dp
    sampledw = 0.0_dp ; sampledc = 0.0_dp ; sampledj = 0.0_dp ; norm = 0.0_dp

    ! Read in the prefactor matrix to turn on/off site<->bath coupling
    DO i=1, nstate
       READ(30,*) ( cfactor(j,i), j=1, nbath )
    END DO
    CLOSE(30)

    ! Read in the frequency and spectral density for each bath
    OPEN(UNIT=10,FILE=SPECTRALIN)
    DO i=1, npts
        READ(10, *) freq(i), ( sd(j,i), j=1, nbath )
    END DO
    CLOSE(10)

    ! Convert the frequencies and spectral densities to atomic units
    freq = freq * convert('wvnbr','au_ang_freq')
    sd = sd * convert('wvnbr','au_energy')

    dw = freq(2) - freq(1)

    !> Sample bath according to the method described by
    !! Wang et al JCP 110, 4828 (1999)
    IF ( freq(1) /= 0.0_dp ) lambda(:,1) = dw * sd(:,1) / freq(1)
    DO i=1, nbath
       DO j=2, npts
          lambda(i,j) = lambda(i,j-1) + dw * sd(i,j) / freq(j)
       END DO
       norm(i) = DBLE(nosc) / lambda(i,npts)
       lambda(i,:) = lambda(i,:) * norm(i)
    END DO

   ! Start Sampling
    DO i=1, nbath
       smode = 1
       DO imode=1, nosc-1
          DO WHILE ( lambda(i,smode) <= REAL(imode) )
             smode = smode + 1
          END DO
          sampledw(i,imode) = freq(smode)
          sampledc(i,imode) = DSQRT( 2.0_dp / pi / norm(i) ) * freq(smode)
          sampledj(i,imode) = sd(i,smode)
       END DO
       sampledw(i,nosc) = freq(npts)
       sampledc(i,nosc) = DSQRT( 2.0_dp / pi / norm(i) ) * freq(npts)
       sampledj(i,nosc) = sd(i,npts)
    END DO

    ! Allocate coupling matrix, frequency vector, Bath DOF vectors
    ALLOCATE(x_bath(nbath*nosc), p_bath(nbath*nosc), omega(nbath*nosc), c(nbath*nosc, nstate))
    c = 0.0_dp ; omega = 0.0_dp ; x_bath = 0.0_dp ; p_bath = 0.0_dp

    ! Fill coupling and frequency arrays
    DO i=1, nstate
      DO l=1, nbath
         DO j=1, nosc
            k = ( l - 1 ) * nosc + j
            c(k,i) = cfactor(l,i) * sampledc(l,j)
            omega(k) = sampledw(l,j)
         END DO
      END DO
    END DO


   !> Writes a log of sampled modes into a file called "sampledbath.log"
   OPEN(unit=40,FILE=SPECTRALOUT)
   DO j=1, nosc
      WRITE(40,*) 'Freq', ( sampledw(i,j) * convert('au_ang_freq','wvnbr'), i = 1, nbath ) , &
           'Spec', ( sampledj(i,j) * convert('au_energy','wvnbr'), i=1, nbath )
   END DO
   CLOSE(40)
   
   DEALLOCATE(freq, sd, lambda, norm, cfactor, sampledw, sampledc, sampledj)      

 END SUBROUTINE initialize_bath

    !> Sample the initial thermal distribution for harmonic bath from the
    !! Wigner transform of the Boltzmann operator. This is under the assumption that
    !! the full initil density operator can be written as \f$ \hat{\rho}_s\hat{\rho}_b^{eq} \f$.
 SUBROUTINE sample_thermal_wigner(x_init, p_init, beta)
   USE random_numbers
   USE unit_conversions
   USE input_output
   USE kinds
   IMPLICIT NONE
   REAL(dp), INTENT(in) :: beta ! For Boltzmann operator
   REAL(dp), INTENT(out) :: x_init(nbath*nosc), p_init(nbath*nosc) ! Initial bath DOFs

   ! Calcuate width of distribution for x_init, p_init
   p_init(:) = DSQRT( omega(:) / ( 2.0_dp * DTANH( 0.5_dp * beta * omega(:) ) ) )
   x_init(:) = p_init(:) / omega(:)

   ! Randomly sample the Distribution for x_init, p_init
   p_init(:) = p_init(:) * gaussian_rn(p_init)   
   x_init(:) = x_init(:) * gaussian_rn(x_init)

 END SUBROUTINE sample_thermal_wigner

 ! Deallocate the bath arrays
 SUBROUTINE finalize_bath
   IMPLICIT NONE
   
   DEALLOCATE(x_bath, p_bath, omega, c)

 END SUBROUTINE finalize_bath

 !> Calculate bilinear coupling force for PLDM
 FUNCTION bilinear_harmonic_force_pldm(x, x_map, p_map, xt_map, pt_map) RESULT(res)
   USE kinds
   USE input_output
   IMPLICIT NONE
   INTEGER :: i
   REAL(dp) :: res(nbath*nosc)
   REAL(dp), INTENT(in) :: x(nbath*nosc)
   REAL(dp), INTENT(in) :: x_map(nstate), p_map(nstate), xt_map(nstate), pt_map(nstate)
   res(:) = - omega(:) ** 2 * x(:)
   DO i=1, nstate
      res(:) = res(:) - c(:,i) * 0.25_dp * ( x_map(i)**2 + p_map(i)**2 + xt_map(i)**2 + pt_map(i)**2 )
      ! 0.5 from the mapping variables and 0.5 from average of fwd & bkwd
   END DO
 
 END FUNCTION bilinear_harmonic_force_pldm

    !> Calculate the bilinear coupling force for TWA and SQC calculations
 FUNCTION bilinear_harmonic_force_twa(x, x_map, p_map) RESULT(res)
   USE kinds
   USE input_output
   IMPLICIT NONE
   INTEGER :: i
   REAL(dp) :: res(nbath*nosc)
   REAL(dp), INTENT(in) :: x(nbath*nosc)
   REAL(dp), INTENT(in) :: x_map(nstate), p_map(nstate)
   
   res(:) = - omega(:) ** 2 * x(:)
   DO i=1, nstate
      res(:) = res(:) - c(:,i) * 0.5_dp * ( x_map(i)**2 + p_map(i)**2 - 2.0_dp * zpe )
      ! 0.5 from the mapping variables
   END DO

 END FUNCTION bilinear_harmonic_force_twa


END MODULE harmonic_bath
