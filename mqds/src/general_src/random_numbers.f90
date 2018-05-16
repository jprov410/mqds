!> This module contains the necessary functions
!! for random number generation.
MODULE random_numbers
  USE kinds
  IMPLICIT NONE
  REAL(dp), PRIVATE :: cdf_g(0:10000), cdf_e(0:10000)
  REAL(dp), ALLOCATABLE :: radial_f(:,:), radial_b(:,:)
  !> cdf_g, cdf_e are radial cumulative distribution functions
  !! for the ground and first excited harmonic oscillator states

  !> radial_f, radial_b are the full radial CDFs for forward
  !! and backward mapping variables based on linear combinations of
  !! cdf_g, cdf_e where coefficients come from final mapping variable
  !! values from previous propagation. These are what the mapping variables
  !! are sampled from

  CONTAINS

      !> Inserts the random number generation seed based
      !! upon the current time. If the calculation is being
      !! done with MPI-parallelization, this subroutine takes
      !! the processing element ID number as an argument to ensure
      !! that the seeds are unique.
    SUBROUTINE initialize_rn(my_pe)
      ! Initialize the random number generator
      IMPLICIT NONE
      INTEGER :: time(12), seedsize
      INTEGER, ALLOCATABLE :: seed(:)
      INTEGER, OPTIONAL, INTENT(in) :: my_pe
      time = 0

      CALL RANDOM_SEED(SIZE=seedsize)

      ALLOCATE( seed(seedsize) )
      seed(:) = 0

      CALL DATE_AND_TIME(VALUES=time)

      IF ( seedsize >= size(time) ) THEN
        seed = seed + time(8)
        seed( 1 : size(time) ) = seed( 1 : size(time) ) + time(:)
      ELSE
        seed(:) = seed(:) + time( 1 : seedsize )
      END IF

      IF ( PRESENT(my_pe) ) seed(:) = seed(:) + 5 * my_pe

      CALL RANDOM_SEED(PUT=seed)

      DEALLOCATE( seed )

    END SUBROUTINE initialize_rn

      !> sample a uniform [0,1) random number
    FUNCTION uniform_rn(input) RESULT(res)
      USE parameters
      IMPLICIT NONE
      REAL(dp) :: input(:)
      REAL(dp) :: res(SIZE(input))
      
      CALL RANDOM_NUMBER(HARVEST=res)
      
    END FUNCTION uniform_rn

      !> Gaussian random number generation using
      !! the Box-Mueller process with two uniform
      !! [0,1) random number
    FUNCTION gaussian_rn(input) RESULT(res)
      USE parameters
      IMPLICIT NONE
      REAL(dp) :: input(:)
      REAL(dp) :: res(SIZE(input))
      REAL(dp) :: make_grn(2, SIZE(input))

      CALL RANDOM_NUMBER(HARVEST=make_grn)
      res(:) = DSQRT( -2.0_dp * DLOG(make_grn(1,:)) ) * DCOS( 2.0_dp * pi * make_grn(2,:))
    END FUNCTION gaussian_rn

    !> This is a subroutine that builds the normalized cumulative
    !! distribution functions for the radial portion of
    !! the ground and first excited state harmonic states
    !! in polar coordinates (\f$ (x - i p) \rightarrow r e^{-i \theta} \f$)
    SUBROUTINE build_current_cdfs(c_f,c_b)
    !SUBROUTINE build_current_cdfs(c, islice)
      USE input_output
      IMPLICIT NONE
      INTEGER :: i, j, k, l
      REAL(dp) ::  a_f, a_b
      REAL(dp) :: a
      COMPLEX(dp), INTENT(in) :: c_f( nstate ), c_b( nstate )
      !COMPLEX(dp), INTENT(in) :: c( nstate, nstate )

      ! Build the current radial CDF based on linear combination of
      ! ground and first excited state radial harmonic oscillator CDFs
      DO i = 1, nstate
        ! prepare normalized coefficients for current term
        a_f = SQRT( c_f(i) * CONJG(c_f(i)) )
        a_b = SQRT( c_b(i) * CONJG(c_b(i)) )
        ! add component of excited state for current term
        radial_f(i,:) = radial_f(i,:) + a_f * cdf_e(:)
        radial_b(i,:) = radial_b(i,:) + a_b * cdf_e(:)
        DO j=1,nstate
          IF ( j /= i ) THEN
            ! add component of ground state for current term
            radial_f(j,:) = radial_f(j,:) + a_f * cdf_g(:)
            radial_b(j,:) = radial_b(j,:) + a_b * cdf_g(:)
          END IF
        END DO
      END DO

      !DO i = 1, nstate
        !DO j = 1, nstate
          ! prepare normalized coefficients for current term
      !  j=i

          ! TRY SQRT(A) ALSOOOOO
      !    a = SQRT( c(i,j) * CONJG(c(i,j)) )

      !    a = SQRT( a )

          ! add component of excited state for current term
      !    radial_f(i,:) = radial_f(i,:) + a * cdf_e(:)
      !    radial_b(j,:) = radial_b(j,:) + a * cdf_e(:)
      !    DO k=1,nstate
      !      IF ( k /= i ) THEN
      !        ! add component of ground state for current term
      !        radial_f(k,:) = radial_f(k,:) + a * cdf_g(:)
      !      END IF
      !      IF ( k /= j ) THEN
      !        ! add component of ground state for current term
      !        radial_b(k,:) = radial_b(k,:) + a * cdf_g(:)
      !      END IF
      !    END DO
      !  !END DO
      !END DO

      !normalize
      DO i = 1, nstate
        radial_f(i,:) = radial_f(i,:) / radial_f(i, SIZE(radial_f(i,:)) - 1)
        radial_b(i,:) = radial_b(i,:) / radial_b(i, SIZE(radial_b(i,:)) - 1)
      END DO

    END SUBROUTINE build_current_cdfs

    SUBROUTINE initialize_radial_distribution_functions
      USE input_output
      IMPLICIT NONE
      INTEGER :: i
      REAL(dp) :: dr = 0.0005_dp

      ! allocate memory for radial distributions
      IF ( ALLOCATED( radial_f ) .EQV. .FALSE. ) THEN
        ALLOCATE( radial_f( nstate, 0 : SIZE(cdf_e) -1 ) )
      END IF
      IF ( ALLOCATED( radial_b ) .EQV. .FALSE. ) THEN
        ALLOCATE( radial_b( nstate, 0 : SIZE(cdf_e) -1 ) )
      END IF
      radial_f = 0.0_dp ; radial_b = 0.0_dp

      cdf_g = 0.0_dp ; cdf_e = 0.0_dp
        DO i = 1, SIZE(cdf_e)-1
          cdf_g(i) = cdf_g(i-1) + dr * ( ( i * dr ) * EXP( - (i * dr)**2 / 2.0_dp ) )
          cdf_e(i) = cdf_e(i-1) + dr * ( ( i * dr )**2 * EXP( - (i * dr)**2 / 2.0_dp ) )
        END DO
    END SUBROUTINE initialize_radial_distribution_functions

    !> subroutine to sample the current cumulative distribution
    !! functions for the radial portion of the harmonic oscillator
    !! states
    SUBROUTINE sample_current_cdfs( r, rt )
      USE input_output
      IMPLICIT NONE
      REAL(dp), INTENT(out) :: r(nstate), rt(nstate)
      REAL(dp) :: dr = 0.0005_dp
      INTEGER :: sval,svalt, i ! sampled index values

      CALL RANDOM_NUMBER( HARVEST = r )
      CALL RANDOM_NUMBER( HARVEST = rt )

      DO i = 1, nstate
        sval = MINLOC( ABS( radial_f(i,:) - r(i) ), 1 )
        svalt = MINLOC( ABS( radial_b(i,:) - rt(i) ), 1 )
        r(i) = sval * dr
        rt(i) = svalt * dr
      END DO

    END SUBROUTINE sample_current_cdfs

    !> This function assigns a random number distributed according to the
    !! radial portion of the ground harmonic oscillator state in polar coordinates
    FUNCTION ground_rn() RESULT( res )
      IMPLICIT NONE
      REAL(dp) :: res
      REAL(dp) :: dr = 0.0005_dp
      INTEGER :: sval

      CALL RANDOM_NUMBER( HARVEST = res )

      sval = MINLOC( ABS( cdf_g - res ), 1 )

      res = sval * dr

    END FUNCTION ground_rn

    !> This function assigns a random number distributed according to the
    !! radial portion of the first excited harmonic oscillator state in polar coordinates
    FUNCTION excited_rn() RESULT( res )
      IMPLICIT NONE
      REAL(dp) :: res
      REAL(dp) :: dr = 0.0005_dp
      INTEGER :: sval

      CALL RANDOM_NUMBER( HARVEST = res )

      sval = MINLOC( ABS( cdf_e - res ), 1 )

      res = sval * dr

    END FUNCTION excited_rn

END MODULE random_numbers
