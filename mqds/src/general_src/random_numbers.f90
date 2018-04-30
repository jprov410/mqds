!> This module contains the necessary functions
! for random number generation
MODULE random_numbers
  USE kinds
  IMPLICIT NONE
  REAL(dp), PRIVATE :: cdf_g(0:40000), cdf_e(0:40000)

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
    SUBROUTINE build_ge_cdfs
      IMPLICIT NONE
      INTEGER :: i
      REAL(dp) :: dr = 0.0001_dp
      cdf_g = 0.0_dp ; cdf_e = 0.0_dp

      ! Build and normalize the cdfs
      DO i = 1, 40000
        cdf_g(i) = cdf_g(i-1) + ( ( i * dr ) * EXP( - (i * dr)**2 / 2.0_dp ) )
        cdf_e(i) = cdf_e(i-1) + ( ( i * dr )**2 * EXP( - (i * dr)**2 / 2.0_dp ) )
      END DO
      cdf_g = cdf_g / cdf_g(40000)
      cdf_e = cdf_e / cdf_e(40000)
    END SUBROUTINE build_ge_cdfs

    !> This function assigns a random number distributed according to the
    !! radial portion of the ground harmonic oscillator state in polar coordinates
    FUNCTION ground_rn() RESULT( res )
      IMPLICIT NONE
      REAL(dp) :: res
      REAL(dp) :: dr = 0.0001_dp
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
      REAL(dp) :: dr = 0.0001_dp
      INTEGER :: sval

      CALL RANDOM_NUMBER( HARVEST = res )

      sval = MINLOC( ABS( cdf_e - res ), 1 )

      res = sval * dr

    END FUNCTION excited_rn

END MODULE random_numbers
