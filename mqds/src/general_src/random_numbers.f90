!> This module contains the necessary functions
! for random number generation
MODULE random_numbers
  IMPLICIT NONE

  PRIVATE
  PUBLIC uniform_rn, gaussian_rn, initialize_rn

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

END MODULE random_numbers
