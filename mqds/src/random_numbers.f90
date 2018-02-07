!
! This module contains the necessary functions
! for random number generation
!
MODULE random_numbers
  IMPLICIT NONE

  PRIVATE
  PUBLIC uniform_rn, gaussian_rn, initialize_rn

  CONTAINS
    SUBROUTINE initialize_rn(my_pe)
      ! Initialize the random number generator
      IMPLICIT NONE
      INTEGER :: time(12)
      INTEGER :: seed(12)
      INTEGER, OPTIONAL, INTENT(in) :: my_pe
      time = 0 ; seed = 0 

      CALL DATE_AND_TIME(VALUES=time)
      seed(:) = time(:)
      IF ( PRESENT(my_pe) ) seed(:) = my_pe + time(:)
      CALL RANDOM_SEED(PUT=seed)

    END SUBROUTINE initialize_rn
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    FUNCTION uniform_rn(input) RESULT(res)
      ! Uniform (0,1] random number generation
      USE parameters
      IMPLICIT NONE
      REAL(dp) :: input(:)
      REAL(dp) :: res(SIZE(input))
      
      CALL RANDOM_NUMBER(HARVEST=res)
      
    END FUNCTION uniform_rn
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    FUNCTION gaussian_rn(input) RESULT(res)
      ! Gaussian random number generation using
      ! the Box-Mueller process
      USE parameters
      IMPLICIT NONE
      REAL(dp) :: input(:)
      REAL(dp) :: res(SIZE(input))
      real(dp) :: make_grn(2,SIZE(input))

      CALL RANDOM_NUMBER(HARVEST=make_grn)
      res(:) = DSQRT( -2.0_dp * DLOG(make_grn(1,:)) ) * DCOS( 2.0_dp * pi * make_grn(2,:))

    END FUNCTION gaussian_rn

END MODULE random_numbers
