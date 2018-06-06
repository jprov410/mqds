!> This file contains the subroutines for computing the one-dimensional
!! and 2-dimensional fourier transforms from the linear and non-linear
!! response functions using. It is to be compiled using f2py (included with numpy)
!! and imported within python analysis scripts as a module. The subroutines
!! are then available within the python script
!~~~~~~~~~~  f2py -c -m analysis f2py_analysis.f90  ~~~~~~~~~~~~~~~!

! One-Dimensional Fourier Transform
SUBROUTINE one_d_ft(real_response, imag_response, time, w, output)
  IMPLICIT NONE
  INTEGER :: iw, it
  REAL(8), DIMENSION(:) :: real_response
  REAL(8), DIMENSION(:) :: imag_response
  REAL(8), DIMENSION(:) :: time
  REAL(8), DIMENSION(:) :: w
  COMPLEX(8), DIMENSION(size(w)) :: output
  COMPLEX(8), DIMENSION(size(real_response)) :: response
  REAL(8) :: dw, dt
  REAL(8), PARAMETER :: pi = 3.1415926536
  REAL(8), PARAMETER :: fs2wvnbr = 33356.41
  COMPLEX(8), PARAMETER :: eye = (0.d0, 1.d0)
  output = (0.d0,0.d0) ; response = (0.d0,0.d0)
  !> Declaration of f2py variables
  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
  !f2py intent(in) tpts
  !f2py intent(in) real_response
  !f2py intent(in) imag_response
  !f2py intent(in) time
  !f2py intent(in) w
  !f2py intent(out) output
  
  CALL write_to_screen('NOTE','PERFORMING ONE-DIMENSIONAL FOURIER TRANSFORM')

  ! set up parameters
  w = w * 2.d0 * pi / fs2wvnbr
  dw = w(2) - w(1)
  dt = time(2) - time(1)

  response = DCMPLX( real_response, imag_response )
  
  DO iw = 1, SIZE(w(:))
     DO it = 1, SIZE(response(:))
        output(iw) = output(iw) &
             + dt * exp( -eye * w(iw) * (it-1) * dt ) * response(it)
     END DO
  END DO
  
    w = w / (2.d0 * pi / fs2wvnbr)

END SUBROUTINE one_d_ft


! Two-Dimensional Fourier Transform
!SUBROUTINE two_d_ft(real_response, imag_response, time1, time2, time3, &
!     w1, w3)
!  IMPLICIT NONE
!  INTEGER :: iw1, it1, iw3, it3
!  REAL(8), DIMENSION(:) :: real_response
!  REAL(8), DIMENSION(:) :: imag_response
!  REAL(8), DIMENSION(:) :: time1, time2, time3
!  REAL(8), DIMENSION(:) :: w
!  COMPLEX(8), DIMENSION(size(w)) :: output
!  COMPLEX(8), DIMENSION(size(real_response)) :: response
!  REAL(8) :: dw, dt
!  REAL(8), PARAMETER :: pi = 3.1415926536
!  REAL(8), PARAMETER :: fs2wvnbr = 33356.41
!  COMPLEX(8), PARAMETER :: eye = (0.d0, 1.d0)
!  output = (0.d0,0.d0) ; response = (0.d0,0.d0)
!  !> Declaration of f2py variables
!  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!  !f2py intent(in) tpts
!  !f2py intent(in) real_response
!  !f2py intent(in) imag_response
!  !f2py intent(in) time
!  !f2py intent(in) w
!  !f2py intent(out) output
!
!  
!  CALL write_to_screen('NOTE','PERFORMING TWO-DIMENSIONAL FOURIER TRANSFORM')
!
!END SUBROUTINE two_d_ft

! Apply Cosine Window function to smooth out end of response function
SUBROUTINE apply_cos_window( response, output )
  IMPLICIT NONE
  REAL(8) :: response(:), output( SIZE(response) )
  INTEGER :: itime
  REAL(8) :: pi=3.1415926536
  !> Declaration of f2py variables
  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
  !f2py intent(in) real_response
  !f2py intent(in) imag_response
  !f2py intent(out) output
  
  CALL write_to_screen('NOTE','APPLYING COSINE WINDOW FUNCTION')
  
DO itime = 1, SIZE( response )
     response(itime) = response(itime) &
          * COS( pi * ( REAL(itime-1))/ (2.0 * ( REAL(SIZE(response)) - 1 ) ) )
  END DO
  
  output = response
END SUBROUTINE apply_cos_window

SUBROUTINE write_to_screen( type, message )
  IMPLICIT NONE
  CHARACTER(*) :: type
  CHARACTER (*):: message

  write(*,*) '***',type,'*** : ', message

END SUBROUTINE write_to_screen
