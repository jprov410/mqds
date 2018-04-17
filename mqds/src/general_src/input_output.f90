!> This is a module that contains input/output
!! filenames, variables, and subroutines for the mqds package
MODULE input_output
  USE kinds
  IMPLICIT NONE
  ! Filenames
  CHARACTER(20), PARAMETER :: INPUT='run.in'
  CHARACTER(20), PARAMETER :: INPHEL='hel.in'
  CHARACTER(20), PARAMETER :: COUPLINGS='couplings.in'
  CHARACTER(20), PARAMETER :: SPECTRALIN='continuumsd.in'
  CHARACTER(20), PARAMETER :: SPECTRALOUT='sampledbath.log'
  CHARACTER(20), PARAMETER :: DIPOLEIN='dipole.in'
  CHARACTER(20), PARAMETER :: ERRORLOG='error.log'
  CHARACTER(20), PARAMETER :: EQUILIBRIUM='equilibrium_site.out'

  ! General variables from input file
  INTEGER :: ntraj, nbstep, nlit, dump, nstate, initstate, initstatet
  INTEGER :: ndof, tdim, ndisc, nslice
  INTEGER :: nbath, nosc
  REAL(dp) :: runtime, temperature
  REAL(dp) :: zpe, window
  
  ! Spectroscopy variables from input file
  INTEGER :: nbstep1, nbstep2, nbstep3
  INTEGER :: branch1, branch2, branch3
  REAL(dp) :: tdelay1, tdelay2, tdelay3
  
  ! Calculation type i.e.) redmat, nonlinear spectra ...
  CHARACTER(20) :: calculation

  ! Dynamics method
  CHARACTER(20) :: method

  ! SQC window shape
  CHARACTER(20) :: windowshape

  ! Basis set for system dynamics
  CHARACTER(20) :: basis
  
  ! Bath type
  CHARACTER(20) :: bath
  
  ! Dummy character
  CHARACTER(1) :: cdum
  
CONTAINS
  !SUBROUTINE read_input
  !  !> Directly reads input from variables in alphabetical order
  !  !! from the "processed_run.in" file
  !  USE kinds
  !  IMPLICIT NONE

  !  OPEN(UNIT=10,FILE=INPUT)
  !  READ(10,*) cdum, basis
  !  READ(10,*) cdum, bath
  !  READ(10,*) cdum, branch1
  !  READ(10,*) cdum, branch2
  !  READ(10,*) cdum, branch3
  !  READ(10,*) cdum, calculation
  !  READ(10,*) cdum, dump
  !  READ(10,*) cdum, initstate
  !  READ(10,*) cdum, initstatet
  !  READ(10,*) cdum, method
  !  READ(10,*) cdum, nbstep
  !  READ(10,*) cdum, nlit
  !  READ(10,*) cdum, nslice
  !  READ(10,*) cdum, nstate
  !  READ(10,*) cdum, nbstep1
  !  READ(10,*) cdum, nbstep2
  !  READ(10,*) cdum, nbstep3
  !  READ(10,*) cdum, ntraj
  !  READ(10,*) cdum, runtime
  !  READ(10,*) cdum, tdelay1
  !  READ(10,*) cdum, tdelay2
  !  READ(10,*) cdum, tdelay3
  !  READ(10,*) cdum, temperature
  !  READ(10,*) cdum, window
  !  READ(10,*) cdum, windowshape
  !  READ(10,*) cdum, zpe
  !  CLOSE(10)

  !END SUBROUTINE read_input

  !> Writes the reduced density matrix using the method,
  !! redmat, and size of timestep between prints (fs)
  SUBROUTINE write_redmat(method, redmat, printstep)
    USE kinds
    IMPLICIT NONE
    INTEGER :: i, j, itime
    INTEGER :: nstate, tdim, namesize
    COMPLEX(dp), INTENT(in) :: redmat(:,:,:)
    REAL(dp), INTENT(in) :: printstep
    CHARACTER(20) :: method
    CHARACTER(17) :: writeformat
    CHARACTER(50) :: filename
    tdim = SIZE(redmat(1,1,:))
    nstate = SIZE(redmat(:,1,1))

    ! Figure out how long the name is going to be
    namesize=LEN_TRIM(method)+2*(INT(LOG10(REAL(nstate)))+1)+1

    ! filename format based on the number of states
    if(nstate >= 1 .AND. nstate < 10) writeformat='(A,A,I1.1,A,I1.1)'
    if(nstate >= 10 .AND. nstate < 100) writeformat='(A,A,I2.2,A,I2.2)'
    if(nstate >= 100 .AND. nstate < 1000) writeformat='(A,A,I3.3,A,I3.3)'
    if(nstate >= 1000 .AND. nstate < 10000) writeformat='(A,A,I4.4,A,I4.4)'
    if(nstate >= 10000 .AND. nstate < 100000) writeformat='(A,A,I5.5,A,I5.5)'

    ! create file for each redmat element and write
    DO i=1, nstate
      DO j=1, nstate
        WRITE(filename, writeformat) TRIM(method),'.',i,'-',j
        OPEN(UNIT=10, FILE=filename)
        DO itime=1, tdim
          WRITE(10,'(3(E13.6,2X))') printstep * (itime-1), redmat(i, j, itime)
        END DO
        CLOSE(10)
      END DO
    END DO


  END SUBROUTINE write_redmat

  !> Writes the linear response function using the method,
  !! response function, and size of timestep between prints (fs)
  SUBROUTINE write_linear_response(method, resp_func, printstep)
    USE kinds
    IMPLICIT NONE
    INTEGER :: i, j, itime
    INTEGER :: tdim, namesize
    COMPLEX(dp), INTENT(in) :: resp_func( 0 : nbstep / dump )
    REAL(dp), INTENT(in) :: printstep
    CHARACTER(20) :: method
    CHARACTER(17) :: writeformat
    CHARACTER(50) :: filename

    ! Figure out how long the name is going to be
    namesize=LEN_TRIM(method)+16
    writeformat='(A,A)'

    ! create file for each redmat element and write
    WRITE(filename, writeformat) TRIM(method),'_linrespfunc.out'
    OPEN(UNIT=10, FILE=filename)
    DO itime=0, ( nbstep / dump )
      WRITE(10,'(3(E13.6,2X))') printstep * (itime), REAL(resp_func(itime)), AIMAG(resp_func(itime))
    END DO
    CLOSE(10)

  END SUBROUTINE write_linear_response

  SUBROUTINE write_equilibrium_site_populations( populations )
    USE kinds
    IMPLICIT NONE
    INTEGER :: istate
    REAL(dp), INTENT(in) :: populations(:)


    OPEN(UNIT=10, FILE=EQUILIBRIUM)
    DO istate = 1, nstate
      WRITE(10, '(E13.6)') populations(istate)
    END DO
    CLOSE(10)

  END SUBROUTINE write_equilibrium_site_populations


  SUBROUTINE write_nonlinear_response( method, response, dt1, dt2, dt3 )
    USE kinds
    IMPLICIT NONE
    INTEGER :: it1, it2, it3
    INTEGER :: namesize
    REAL(dp), INTENT(in) :: dt1, dt2, dt3
    COMPLEX(dp), INTENT(in) :: response( 8, 0 : INT(nbstep1/branch1),&
            0 : INT(nbstep2/branch2), 0 : INT(nbstep3/branch3) )
    CHARACTER(20) :: method
    CHARACTER(17) :: writeformat
    CHARACTER(50) :: filename
    ! Figure out how long the name is going to be
    namesize=LEN_TRIM(method)+14
    writeformat='(A,A)'


    !dimension 1,2 --> K1(-++),K1(+--)
    WRITE(filename, writeformat) TRIM(method),'_nonlin-++.out'
    OPEN(UNIT=10, FILE=filename)
    DO it2 = 0, INT(nbstep2/branch2)
      DO it1 = 0, INT(nbstep1/branch1)
        DO it3 = 0, INT(nbstep3/branch3)
          WRITE(10,'(5(E13.6,2X))') it2 * dt1, it2 * dt2, it3 * dt3, REAL(response(1,it1,it2,it3)), AIMAG(response(1,it1,it2,it3))
        END DO
      END DO
    END DO
    CLOSE(10)

    WRITE(filename, writeformat) TRIM(method),'_nonlin+--.out'
    OPEN(UNIT=10, FILE=filename)
    DO it2 = 0, INT(nbstep2/branch2)
      DO it1 = 0, INT(nbstep1/branch1)
        DO it3 = 0, INT(nbstep3/branch3)
          WRITE(10,'(5(E13.6,2X))') it2 * dt1, it2 * dt2, it3 * dt3, REAL(response(2,it1,it2,it3)), AIMAG(response(2,it1,it2,it3))
        END DO
      END DO
    END DO
    CLOSE(10)

    !dimension 3,4 --> K2(+-+),K2(-+-)
    WRITE(filename, writeformat) TRIM(method),'_nonlin+-+.out'
    OPEN(UNIT=10, FILE=filename)
    DO it2 = 0, INT(nbstep2/branch2)
      DO it1 = 0, INT(nbstep1/branch1)
        DO it3 = 0, INT(nbstep3/branch3)
          WRITE(10,'(5(E13.6,2X))') it2 * dt1, it2 * dt2, it3 * dt3, REAL(response(3,it1,it2,it3)), AIMAG(response(3,it1,it2,it3))
        END DO
      END DO
    END DO
    CLOSE(10)

    WRITE(filename, writeformat) TRIM(method),'_nonlin-+-.out'
    OPEN(UNIT=10, FILE=filename)
    DO it2 = 0, INT(nbstep2/branch2)
      DO it1 = 0, INT(nbstep1/branch1)
        DO it3 = 0, INT(nbstep3/branch3)
          WRITE(10,'(5(E13.6,2X))') it2 * dt1, it2 * dt2, it3 * dt3, REAL(response(4,it1,it2,it3)), AIMAG(response(4,it1,it2,it3))
        END DO
      END DO
    END DO
    CLOSE(10)

    !dimension 5,6 --> K3(++-),K3(--+)
    WRITE(filename, writeformat) TRIM(method),'_nonlin++-.out'
    OPEN(UNIT=10, FILE=filename)
    DO it2 = 0, INT(nbstep2/branch2)
      DO it1 = 0, INT(nbstep1/branch1)
        DO it3 = 0, INT(nbstep3/branch3)
          WRITE(10,'(5(E13.6,2X))') it2 * dt1, it2 * dt2, it3 * dt3, REAL(response(5,it1,it2,it3)), AIMAG(response(5,it1,it2,it3))
        END DO
      END DO
    END DO
    CLOSE(10)

    WRITE(filename, writeformat) TRIM(method),'_nonlin--+.out'
    OPEN(UNIT=10, FILE=filename)
    DO it2 = 0, INT(nbstep2/branch2)
      DO it1 = 0, INT(nbstep1/branch1)
        DO it3 = 0, INT(nbstep3/branch3)
          WRITE(10,'(5(E13.6,2X))') it2 * dt1, it2 * dt2, it3 * dt3, REAL(response(6,it1,it2,it3)), AIMAG(response(6,it1,it2,it3))
        END DO
      END DO
    END DO
    CLOSE(10)

    !dimension 7,8 --> K4(+++),K4(---)
    WRITE(filename, writeformat) TRIM(method),'_nonlin+++.out'
    OPEN(UNIT=10, FILE=filename)
    DO it2 = 0, INT(nbstep2/branch2)
      DO it1 = 0, INT(nbstep1/branch1)
        DO it3 = 0, INT(nbstep3/branch3)
          WRITE(10,'(5(E13.6,2X))') it2 * dt1, it2 * dt2, it3 * dt3, REAL(response(7,it1,it2,it3)), AIMAG(response(7,it1,it2,it3))
        END DO
      END DO
    END DO
    CLOSE(10)

    WRITE(filename, writeformat) TRIM(method),'_nonlin---.out'
    OPEN(UNIT=10, FILE=filename)
    DO it2 = 0, INT(nbstep2/branch2)
      DO it1 = 0, INT(nbstep1/branch1)
        DO it3 = 0, INT(nbstep3/branch3)
          WRITE(10,'(5(E13.6,2X))') it2 * dt1, it2 * dt2, it3 * dt3, REAL(response(8,it1,it2,it3)), AIMAG(response(8,it1,it2,it3))
        END DO
      END DO
    END DO
    CLOSE(10)

  END SUBROUTINE write_nonlinear_response

  SUBROUTINE read_input
    USE kinds
    IMPLICIT NONE
    INTEGER, PARAMETER :: numkeywords=26
    CHARACTER( LEN=11 ) :: keyword( numkeywords )
    CHARACTER( LEN=1 ) :: cdum
    CHARACTER( LEN=80 ) :: line
    INTEGER :: stat( numkeywords ), location( numkeywords )
    INTEGER :: iword, iline, i, j
    CHARACTER( LEN=: ), ALLOCATABLE :: trimword
    !=============================INITIALIZE=============================!
    method = 'pldm' ; calculation = 'redmat' ; basis = 'diabatic'
    bath = 'harmonic' ; branch1 = 1 ; branch2 = 1 ; branch3 = 1
    dump = 1 ; initstate = 1 ; initstatet = 1 ; nbstep = 1000
    nlit = 20 ; nslice = 1 ; nstate = 2
    nbstep1 = 500 ; nbstep2 = 0 ; nbstep3 = 500
    ntraj = 1000 ; runtime = 1000.0_dp
    tdelay1 = 500.0_dp ; tdelay2 = 0.0_dp ; tdelay3 = 500.0_dp
    temperature = 77.0_dp ; window = 0.366_dp ; windowshape = 'square'
    zpe = 0.366_dp
    !=============================INITIALIZE=============================!

    !===GO TO KEYWORD LIBRARY===!
    GO TO 123
    321 CONTINUE
    !===GO TO KEYWORD LIBRARY===!

    OPEN( UNIT=10, file=INPUT )
    DO iword = 1, numkeywords

      REWIND(10) ; iline=0

      stat( iword ) = 0

      IF ( allocated( trimword ) .EQV. .TRUE. ) THEN
        DEALLOCATE( trimword )
      END IF

      ALLOCATE( CHARACTER( LEN=len_trim( keyword(iword) ) ) :: trimword )

      trimword = trim( keyword(iword) )

      20 CONTINUE !CHECK EACH LINE FOR THIS KEYWORD

      iline=iline+1

      READ( 10, '(a)', err=30, end=30 ) line

      stat( iword ) = index( line, trimword )

      IF ( stat( iword ) /= 0 ) THEN
        location( iword ) = iline
        GO TO 30
      END IF
      GO TO 20
      30 CONTINUE
    END DO

    !===ASSIGN INPUT VALUES===!
    GO TO 234
    432 CONTINUE
    !===ASSIGN INPUT VALUES===!

    CLOSE(10)

    !=============END OF SUBROUTINE=============!
    RETURN
    !=============================================!
    !===============KEYWORD LIBRARY===============!
    !=============================================!
    123 CONTINUE
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !  IF THIS IS CHANGED INPUT VALUE ASSIGNMENT  !
    !  SECTION MUST ALSO BE CHANGED ACCORDINGLY   !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !CHARACTERS
    keyword(1) ='basis      '
    keyword(2) ='bath       '
    keyword(3) ='branch1    '
    keyword(4) ='branch2    '
    keyword(5) ='branch3    '
    keyword(6) ='calculation'
    keyword(7) ='dump       '
    keyword(8) ='initstate  '
    keyword(9) ='initstatet '
    keyword(10)='method     '
    keyword(11)='nbstep     '
    keyword(12)='nlit       '
    keyword(13)='nslice     '
    keyword(14)='nstate     '
    keyword(15)='nbstep1    '
    keyword(16)='nbstep2    '
    keyword(17)='nbstep3    '
    keyword(18)='ntraj      '
    keyword(19)='runtime    '
    keyword(20)='tdelay1    '
    keyword(21)='tdelay2    '
    keyword(22)='tdelay3    '
    keyword(23)='temperature'
    keyword(24)='window     '
    keyword(25)='windowshape'
    keyword(26)='zpe        '
    GO TO 321                                     !
    !=============================================!

    !====================================================!
    !===============INPUT VALUE ASSIGNMENT===============!
    !====================================================!
    234 CONTINUE
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !    IF THIS IS CHANGED THEN THE KEYWORD LIBRARY     !
    !     SECTION MUST ALSO BE CHANGED ACCORDINGLY       !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !==============Basis Type===============!
    REWIND(10)
    IF ( location(1) /= 0 ) THEN
      DO i=1,location(1)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, basis
    END IF
    !=======================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !=============Bath Type============!
    REWIND(10)
    IF (location(2)/=0) THEN
      DO i=1,location(2)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, bath
    END IF
    !===================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !=================Branch 1================!
    REWIND(10)
    IF ( location(3) /= 0 ) THEN
      DO i=1,location(3)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, branch1
    END IF
    !========================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !=================Branch 2===============!
    REWIND(10)
    IF ( location(4) /= 0 ) THEN
      DO i=1,location(4)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, branch2
    end if
    !========================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !=================Branch 3===============!
    REWIND(10)
    IF ( location(5) /=0 ) THEN
      DO i=1,location(5)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, branch3
    END IF
    !========================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !=============Calculation Type================!
    REWIND(10)
    IF ( location(6) /= 0 ) THEN
      DO i=1,location(6)-1
        READ(10,*) cdum
      END DO
      read(10,*) cdum, calculation
    END IF
    !============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !==========number bath steps per observation==============!
    REWIND(10)
    IF ( location(7) /= 0 ) THEN
      DO i=1,location(7)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, dump
    END IF
    !============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !=================initial fwd state====================!
    REWIND(10)
    IF ( location(8) /= 0 ) THEN
      DO i=1,location(8)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, initstate
    END IF
    !============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !=================initial BKWD state====================!
    REWIND(10)
    IF ( location(9) /= 0 ) THEN
      DO i=1,location(9)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, initstatet
    END IF
    !============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !================Method for Calculation====================!
    REWIND(10)
    IF ( location(10) /= 0 ) THEN
      DO i=1,location(10)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, method
    END IF
    !============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !===============Number of bath steps====================!
    REWIND(10)
    IF ( location(11) /= 0) THEN
      DO i=1,location(11)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, nbstep
    END IF
    !===============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !=======Number of system per bath steps=========!
    REWIND(10)
    IF ( location(12) /= 0 ) THEN
      DO i=1,location(12)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, nlit
    END IF
    !===============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !=========Number of slices during time==========!
    REWIND(10)
    IF ( location(13) /= 0 ) THEN
      DO i=1,location(13)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, nslice
    END IF
    !===============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !===========Number of system states=============!
    REWIND(10)
    IF ( location(14) /= 0 ) THEN
      DO i=1,location(14)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, nstate
    END IF
    !===============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !==========Number of bath steps in t1===========!
    REWIND(10)
    IF ( location(15) /= 0 ) THEN
      DO i=1,location(15)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, nbstep1
    END IF
    !===============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !==========Number of bath steps in t2===========!
    REWIND(10)
    IF ( location(16) /= 0 ) THEN
      DO i=1,location(16)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, nbstep2
    END IF
    !===============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !==========Number of bath steps in t3===========!
    REWIND(10)
    IF ( location(17) /= 0 ) THEN
      DO i=1,location(17)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, nbstep3
    END IF
    !===============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !============Number of trajectories=============!
    REWIND(10)
    IF ( location(18) /= 0 ) THEN
      DO i=1,location(18)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, ntraj
    END IF
    !===============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !===========runtime in femtoseconds=============!
    REWIND(10)
    IF ( location(19) /= 0 ) THEN
      DO i=1,location(19)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, runtime
    END IF
    !===============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !==============t1 interval time=================!
    REWIND(10)
    IF ( location(20) /= 0 ) THEN
      DO i=1,location(20)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, tdelay1
    END IF
    !===============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !==============t2 interval time=================!
    REWIND(10)
    IF ( location(21) /= 0 ) THEN
      DO i=1,location(21)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, tdelay2
    END IF
    !===============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !==============t3 interval time=================!
    REWIND(10)
    IF ( location(22) /= 0 ) THEN
      DO i=1,location(22)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, tdelay3
    END IF
    !===============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !===============bath temperature=================!
    REWIND(10)
    IF ( location(23) /= 0 ) THEN
      DO i=1,location(23)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, temperature
    END IF
    !===============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !=========half of square window width===========!
    REWIND(10)
    IF ( location(24) /= 0 ) THEN
      DO i=1,location(24)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, window
    END IF
    !===============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !=========half of square window width===========!
    REWIND(10)
    IF ( location(25) /= 0 ) THEN
      DO i=1,location(25)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, windowshape
    END IF
    !===============================================!

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    !=========zero point energy to include==========!
    REWIND(10)
    IF ( location(26) /= 0 ) THEN
      DO i=1,location(26)-1
        READ(10,*) cdum
      END DO
      READ(10,*) cdum, zpe
    END IF
    !===============================================!

    GO TO 432

  END SUBROUTINE read_input

END MODULE input_output
