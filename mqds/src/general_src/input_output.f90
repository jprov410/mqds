MODULE input_output
  USE kinds
  IMPLICIT NONE
!
! This is a module that contains input/output
! filenames, variables, and subroutines for the mqds package
!
  ! Filenames
  CHARACTER(20), PARAMETER :: INPUT='processed_run.in'
  CHARACTER(20), PARAMETER :: INPHEL='hel.in'
  CHARACTER(20), PARAMETER :: COUPLINGS='couplings.in'
  CHARACTER(20), PARAMETER :: SPECTRALIN='continuumsd.in'
  CHARACTER(20), PARAMETER :: SPECTRALOUT='sampledbath.log'
  CHARACTER(20), PARAMETER :: DIPOLEIN='dipole.in'
  CHARACTER(20), PARAMETER :: ERRORLOG='error.log'

  ! General variables from input file
  INTEGER :: ntraj, nbstep, nlit, dump, nstate, initstate, initstatet
  INTEGER :: ndof, tdim, ndisc, nslice
  INTEGER :: nbath, nosc
  REAL(dp) :: runtime, temperature
  REAL(dp) :: zpe, window
  
  ! Spectroscopy variables from input file
  INTEGER :: nstep1, nstep2, nstep3, tdim1, tdim2, tdim3
  INTEGER :: branch1, branch2, branch3
  REAL(dp) :: tdelay1, tdelay2, tdelay3
  
  ! Calculation type i.e.) redmat, nonlinear spectra ...
  CHARACTER(20) :: calculation
  
  ! Dynamics method
  CHARACTER(20) :: method
  
  ! Basis set for system dynamics
  CHARACTER(20) :: basis
  
  ! Bath type
  CHARACTER(20) :: bath
  
  ! Dummy character
  CHARACTER(1) :: cdum
  
CONTAINS
   
  SUBROUTINE read_input
    !
    ! Directly reads input from variables in alphabetical order
    !
    USE kinds
    IMPLICIT NONE
    
    OPEN(UNIT=10,FILE=INPUT)
    READ(10,*) cdum, basis
    READ(10,*) cdum, bath
    READ(10,*) cdum, branch1
    READ(10,*) cdum, branch2
    READ(10,*) cdum, branch3
    READ(10,*) cdum, calculation
    READ(10,*) cdum, dump
    READ(10,*) cdum, initstate
    READ(10,*) cdum, initstatet
    READ(10,*) cdum, method
    READ(10,*) cdum, nbstep
    READ(10,*) cdum, nlit
    READ(10,*) cdum, nslice
    READ(10,*) cdum, nstate
    READ(10,*) cdum, nstep1
    READ(10,*) cdum, nstep2
    READ(10,*) cdum, nstep3
    READ(10,*) cdum, ntraj
    READ(10,*) cdum, runtime
    READ(10,*) cdum, tdelay1
    READ(10,*) cdum, tdelay2
    READ(10,*) cdum, tdelay3
    READ(10,*) cdum, temperature
    READ(10,*) cdum, window
    READ(10,*) cdum, zpe
    CLOSE(10)
    
  END SUBROUTINE read_input
  
  ! Writes the reduced density matrix using the method
  ! redmat and size of timestep between prints (fs)
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
  
END MODULE input_output
