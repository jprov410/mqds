!
! This module contains the necessary subroutines
! and functions to deal with the diabatic basis set
! and will be extended to deal with adiabatic basis
!
MODULE hamiltonians
  USE kinds
  IMPLICIT NONE
  REAL(dp), ALLOCATABLE :: hel(:,:)
  
CONTAINS
  
  ! Allocate the electronic Hamiltonian
  SUBROUTINE initialize_hel
    USE input_output
    IMPLICIT NONE
    
    ALLOCATE( hel(nstate,nstate) )
    
  END SUBROUTINE initialize_hel
  
  ! Read the input electronic Hamiltonian
  SUBROUTINE read_hel
    USE input_output
    IMPLICIT NONE
    INTEGER :: i,j
    
    OPEN(UNIT=10, FILE=INPHEL)
    DO i=1,nstate
       READ(10,*) ( hel(j,i), j=1, nstate)
    END DO
    CLOSE(10)
    
  END SUBROUTINE read_hel
  
  ! Deallocate the electronic Hamiltonian
  SUBROUTINE finalize_hel
    USE input_output
    IMPLICIT NONE
    
    DEALLOCATE( hel )
    
  END SUBROUTINE finalize_hel

  ! Takes the bath position and uses electronic hamiltonian to calculate hel + V(x_bath)
  ! leaves out the 1/2 m w^2 x^2 term because constant on diagonal
  FUNCTION diabatic_bilinear_coupling_hamiltonian(x_bath, coupling_matrix) RESULT(res)
    USE kinds
    USE input_output
    INTEGER :: i,j
    REAL(dp), INTENT(in) :: x_bath(:)
    REAL(dp), INTENT(in) :: coupling_matrix(:,:)
    REAL(dp) :: res(nstate, nstate)

    res = hel
    DO i=1, nstate
       DO j=1, nbath*nosc
          res(i, i) = res(i, i) + coupling_matrix(j,i) * x_bath(j)
       END DO
    END DO
    
  END FUNCTION diabatic_bilinear_coupling_hamiltonian


END MODULE hamiltonians
