!> Tests the Hamiltonians module
PROGRAM test_hamiltonians
    USE kinds
    USE input_output
    USE hamiltonians
    USE harmonic_bath
    IMPLICIT NONE
    REAL(dp) :: ham(2,2)
    REAL(dp) :: x(4)
    nstate = 2
    nbath = 2
    nosc = 2
    x = 0.0_dp
    ham = 0.0_dp
    ALLOCATE( c(4,2) )
    c = 0.0_dp

    CALL initialize_hel
    hel = 0.0_dp

    CALL read_hel

    ham = diabatic_bilinear_coupling_hamiltonian(x,c)
    IF ( ham(1,1) /= hel(1,1) ) THEN
        WRITE(*,*) 'diabatic_bilinear_coupling_hamiltonian not working correctly, &
                expected', hel(1,1), 'but got', ham(1,1)
        STOP 1
    END IF

    CALL finalize_hel
    DEALLOCATE( c )
END PROGRAM test_hamiltonians