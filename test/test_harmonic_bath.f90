!> Tests the Harmonic Bath module.
PROGRAM test_harmonic_bath
    USE harmonic_bath
    USE kinds
    USE input_output
    USE unit_conversions
    IMPLICIT NONE
    REAL(dp) :: beta
    REAL(dp), ALLOCATABLE :: test(:)
    REAL(dp), ALLOCATABLE :: x(:),p(:),xt(:),pt(:)
    nstate = 2 ; temperature = 77.0_dp
    beta = 1.0_dp / ( temperature * convert('kelvin','au_energy') )

    CALL initialize_bath
    ALLOCATE( test(nbath*nosc), x(nstate), p(nstate), xt(nstate),pt(nstate) )
    test=0.0_dp ; x=0.0_dp ; p=0.0_dp ; xt=0.0_dp ; pt=0.0_dp

    CALL sample_thermal_wigner(x_bath,p_bath,beta)

    IF (x_bath(1) == 0) THEN
        WRITE(*,*) 'Thermal bath not being sampled correctly!'
        STOP 1
    END IF

    test = bilinear_harmonic_force_pldm(x_bath,x,p,xt,pt)
    IF (test(1) == 0) THEN
        WRITE(*,*) 'Bilinear harmonic force pldm returning incorrect value'
        STOP 1
    END IF

    test = bilinear_harmonic_force_twa(x_bath,x,p)
    IF (test(1) == 0) THEN
        WRITE(*,*) 'Bilinear harmonic force twa returning incorrect value'
        STOP 1
    END IF


    CALL finalize_bath
    DEALLOCATE( test, x, p, xt, pt )

END PROGRAM test_harmonic_bath