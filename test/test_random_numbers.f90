PROGRAM test_random_numbers
    USE random_numbers
    USE kinds
    IMPLICIT NONE
    REAL(dp) :: unif1(1), unif2(1)
    REAL(dp) :: gran1(1), gran2(1)

    call initialize_rn
    ! Check out the uniform random numbers
    unif1 = uniform_rn(unif1)
    unif2 = uniform_rn(unif2)
    IF ( unif1(1) == unif2(1) ) THEN
        WRITE(*,*) 'uniform random numbers are identical!'
        STOP 1
    END IF
    ! Check out the gaussian random numbers
    gran1 = gaussian_rn(gran1)
    gran2 = gaussian_rn(gran2)
    IF ( gran1(1) == gran2(1) ) THEN
        WRITE(*,*) 'gaussian random numbers are identical!'
        STOP 1
    END IF

    call initialize_rn
    unif2 = uniform_rn(unif2)
    IF ( unif1(1) == unif2(1) ) THEN
        WRITE(*,*) 'subsequent seeds are identical!'
        STOP 1
    END IF


END PROGRAM test_random_numbers