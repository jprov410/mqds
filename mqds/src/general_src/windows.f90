!> Module that contains the subroutines and functions
!! to perform SQC windowing calculations.
MODULE windows
    USE kinds
    USE input_output
    USE random_numbers
    USE mapping_variables
    IMPLICIT NONE
    REAL(dp), ALLOCATABLE, PRIVATE :: toll(:), od_toll(:)


CONTAINS

    !> Allocate the windowing arrays and set initial tolls
    !! (for normalizing the density matrix) to 0.
    SUBROUTINE initialize_windows
        IMPLICIT NONE

        ALLOCATE( toll( 0 : nbstep / dump ), od_toll( 0 : nbstep / dump ) )
        toll = 0.0_dp ; od_toll = 0.0_dp

    END SUBROUTINE initialize_windows

    !> Deallocate the windowing arrays.
    SUBROUTINE finalize_windows
        IMPLICIT NONE

        DEALLOCATE( toll, od_toll )

    END SUBROUTINE finalize_windows

    !> Sample the initial distribution of mapping variables for
    !! truncated wigner approximation in action-angle variables
    !! with broadened windows.
    SUBROUTINE sample_sqc_map(x_init, p_init)
    USE random_numbers
    USE parameters
    USE input_output
    USE kinds
    IMPLICIT NONE
    INTEGER :: istate
    REAL(dp), INTENT(out) :: x_init(nstate), p_init(nstate)
    REAL(dp) :: n(nstate), q(nstate), sum_action
    n = 0.0_dp ; q = 0.0_dp ; sum_action = 0.0_dp
    prod = (1.0_dp, 0.0_dp)

    IF ( windowshape /= 'square' .AND. windowshape /= 'triangle' ) THEN
        OPEN(UNIT=10, FILE=ERRORLOG)
        WRITE(10,*) 'You have selected a nonexistent window shape, the &
                input should read e.g.) "windowshape square"'
        CLOSE(10)
        STOP
    END IF


    IF ( windowshape == 'square' ) THEN
        ! Sample rectangular windows
        n = ( 1.0_dp - 2.0_dp * uniform_rn(n) ) * window
        n(initstate) = n(initstate) + 0.5_dp
        n(initstatet) = n(initstatet) + 0.5_dp
    ELSE IF ( windowshape == 'triangle' ) THEN
        window = 1.0_dp / 3.0_dp ; zpe = 1.0_dp / 3.0_dp
        ! Sample triangular windows
        100 continue
        n = uniform_rn( n )
        sum_action = SUM( n )
        IF ( sum_action > 1.0_dp ) THEN
            GO TO 100
        END IF
        n = n - window
        n(initstate) = n(initstate) + 0.5_dp
        n(initstatet) = n(initstatet) + 0.5_dp
    END IF


    ! Sample uniform angle from 0 to 2pi
    q(:) = 2.0_dp * pi * uniform_rn(q)

    IF ( initstate /= initstatet ) THEN
        prod =  EXP( eye * ( q(initstate) - q(initstatet) ) )
    END IF
    ! treat the zero-point energy and window as two different parameters
    ! that are set in the input file
    p_init = -DSQRT( 2.0_dp * ( n + zpe ) ) * DSIN( q )
    x_init = DSQRT( 2.0_dp * ( n + zpe ) ) * DCOS( q )

    END SUBROUTINE sample_sqc_map

    !> Calculate the reduced density matrix using SQC
    FUNCTION sqc_redmat(x, p, itime) RESULT( res )
        USE kinds
        USE parameters
        USE input_output
        IMPLICIT NONE
        INTEGER :: i, j, k, itime, icount
        REAL(dp), INTENT(in) :: x(nstate), p(nstate)
        REAL(dp) :: n(nstate), q(nstate)
        REAL(dp) :: sum_action
        COMPLEX(dp) :: res( nstate, nstate )
        res = (0.0_dp, 0.0_dp) ; sum_action = 0.0_dp
        n = 0.5_dp * ( p**2 + x**2 ) - zpe
        q = -ATAN2( p, x )
        icount = 0

        ! loop over states to check for populations and coherences
        IF ( windowshape == 'square' ) THEN
            DO i=1, nstate
            !~~~~POPULATIONS~~~~!
            ! Check if state i is close enough to action = 1
            IF ( n(i) < 1.0_dp + window .AND. n(i) >= 1.0_dp - window ) THEN
                DO j=1, nstate
                    ! Check if all other states are close enough to action = 0
                    IF ( j /= i ) THEN
                        IF ( n(j) > window .OR. n(j) <= -window ) EXIT
                    END IF
                END DO
                res(i,i) = 1.0_dp * prod
                toll(itime) = toll(itime) + 1.0_dp
            END IF

            !~~~~~COHERENCES~~~~!
            ! check if state i is close enough to action = 0.5
            IF ( n(i) < 0.5_dp + window .AND. n(i) >= 0.5_dp - window ) THEN
                DO j=1, nstate
                    IF ( j /= i ) THEN
                        ! Make sure state j is close enough to action = 0.5
                        IF ( n(j) > 0.5_dp + window .OR. n(j) <= 0.5_dp - window ) CYCLE
                        DO k=1, nstate
                            IF ( k /= i .AND. k /= j ) THEN
                                ! Make sure rest of states are close enough to action = 0
                                IF ( n(k) > window .OR. n(k) <= -window ) GO TO 99
                            END IF
                        END DO
                        res(i,j) = EXP( eye * ( q(j) - q(i) ) ) * prod
                        icount = icount + 1
                        IF ( icount == 1 ) od_toll(itime) = od_toll(itime) + 2.0_dp
                        99 CONTINUE
                    END IF
                END DO
            END IF
        END DO
        END IF

        IF ( windowshape == 'triangle' ) THEN
            sum_action = SUM( n )
            IF( sum_action <= 2.0_dp - nstate * window ) THEN
                DO i = 1, nstate
                    ! Check for populations
                    IF ( n(i) > 1.0_dp - window ) THEN
                        DO j = 1, nstate
                            IF ( j /= i ) THEN
                                IF ( n(j) <= -window ) THEN
                                    EXIT
                                END IF
                            END IF
                        END DO
                        res(i,i) = 1.0_dp * prod
                        toll(itime) = toll(itime) + 1.0_dp
                    END IF

                    ! Check for coherences
                    IF ( n(i) > 0.5_dp - window) THEN
                        DO j = 1, nstate
                            IF ( j /= i ) THEN
                                IF ( n(j) <= 0.5_dp - window ) CYCLE
                                DO k = 1, nstate
                                    IF ( k /= i .AND. k /= j ) THEN
                                        IF ( n(k) <= - window ) GO TO 88
                                    END IF
                                END DO
                                res(i,j) = EXP( eye * ( q(j) - q(i) ) ) * prod
                                icount = icount + 1
                                IF ( icount == 1 ) od_toll(itime) = od_toll(itime) + 2.0_dp
                                88 CONTINUE
                            END IF
                        END DO
                    END IF
                END DO
            END IF
        END IF


    END FUNCTION sqc_redmat

    !> Normalize the reduced density matrix using toll/od_toll from sqc_redmat
    SUBROUTINE normalize_sqc_redmat( redmat )
        IMPLICIT NONE
        COMPLEX(dp) :: redmat( nstate, nstate, 0 : nbstep / dump )
        INTEGER :: i, j, itime

        DO itime=0, nbstep / dump, 1
            IF ( toll(itime) == 0 ) toll(itime) = 1.0_dp
            IF ( od_toll(itime) == 0 ) od_toll(itime) = 1.0_dp
        END DO

        DO i=1, nstate
            DO j=1, nstate
                IF ( i == j ) THEN
                    redmat(j, i, :) = redmat(j, i, :) / toll(:)
                ELSE
                    redmat(j, i, :) = redmat(j, i, :) / od_toll(:)
                END IF
            END DO
        END DO

    END SUBROUTINE normalize_sqc_redmat


END MODULE windows