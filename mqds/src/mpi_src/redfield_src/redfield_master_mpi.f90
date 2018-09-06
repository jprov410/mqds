!> This is the Redfield Theory master subroutine
!! that determines which Redfield-based calculation to
!! carry out.
SUBROUTINE redfield_master_mpi
    USE input_output
    USE hamiltonians
    USE harmonic_bath

    IF (bath == 'harmonic') THEN
        CALL initialize_bath
        IF (basis == 'diabatic') THEN
            CALL initialize_hel
            IF (calculation == 'redmat') THEN
                PRINT*, 'no subroutine for redmat with redfield yet :('
            END IF
            IF (calculation == 'absorption') THEN
                PRINT*, 'no subroutine for absorption with redfield yet :('
            END IF
            IF (calculation == 'rates') THEN
                CALL calculate_redfield_rates_mpi
            END IF
        END IF
    END IF

END SUBROUTINE redfield_master_mpi