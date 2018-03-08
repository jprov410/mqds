PROGRAM test_mapping_variables
    USE mapping_variables
    USE kinds
    USE input_output
    IMPLICIT NONE
    REAL(dp), ALLOCATABLE :: test(:,:)
    REAL(dp) :: prop_check
    nstate = 2
    nlit=10


    ALLOCATE( test(nstate,nstate) )
    test = 0.0_dp

    !! TEST PLDM-SPECIFIC SUBROUTINES/FUNCTIONS AND PROPAGATION

    CALL initialize_pldm_map

    CALL sample_pldm_map( x_map, p_map, xt_map, pt_map )
    test = pldm_redmat( x_map, p_map, xt_map, pt_map )
    IF( test(1,1) == 0.0_dp ) THEN
        WRITE(*,*) 'pldm_redmat function not working correctly!'
        STOP 1
    END IF

    test = 1.0_dp ; test(nstate,nstate) = 0.0_dp
    prop_check = x_map(1)
    CALL verlet_mapping_variables( x_map, p_map, test, 10.0_dp )
    IF ( x_map(1) == prop_check ) THEN
        WRITE(*,*) 'verlet mapping propagation not working correctly!'
        STOP 1
    END IF

    CALL finalize_pldm_map

    !! TEST TWA-SPECIFIC SUBROUTINES/FUNCTIONS

    CALL initialize_twa_map

    CALL sample_twa_map( x_map, p_map )
    test = twa_redmat( x_map, p_map )
    IF( test(1,1) == 0.0_dp ) THEN
        WRITE(*,*) 'twa_redmat function not working correctly!'
        STOP 1
    END IF

    CALL finalize_twa_map

    DEALLOCATE( test )
END PROGRAM test_mapping_variables