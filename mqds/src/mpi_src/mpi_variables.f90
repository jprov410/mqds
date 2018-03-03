MODULE mpi_variables
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER :: ierr, mype, npes, ipe
    INTEGER :: result_size
    INTEGER :: tag = 123
    INTEGER :: status(MPI_STATUS_SIZE)
END MODULE mpi_variables