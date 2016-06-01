PROGRAM main
  IMPLICIT NONE
#ifdef MPI
  INCLUDE 'mpif.h'
#endif
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(12,256)

  !-- Parameters of input matrices
  INTEGER, PARAMETER :: numberOfMatrices = 24
  INTEGER, PARAMETER :: sizeOfSquareMatrices = 500 !-- 500 * 500
  
  REAL(r8) :: inputMatrix(numberOfMatrices, sizeOfSquareMatrices, sizeOfSquareMatrices)
  REAL(r8) :: localMatrix(sizeOfSquareMatrices, sizeOfSquareMatrices), outputMatrix(sizeOfSquareMatrices, sizeOfSquareMatrices)

  CHARACTER(LEN=10) :: dummyString !-- some string we will use later

  !-- Define some indices
  INTEGER :: i, j, k, begin_index, end_index, first_matrix_index, second_matrix_index

  !-- MPI parameters
  INTEGER :: nprocs, rank, ierr

#ifdef MPI
  INTEGER status(MPI_STATUS_SIZE)
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
#else
  nprocs = 1
  rank = 0
#endif

  DO i = 1, numberOfMatrices
     WRITE(dummyString, "(i4)") i - 1 !-- Convert int to string

     OPEN(UNIT = 1, FILE = 'matrix_' // TRIM(ADJUSTL(dummyString)) // '.dat', ACTION='READ')
     DO j = 1, sizeOfSquareMatrices
        READ (1, *) inputMatrix(i, j, :) !-- Read a line
     END DO
     CLOSE(UNIT = 1)
  END DO

  !-- Note that there will be numberOfMatrices times numberOfMatrices matrix multiplications to do
  begin_index = (rank * numberOfMatrices * numberOfMatrices + nprocs - 1) / nprocs + 1
  end_index = ( (rank + 1) * numberOfMatrices * numberOfMatrices + nprocs - 1 ) / nprocs
  
  localMatrix = 0.0_r8
  outputMatrix = 0.0_r8
  DO k = begin_index, end_index
     first_matrix_index = (k - 1) / numberOfMatrices + 1
     second_matrix_index = MOD( k - 1, numberOfMatrices ) + 1
     WRITE(*, *) "Core ", rank, " is calculating Matrix ", first_matrix_index - 1, " times ", second_matrix_index - 1
     localMatrix = localMatrix + MATMUL( inputMatrix(first_matrix_index, :, :), inputMatrix(second_matrix_index, :, :) )
  END DO

#ifdef MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(localMatrix, outputMatrix, sizeOfSquareMatrices * sizeOfSquareMatrices, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#else
  outputMatrix = localMatrix
#endif

  IF (rank .EQ. 0 ) THEN
     WRITE(dummyString, "(i4)") sizeOfSquareMatrices !-- Convert int to string

     OPEN(UNIT = 1, FILE = 'matrix_output.dat')
     DO i = 1, sizeOfSquareMatrices
        WRITE(1, "(" // TRIM(ADJUSTL(dummyString)) // "F10.2)") outputMatrix(i, :)
     END DO
  END IF
  
#ifdef MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  CALL MPI_FINALIZE(ierr)
#endif
  
  STOP
END PROGRAM main
