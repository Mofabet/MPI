program TEST
implicit none
include 'mpif.h'
INTEGER :: ERR, SIZE, RANK, ST(MPI_STATUS_SIZE)
integer :: i, j, n1, n2, m1, m2, c_1, c_2, c_3, c_4, duck_1, duck_2, iter
integer, allocatable ::  MBAND(:), disp(:)
double precision, allocatable :: A(:,:), B(:,:), C(:,:), D(:,:), E(:,:), iBAND(:,:)
double precision, allocatable :: X(:), X_0(:), X_1(:), g(:), g_0(:), tmp_alloc(:)
double precision :: tmp, eps, err_0, error

call MPI_INIT(ERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE, ERR)
call MPI_COMM_RANK(MPI_COMM_WORLD, RANK, ERR)
