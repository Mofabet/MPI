program Laplace

implicit none
include 'mpif.h'
INTEGER :: ERR, SIZE, RANK, ST(MPI_STATUS_SIZE)
integer :: i, j, n1, n2, m1, m2, c_1, c_2, c_3, c_4, duck_1, disp, iter
double precision, allocatable :: A(:,:), B(:,:), C(:,:), D(:,:), E(:,:), MBAND(:,:), iBAND(:,:), X(:,:), X_0(:,:), X_1(:,:), g(:,:), g_0(:,:)
double precision :: tmp, tmp_2, eps, err_0, error,

call MPI_INIT(ERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE, ERR)
call MPI_COMM_RANK(MPI_COMM_WORLD, RANK, ERR)

if (RANK .eq. 0) then
  write(6,*)'Точность'
  read(*,*)eps

  open(10, file = 'A', form = 'formatted', status = 'unknown')
  read(10,*)n1
  read(10,*)m1

  do i=1,n1
      read(10,*)(A(i,j), j=1,m1)
  enddo
