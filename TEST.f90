program TEST
implicit none
include 'mpif.h'
INTEGER :: ERR, SIZE, RANK, ST(MPI_STATUS_SIZE)
integer :: i, j, n1, n2, m1, m2, c_1, c_2, c_3, c_4, duck_1, duck_2, iter
integer, allocatable ::  MBAND(:), disp(:)
real :: time, time_1, time_2
double precision, allocatable :: A(:,:), B(:,:), C(:,:), D(:,:), E(:,:), iBAND(:,:)
double precision, allocatable :: X(:), X_0(:), X_1(:), g(:), g_0(:), tmp_alloc(:)
double precision :: tmp, eps, err_0, error

!  3. Прочитать на «0» процессе из файла А массив NxN (размерность массива N записана в первой строке файла А).
!  Переслать массив с «0» процесса на «1» процесс.
!  На первом процессе вывести на экран полученный массив.
!  Измерить время выполнения параллельной части программы.


call MPI_INIT(ERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE, ERR)
call MPI_COMM_RANK(MPI_COMM_WORLD, RANK, ERR)

if (RANK .eq. 0)  then
time_1 = MPI_WTIME(ERR)
open(10, file = 'A', form = 'formatted', status = 'unknown')
read(10,*)n1
read(10,*)m1

allocate(A(n1,m1), B(n1,m1))

do i=1,n1
    read(10,*)(A(i,j), j=1,n1)
enddo

call MPI_SEND(A, n1*m1,  MPI_DOUBLE_PRECISION, 1, 11, MPI_COMM_WORLD, ERR)
deallocate(A)
endif

call MPI_BCAST(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)
call MPI_BCAST(m1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)

if (RANK .eq. 1) then
allocate(B(n1,m1))
call MPI_RECV(B, n1*m1, MPI_DOUBLE_PRECISION, 0, 11, MPI_COMM_WORLD, ST,ERR)

write(6,*)'Matrix A = '
do c_1 = 1, n1
  write(6,'(100f8.3)')(B(c_1, c_2), c_2 =1, m1)
enddo
deallocate(B)
endif

if (RANK .eq. 0)  then
time_2 = MPI_WTIME(ERR)
time   = time_2 - time_1
write(6,*)'TIME = ', time
endif

call MPI_FINALIZE(ERR)
end
