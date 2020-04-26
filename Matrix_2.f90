program Matrix_1
implicit none
include 'mpif.h'
INTEGER :: ERR, SIZE, RANK, ST(MPI_STATUS_SIZE)
integer :: n, m, i, j, hlf, z
double precision, allocatable :: A(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:), AT(:,:), ATI(:,:), T(:,:)
double precision :: tmp

call MPI_INIT(err)
call MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE, ERR)
call MPI_COMM_RANK(MPI_COMM_WORLD, RANK, ERR)

if (RANK .eq. 0) then
open(10, file = 'A', form = 'formatted', status = 'unknown')
read(10,*)n
read(10,*)n

allocate (A(n,n))
hlf = n/2
Z = hlf*hlf
do i = 1,n
    read(10,*)(A(i,j), j=1,n)
enddo

write(6,*)'Matrix A = '
do i=1,n
write(6,'(100f8.1)')( A(i,j), j=1,n)
enddo

endif

call MPI_BCAST(hlf, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)
call MPI_BCAST(Z, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)

allocate(A1(hlf,hlf), A2(hlf,hlf), A3(hlf,hlf), A4(hlf,hlf), AT(hlf,hlf), T(hlf,hlf))

if(rank .eq. 0) then
do i = 1, hlf
  do j = 1, hlf
    A1(i,j) = A(i, j)
    A2(i,j) = A(i+hlf, j)
    A3(i,j) = A(i, j+hlf)
    A4(i,j) = A(i+hlf, j+hlf)
  enddo
enddo

call MPI_SEND(A2, Z,  MPI_DOUBLE_PRECISION, 1, 11, MPI_COMM_WORLD, ERR)
write(6,*)'CORE_0_SEND1: IM OK'
call MPI_SEND(A3, Z,  MPI_DOUBLE_PRECISION, 2, 12, MPI_COMM_WORLD, ERR)
write(6,*)'CORE_0_SEND2: IM OK'
call MPI_SEND(A4, Z,  MPI_DOUBLE_PRECISION, 3, 13, MPI_COMM_WORLD, ERR)
write(6,*)'CORE_0_SEND3: IM OK'

A1 = transpose(A1)
!остальные проц
else
  write(6,*)'CORE_1/2/3_RECV0: IM OK',rank

  call MPI_RECV(T, Z, MPI_DOUBLE_PRECISION, 0, 10+rank, MPI_COMM_WORLD, ST,ERR)
  write(6,*)'CORE_1/2/3_RECV1: IM OK',rank
  endif

  write(6,*)'AFTERIF: IM OK',rank
  AT = transpose(T)
  write(6,*)'Transpose: IM OK',rank

  call MPI_SEND(AT, Z,  MPI_DOUBLE_PRECISION, 0, 20+rank, MPI_COMM_WORLD, ERR)
  write(6,*)'SEND4: IM OK'

  if(rank .eq. 0) then
  call MPI_RECV(A2, Z, MPI_DOUBLE_PRECISION, 1, 21, MPI_COMM_WORLD, ST,ERR)
  call MPI_RECV(A3, Z, MPI_DOUBLE_PRECISION, 2, 22, MPI_COMM_WORLD, ST,ERR)
  call MPI_RECV(A4, Z, MPI_DOUBLE_PRECISION, 3, 23, MPI_COMM_WORLD, ST,ERR)

  do i = 1, hlf
    do j = 1, hlf
      A(i,j) = A1(i,j)
      A(i+hlf,j) = A3(i,j)
      A(i,j+hlf) = A2(i,j)
      A(i+hlf,j+hlf) = A4(i,j)
    enddo
  enddo

  write(6,*)'Transpouse Matrix = '
  do i=1,n
  write(6,'(100f8.3)')( A(i,j),j=1,n)
  enddo
  write(6,*)'Matrix has benn transposed!'
  endif

      CALL MPI_FINALIZE(ERR)


      end
