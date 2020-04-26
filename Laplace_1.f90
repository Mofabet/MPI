program Laplace

implicit none
include 'mpif.h'
INTEGER :: ERR, SIZE, RANK, ST(MPI_STATUS_SIZE)
integer :: i, j, n, m, c_1, c_2, c_3, c_4, duck_1, disp, iter
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

  allocate(A(n,m))

  do i=1,n
      read(10,*)(A(i,j), j=1,m)
  enddo

  MBAND(:) = N1/SIZE

  do i = 1, mod(m, SIZE)
      disp(i) = disp(i)+1
  enddo

  MBAND(1) = disp(1) + 1
  MBAND(SIZE) = disp(SIZE) + 1

  do c_1 = 2, SIZE - 1
      MBAND(i) = disp(i) + 2
  enddo
endif


  CALL MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)
  CALL MPI_BCAST(m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)
  CALL MPI_BCAST(eps, 1, MPI_REAL, 0, MPI_COMM_WORLD, ERR)
  CALL MPI_BCAST(MBAND, SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ERR)

  if (RANK .eq. 0) then
    do c_1 = 1, SIZE - 1
      int = MBAND(i + 1)
      allocate(iBAND(MBAND(c_1+1),m1))
      do c_2 = 1, MBAND(c_1)!m     !do c_1 = 1,iBAND
        do c_3 = 1, m1
            iBAND(c_2,c_3) = E(disp(с_1 + 1) + c_2 - 1, c_3) !заполнение, c_1 i+j
        enddo
      enddo

      call MPI_SEND(iBAND, int*n,MPI_DOUBLE_PRECISION,c_1,20+c_1,MPI_COMM_WORLD,ERR)
      deallocate (iBAND)!---
      !duck_2 = duck_2 + 1
      !write(6,*),duck_2,RANK
  enddo !68 79

  int = MBAND(1)
  allocate(iBAND(n,int))
  do c_1 = 1, n
    do c_2 = 1, int
      iBAND(c_1,c_2) = A(c_1,c_2)
    enddo
  enddo

  !do
  !enddo
else
  int = MBAND(RANK + 1)
  CALL MPI_RECV(iBAND, int*n,MPI_DOUBLE_PRECISION,0,20+RANK,MPI_COMM_WORLD,ST,ERR)
  endif

  
      CALL MPI_SENDRECV(x(:,2), n, MPI_DOUBLE_PRECISION, left,(RANK)+(k*k+k)*99, &
      band(:,q), n, MPI_DOUBLE_PRECISION, right, &
      (RANK+1)+(k*k+k)*99, MPI_COMM_WORLD, ST, ERR)

      CALL MPI_SENDRECV(x(:,q-1), n, MPI_DOUBLE_PRECISION, right, (RANK)+(k*k+k)*909, &
      band(:,1), n, MPI_DOUBLE_PRECISION, left, (RANK-1)+(k*k+k)*909, &
      MPI_COMM_WORLD, ST, ERR)

      CALL MPI_ALLREDUCE(error, sumerr, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ERR)
      sumerr = sqrt(sumerr)
