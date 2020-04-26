program SLAU
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

allocate(MBAND(SIZE),disp(SIZE))

if (RANK .eq. 0) then
  write (6,*)'Eps'
  read(*,*)Eps

open(10, file = 'A', form = 'formatted', status = 'unknown')
read(10,*)n1
read(10,*)m1

allocate(A(n1,m1), C(n1,m1), D(n1,m1)) !Core_0
allocate(X_0(n1),X(n1))

do i=1,n1
    read(10,*)(A(i,j), j=1,n1)
enddo

open(10, file = 'B', form = 'formatted', status = 'unknown')
read(10,*)n2 !=n1
read(10,*)m2 !=1

!allocate(A(n1,m1), B(n2,m2))
  !???
  allocate(B(n2,m2), g(n2))

  do i=1,n2
      read(10,*)(B(i,j), j=1,m2)
  enddo

  ! A*x=b
  !x=bx+g
  C(:,:) = 0.d0
  do i = 1, n1
    C(i,i)=1/A(i,i)
  enddo

  do i = 1,n1
    do j = 1,m1
          C(i,j)=A(i,j)/A(i,i) !poluchaem obrat-yu
          g(i)=B(i,1)/A(i,i)
      enddo
  enddo

  D= -C
  X_0 = g

  do i = 1, n1
    !do j = 1, m1
      D(i,i) = D(i,i) + 1.d0
    !X_0
  !enddo
  enddo
!allocate(MBAND(SIZE))-----------------------
  MBAND(:) = N1/SIZE

  do i = 1, mod(n1, SIZE)
      MBAND(i) = MBAND(i) + 1
  enddo
  disp(1) = 0
  do c_1 =2,SIZE
      disp(c_1) = disp(c_1 - 1) + MBAND(c_1 - 1)
  enddo
!  disp = MBAND(1)-1 !пока 1
endif ! 16

!!!
CALL MPI_BCAST(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)
CALL MPI_BCAST(m1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)
CALL MPI_BCAST(n2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)
CALL MPI_BCAST(m2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)
CALL MPI_BCAST(eps, 1, MPI_REAL, 0, MPI_COMM_WORLD, ERR)
CALL MPI_BCAST(MBAND, SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ERR)
!!!
 !---
!allocate(g(n1),X_0(n1))
!allocate(g_0(MBAND(RANK+1)))

if (RANK .eq. 0) then

  !  iBAND = MBAND(i+1)
!--------------------------------------------------------------
  do c_1 = 1,SIZE-1
    allocate(iBAND(MBAND(c_1 + 1),m1),g_0(MBAND(c_1 + 1)))
    do c_2 = 1, MBAND(c_1 + 1)!m     !do c_1 = 1,iBAND
      g_0(c_2) = g(c_2 + disp(c_1 + 1))
      do c_3 = 1, m1
          iBAND(c_2,c_3) = D(disp(c_1 + 1) + c_2, c_3) !заполнение, c_1 i+j
      enddo
    enddo
      !send& и нужен деалок
      call MPI_SEND(iBAND, MBAND(c_1 + 1)*m1,MPI_DOUBLE_PRECISION,c_1,20+RANK+c_1,MPI_COMM_WORLD,ERR)
      call MPI_SEND(g_0, MBAND(c_1 + 1),MPI_DOUBLE_PRECISION, c_1, 30+RANK+c_1, MPI_COMM_WORLD, ERR) !-----
      deallocate (iBAND, g_0)!---
      !duck_2 = duck_2 + 1
      !write(6,*),duck_2,RANK
  enddo !68 79
endif
  !allocate(iBAND(:,:)) !0
!--------------------------------------------------------------
  !  do c_1 = 1, SIZE-1
  !    allocate(g_0(MBAND(c_1)))
  !    do c_2 = 1, MBAND(c_1)
  !      g_0(c_2) = g(c_1+disp(c_1) + 1)  !VEKTOR G
  !    enddo
  !      deallocate(g_0)                                                                                !|
  !      !                                                                                              !|
  !  enddo!!!!!!!!!!!!!!! теперь нет айбенда и g                                                        !|
    !allocate(g_0(:,:))   !0                                                                           !|
 !--------------------------------------------------------------                                       !|
 !deallocate(?)
 if(RANK.eq.0) then
 allocate(iBAND(MBAND(1),m1), g_0(MBAND(1)))
 do c_2 = 1, MBAND(1)
 g_0(c_2) = g(c_2)
 do c_3 = 1, m1
       iBAND(c_2,c_3) = D(c_2, c_3)
 enddo
 enddo
                                                                                       !|
else !65 74
  allocate(X_0(n1),X(n1))                                                                                          !|
  allocate(g_0(MBAND(RANK+1))) !COL                                                                    !|
  allocate(iBAND(MBAND(RANK+1),m1))          ! ne 0                                                    !|
    CALL MPI_RECV(iBAND,MBAND(RANK+1)*m1,MPI_DOUBLE_PRECISION,0,20+RANK,MPI_COMM_WORLD,ST,ERR)       !|
    CALL MPI_RECV(g_0,MBAND(RANK+1),MPI_DOUBLE_PRECISION,0,30+RANK,MPI_COMM_WORLD,ST,ERR)          !<-----
endif !65

!all 2
          !  =0.d0
          !alloc?
!allocate(X_0(:,:)) !все?
!allocate(X_0(m1))
!X_0=0.d0
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc                                                                        !ccc

  do iter = 1,20 !raws   -----ITER C
    allocate(X_1(MBAND(RANK+1)))
    CALL MPI_BCAST(X_0, n1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ERR)
    err_0 = 0
    error = 0

   do c_1 = 1 , MBAND(RANK+1)
      tmp = 0.d0
      do c_2 = 1,m1
         tmp = tmp + iBAND(c_1,c_2)*X_0(c_2)
      enddo
      X_1(c_1) = tmp + g_0(c_1)
      !write(6,*),X(c_1),RANK
   enddo !+
   !write(6,*),X(c_1),RANK
   !-----------------ERR

  !duck_1 = duck_1 + 1
  !write(6,*),duck_1,RANK
  do c_1 = 1 , MBAND(RANK+1)
     err_0 = err_0 + ((X_1(c_1)-X_0(c_1 + disp(RANK+1)))**2)
  enddo
!MPI_REDUCE(SBUF, RBUF, COUNT, DATATYPE, OP, ROOT, COMM, IERR)
!C-C-C
  call MPI_ALLREDUCE(err_0, error, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,ERR)
  error = sqrt(error)

  if(RANK .eq. 0) then
    write(6,*)'ERROR = ', error
    do c_1 = 1, MBAND(1)
      X_0(c_1) = X_1(c_1)
    enddo
    deallocate(X_1)
  else
    call MPI_SEND(X_1, MBAND(RANK + 1),MPI_DOUBLE_PRECISION,0,40 + RANK + 100*iter,MPI_COMM_WORLD,ERR)
    deallocate(X_1)
  endif

  if(RANK .eq. 0) then
    do c_2 = 1, SIZE - 1
      allocate(X_1(MBAND(c_2 + 1)))
      call MPI_RECV(X_1, MBAND(c_2 + 1),MPI_DOUBLE_PRECISION,c_2,40 + c_2 + 100*iter,MPI_COMM_WORLD, ST, ERR)

    !error = sqrt(error)
    do c_3 = 1, MBAND(c_2 + 1) !???
       X_0(c_3 + disp(c_2 + 1)) = X_1(c_3)
    enddo
    deallocate(X_1)
  enddo
  do c_4 = 1, n1
   X(c_4) = X_0(c_4)
  enddo



!    do c_3 = 1, SIZE-1
!       allocate(tmp_alloc(MBAND(c_3)))
!             do c_4 = 1,MBAND(c_3)
!               X(c_4+disp(c_3))=tmp_alloc(c_4)
!            enddo
!       deallocate(tmp_alloc)
!    enddo
!C-C-C

  endif !131 137
!  call MPI_BCAST(X,m1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ERR)
!  call MPI_BCAST(error,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ERR)
  if(error .lt. eps) goto 1000
enddo
!ccc                                                                        !ccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

if (RANK .eq. 0) then
            write(6,*)'Max number of iter was reached = ', error
                do c_1 = 1, n1
                   write(6,*)'X = ',c_1,' = ',X(c_1)
                enddo
          endif
         goto 1001


1000     if (RANK .eq. 0) then
            write(6,*)'accuracy was reached = ', error
                do c_1 = 1, n1
                   write(6,*)'X = ',c_1,' = ',X(c_1)
                enddo
          endif

1001         call MPI_FINALIZE(ERR)
end
