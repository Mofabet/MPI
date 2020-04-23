program SLAU

implicit none
include 'mpif.h'
INTEGER :: ERR, SIZE, RANK, ST(MPI_STATUS_SIZE)
integer :: i, j, n1, n2, m1, m2, c_1, c_2, c_3, c_4, duck_1, disp, iter
double precision, allocatable :: A(:,:), B(:,:), C(:,:), D(:,:), E(:,:), MBAND(:,:), iBAND(:,:)
double precision, allocatable :: X(:,:), X_0(:,:), X_1(:,:), g(:,:), g_0(:,:)
double precision :: tmp, tmp_2, eps, err_0, error

call MPI_INIT(ERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE, ERR)
call MPI_COMM_RANK(MPI_COMM_WORLD, RANK, ERR)

allocate(A(:,:), B(:,:), C(:,:),D(:,:), E(:,:), MBAND(SIZE),iBAND(:,:))

  open(10, file = 'B', form = 'formatted', status = 'unknown')
  read(10,*)n2 !=n1
  read(10,*)m2 !=1
allocate(A(n1,m1), B(n2,m2))
  !???

  do i=1,n2
      read(10,*)(B(i,j), j=1)
  enddo

  ! A*x=b
  !x=bx+g
  do i = 1,n1
    do j = 1,m1
          C(i,j)=A(i,j)/A(i,i) !poluchaem obrat-yu
      enddo
  enddo

  D= -C

  do i = 1, n1
    do j = 1, m1
      E(i,i) = D(i,i)+1.d0
    !X_0
  enddo
  enddo

  MBAND(:) = N1/SIZE

  do i = 1, mod(n1, SIZE)
      MBAND(i) = MBAND(i)+1
  enddo
  disp(1)=0
  do i =2,SIZE
      disp(i)=disp(i-1)+MBAND(i-1)
  enddo
!  disp = MBAND(1)-1 !пока 1
endif ! 16

!!!
CALL MPI_BCAST(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)
CALL MPI_BCAST(m1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)
CALL MPI_BCAST(eps, 1, MPI_REAL, 0, MPI_COMM_WORLD, ERR)
CALL MPI_BCAST(MBAND, SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ERR)
!!!
 !---
allocate(g(:,:))
allocate(g_0(MBAND(RANK+1)))

if (RANK .eq. 0) then

  !  iBAND = MBAND(i+1)
!--------------------------------------------------------------
  do c_1 = 1,SIZE-1
    allocate(iBAND(MBAND(c_1),m1))
    do c_2 = 1, MBAND(c_1)!m     !do c_1 = 1,iBAND
      do c_3 = 1, m1
          iBAND(i,j) = E(c_1+disp) !заполнение, c_1 i+j
      enddo
    enddo
      !send&
      !call MPI_SEND(iBAND,MPI_DOUBLE_PRECISION,000000,,MPI_COMM_WORLD,ERR) !---
  enddo !68 79
!--------------------------------------------------------------
    do c_1 = 1, SIZE - 1
      do c_2 = 1, MBAND(c_1)
        g_0(c_2) = g(c_1+disp(c_1))  !VEKTOR G
      enddo
        call MPI_SEND(g_0, MBAND(c_1),MPI_DOUBLE_PRECISION, c_1, 30+RANK+c_1, MPI_COMM_WORLD, ERR) !-----
        deallocate(g_0)                                                                                !|
    enddo!!!!!!!!!!!!!!!                                                                               !|
 !--------------------------------------------------------------                                       !|
 !deallocate(?)                                                                                        !|
else !65 74                                                                                            !|
    CALL MPI_RECV(iBAND,BAND(RANK+1)*m1,MPI_DOUBLE_PRECISION,0,20+RANK,MPI_COMM_WORLD,ST,ERR)          !|
    CALL MPI_RECV(g_0,BAND(RANK+1),MPI_DOUBLE_PRECISION,0,30+RANK,MPI_COMM_WORLD,ST,ERR)          !<-----
endif !65

!all 2
          !  =0.d0
          !alloc?
allocate(X_0(:,:))
allocate(X_0(m1))
allocate(X_1(MBAND(RANK+1)))
X_0=0.d0

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc                                                                        !ccc

  do iteration = 1,100 !raws   -----ITER C
   duck_1 = 0.d0
   do c_1 = 1 , MBAND(RANK+1)
      tmp = 0.d0
      do c_2 = 1,m1
         tmp = tmp + MAND(c_1,c_2)*X_0(c_2)
      enddo
      X_1(c_1) = tmp + g_0(c_1)
      !write(6,*),X(c_1),RANK
   enddo !+
   !write(6,*),X(c_1),RANK
   !-----------------ERR

  err_0 = 0.d0
  !duck_1 = duck_1 + 1
  !write(6,*),duck_1,RANK
!C-C-C
  do c_1 = 1 , MBAND(RANK+1)
     err_0 = err_0 + ((X_1(c_1)-X_0(c_1))**2)
  enddo
!MPI_REDUCE(SBUF, RBUF, COUNT, DATATYPE, OP, ROOT, COMM, IERR)
  call MPI_REDUCE(err_0,error,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ERR)

  if(RANK .ne. 0) then
     call MPI_SEND(X_1,MBAND(RANK+1),MPI_DOUBLE_PRECISION,0,40*RANK,MPI_COMM_WORLD, ERR)
  else
    error = sqrt(error)
    do c_2 = 1,MBAND(RANK+1) !???
       X_0(c_2)=X_1(c_2)
    enddo

    do c_3 = 1, SIZE-1
       allocate(tmp_2(MBAND(c_3)))
       call MPI_RECV(tmp_2,MBAND,MPI_DOUBLE_PRECISION,c_3,40*RANK,MPI_COMM_WORLD, ERR)
             do c_4 = 1,MBAND(c_3)
               X(c_4+disp(c_3))=tmp_2(c_4)
            enddo
       deallocate(tmp_2)
    enddo
!C-C-C
  endif !131 137
  call MPI_BCAST(X,m1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ERR)
  call MPI_BCAST(error,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ERR)
  if(error .lt. eps) then
   goto 1000
  endif
       enddo

!ccc                                                                        !ccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

1000        if (RANK .eq. 0) then
                do c_1 = 1, n1
                 write(6,*)'X = ',c_1,' = ',X_0(c_1)
                enddo
              endif
         call MPI_FINALIZE(ERR)
         end
