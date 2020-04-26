program SLAU
implicit none
include 'mpif.h'
INTEGER :: ERR, SIZE, RANK, ST(MPI_STATUS_SIZE)
integer :: i, j, n1, n2, m1, m2, c_1, c_2, c_3, c_4, c_5, iter
integer, allocatable ::  MBAND(:), disp(:)
double precision, allocatable :: A(:,:), B(:,:), C(:,:), D(:,:), iBAND(:,:)
double precision, allocatable :: X(:), XOLD(:), XNEW(:), g(:), g_0(:)
double precision :: tmp, eps, err_0, error

call MPI_INIT(ERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE, ERR)
call MPI_COMM_RANK(MPI_COMM_WORLD, RANK, ERR)

    allocate(MBAND(SIZE),disp(SIZE))

    if (RANK .eq. 0) then

    open(10, file = 'A', form = 'formatted', status = 'unknown')
    read(10,*)n1
    read(10,*)m1
    endif
    CALL MPI_BCAST(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)
    CALL MPI_BCAST(m1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)
    allocate(XOLD(n1), XNEW(n1))


    if(RANK.eq.0) then
    allocate(A(n1,m1), C(n1,m1), D(n1,m1))

    do i=1,n1
    read(10,*)(A(i,j), j=1,n1)
    enddo

    write(6,*)'Matrix A = '
    do c_1 = 1, n1
    write(6,'(100f8.3)')(A(c_1, c_2), c_2 =1, m1)
    enddo

    open(10, file = 'B', form = 'formatted', status = 'unknown')
    read(10,*)n2
    read(10,*)m2
    endif

    CALL MPI_BCAST(n2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)
    CALL MPI_BCAST(m2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)

    if(RANK.eq.0) then
    allocate(B(n2,m2), g(n2))
    do i=1,n2
    read(10,*)(B(i,j), j=1,m2)
    enddo

    write(6,*)'Vector B = '
    do c_1 = 1, n2
    write(6,'(100f8.3)')(B(c_1, c_2), c_2 =1, m2)
    enddo

    C(:,:) = 0.d0
    do i = 1, n1
    C(i,i)=1/A(i,i)
    enddo
    do i = 1,n1
    do j = 1,m1
        C(i,j)=A(i,j)/A(i,i)
        g(i)=B(i,1)/A(i,i)
    enddo
    enddo
    XOLD = g
    D = -C

    do i = 1, n1
    D(i,i) = D(i,i)+1.d0
    enddo

    write(6,*)'Matrix B = '
    do c_1 = 1, n1
    write(6,'(100f8.3)')(D(c_1, c_2), c_2 =1, m1)
    enddo

    write(6,*)'Vector G = '
    do c_1 = 1, n2
    write(6,'(100f8.3)')g(c_1)
    enddo


    MBAND(:) = n1/SIZE

    do i = 1, mod(n1, SIZE)
        MBAND(i) = MBAND(i)+1
    enddo
    disp(1)=0
    do i =2,SIZE
        disp(i)=disp(i-1)+MBAND(i-1)
    enddo
    write (6,*)'Eps'
    read(*,*)Eps

    endif


    CALL MPI_BCAST(eps, 1, MPI_REAL, 0, MPI_COMM_WORLD, ERR)
    CALL MPI_BCAST(MBAND, SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ERR)
    CALL MPI_BCAST(disp, SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ERR)
    if (RANK .eq. 0) then

    do c_1 = 1,SIZE-1
    allocate(iBAND(MBAND(c_1+1),m1), g_0(MBAND(c_1+1)))
    do c_2 = 1, MBAND(c_1+1)
    g_0(c_2) = g(c_2+disp(c_1 + 1))
    do c_3 = 1, m1
          iBAND(c_2,c_3) = D(c_2+disp(c_1+1), c_3)
    enddo
    enddo
    call MPI_SEND(iBAND, MBAND(c_1+1)*m1, MPI_DOUBLE_PRECISION, c_1, c_1*10, MPI_COMM_WORLD, ERR)
    call MPI_SEND(g_0, MBAND(c_1+1) , MPI_DOUBLE_PRECISION, c_1, c_1*11, MPI_COMM_WORLD, ERR)
!    write(6,*)'Sent to proc #', c_1
!    write(6,*)'iBand for proc #', c_1, '='
!    do c_4 = 1, MBAND(c_1+1)
!    write(6, '(100f8.3)')(iBAND(c_4, c_5), c_5=1,m1)
!    enddo

!    write(6,*)'g_0 for proc #', c_1, '='
!    do c_4 = 1, MBAND(c_1+1)
!    write(6, '(100f8.3)')(g_0(c_4))
!    enddo

    deallocate (iBAND, g_0)

    enddo

    endif


    if(RANK.eq.0) then
    allocate(iBAND(MBAND(1),m1), g_0(MBAND(1)))
    do c_2 = 1, MBAND(1)
    g_0(c_2) = g(c_2)
    do c_3 = 1, m1
          iBAND(c_2,c_3) = D(c_2, c_3)
    enddo
    enddo

!    write(6,*)'iBand for proc # 0 ='
!    do c_4 = 1, MBAND(1)
!    write(6, '(100f8.3)')(iBAND(c_4, c_5), c_5=1,m1)
!    enddo

!    write(6,*)'g_0 for proc # 0 ='
!    do c_4 = 1, MBAND(1)
!    write(6, '(100f8.3)')(g_0(c_4))
!    enddo


    else
    allocate(g_0(MBAND(RANK+1)))
    allocate(iBAND(MBAND(RANK+1),m1))
    CALL MPI_RECV(g_0, MBAND(RANK+1), MPI_DOUBLE_PRECISION, 0, RANK*11, MPI_COMM_WORLD, ST, ERR)
    CALL MPI_RECV(iBAND, MBAND(RANK+1)*m1, MPI_DOUBLE_PRECISION, 0, RANK*10, MPI_COMM_WORLD, ST, ERR)
    endif

    do iter = 1,20

    allocate(X(MBAND(RANK+1)))
    CALL MPI_BCAST(XOLD,  n1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ERR)
    err_0 = 0
    error = 0


    do c_1 = 1 , MBAND(RANK+1)
        tmp = 0.d0

        do c_2 = 1,m1
            tmp = tmp + iBAND(c_1,c_2)*XOLD(c_2)
        enddo
        X(c_1) = tmp + g_0(c_1)
    enddo
    !write(6,*)'iBand for proc #', RANK, '=', X



    do c_1 = 1 , MBAND(RANK+1)
        err_0 = err_0 + (X(c_1)-XOLD(c_1 +disp(RANK+1)))**2
    enddo


    call MPI_ALLREDUCE(err_0, error, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,ERR)
    error = sqrt(error)


    if (RANK .eq. 0) then
    write(6,*)'ERROR = ', error
    do c_1 = 1, MBAND(1)
    XOLD(c_1) = X(c_1)
    enddo
    deallocate(X)
    else
    CALL MPI_SEND(X, MBAND(RANK+1), MPI_DOUBLE_PRECISION, 0, RANK*1000+iter, MPI_COMM_WORLD, ERR)
    deallocate(X)
    endif

    if (RANK.eq. 0) then
    do c_1 = 1, SIZE-1
    allocate(X(MBAND(c_1+1)))
    CALL MPI_RECV(X, MBAND(c_1+1), MPI_DOUBLE_PRECISION, c_1, c_1*1000+iter, MPI_COMM_WORLD, ST, ERR)
    do c_2 = 1, MBAND(c_1 +1)
    XOLD(c_2+disp(c_1+1)) = X(c_2)
    enddo
    deallocate(X)
    enddo
    endif



    !write(6,*)'Error', error

    if(error.lt.eps) goto 1000
    enddo


    if (RANK .eq. 0) then
    write(6,*)'Error was reached = ', error
    do c_1 = 1, n1
    write(6,*)'X = ',c_1,' = ',XOLD(c_1)
    enddo
    endif


1000    call MPI_FINALIZE(ERR)




    end
