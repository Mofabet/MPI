program Laplace

  implicit none
  include 'mpif.h'
  INTEGER :: ERR, SIZE, RANK, ST(MPI_STATUS_SIZE)
  integer :: i, j, n, m, s, c_1, c_2, c_3, duck_1, k, top, bottom, iter, iterrations, NE_INT, mrbin
  integer, allocatable ::  MBAND(:), row(:), disp(:)
  real :: eps, err_0, error, dt1, dt2
  double precision, allocatable :: A(:,:), B(:,:), C(:,:), D(:,:), E(:,:), iBAND(:,:)
  double precision, allocatable :: tm(:,:), t(:,:), out(:,:), g_0(:)
  double precision :: tmp

call MPI_INIT(ERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE, ERR)
call MPI_COMM_RANK(MPI_COMM_WORLD, RANK, ERR)

if (RANK .eq. 0) then !0000000000
  write (6,*)'Initialization of the program'
  write (6,*)'Enter the required accuracy:'
  read(*,*)eps
!гениратор
  write(6,*) 'Input number of lines:'
  read(*,*)n
  write(6,*) 'Input number of columns:'
  read(*,*)m
  allocate(tm(n,m))
  tm(:,:)=0.d0
  !чочтавим матрицу темп
  write(6,*) 'Input temperatures on 4 corners'        !(1:m)
  write(6,*) 'Upper left:'                            !(n:1)
  read(*,*) tm(1,1)
  write(6,*) 'Upper right:'
  read(*,*) tm(1,m)
  write(6,*) 'Lower left:'
  read(*,*) tm(n,1)
  write(6,*) 'Lower right:'
  read(*,*) tm(n,m)

  write(6,*) 'Number of iterrations ='
  read(*,*) iterrations
write(6,*)'--------------------------------------------------------------------'
endif !0000000000
  !open(10, file = 'A', form = 'formatted', status = 'unknown')
  !read(10,*)n1
  !read(10,*)m1

  allocate(row(SIZE))

!  do i=1,n
!      read(10,*)(A(i,j), j=1,m)
!  enddo
if (rank.eq.0) then !0000000000

dt1=(tm(1,m)-tm(1,1))/m
dt2=(tm(n,m)-tm(n,1))/m

do i=2, m-1
  tm(1,i)=tm(1,1)+i*dt1
  tm(n,i)=tm(n,1)+i*dt2
enddo

dt1=(tm(n,1)-tm(1,1))/n
dt2=(tm(n,m)-tm(1,m))/n

do i=2, n-1
  tm(i,1)=tm(1,1)+i*dt1
  tm(i,m)=tm(1,m)+i*dt2
enddo

write(6,*)'Your list is'
do c_1 = 1, n
  write(6,'(100f8.3)')(tm(c_1,c_2), c_2=1,m)
enddo

!need alloc
allocate(disp(SIZE))
!allocate(array, stat=err)
!if ( /= 0) print *, ": Allocation request denied"
!
!if (allocated()) deallocate(, stat=)
!if ( /= 0) print *, ": Deallocation request denied"
  disp(:) = N/SIZE

  do c_1 = 1, mod(n, SIZE)
      disp(c_1) = disp(c_1)+1
  enddo

  row(:) = disp(:) + 1
  !MBAND(SIZE) = disp(SIZE) + 1

  do c_2 = 2, SIZE - 1
      row(c_2) = row(c_2) + 1
  enddo

  write(6,*)'disp = ', disp

!  write(6,*)'row = ', row
!-------------------------------------B----------------------------------------!
duck_1 = disp(1)-2
    do c_1 = 1, SIZE - 1
      NE_INT = row(c_1 + 1)
      allocate(iBAND(NE_INT,m))
      do c_2 = 1, NE_INT!m     !do c_1 = 1,iBAND
        do c_3 = 1, m
            iBAND(c_2,c_3) = tm(c_1+c_2+ duck_1, c_3) !заполнение, c_1 i+j !!!!!!!!!!!!
        enddo
      enddo
!2 3 svobodni
duck_1 = duck_1 + disp(c_1 + 1) - 1
      call MPI_SEND(iBAND, NE_INT*m,MPI_DOUBLE_PRECISION,c_1,20+c_1,MPI_COMM_WORLD,ERR)!mmm*c
      write(6,*)'Band on ',c_1,'=='
      do c_2=1, NE_INT
        write(6,'(100f8.3)')(iBAND(c_2,c_3),c_3 =1,m)
      enddo
      deallocate (iBAND)!---
      !duck_2 = duck_2 + 1
      !write(6,*),duck_2,RANK
  enddo !68 79
endif  !0000000000
!if (RANK .ne. 0) then
!  allocate(row(SIZE))
!endif

CALL MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)
CALL MPI_BCAST(m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)
CALL MPI_BCAST(eps, 1, MPI_REAL, 0, MPI_COMM_WORLD, ERR)
CALL MPI_BCAST(iterrations,1, MPI_INTEGER, 0, MPI_COMM_WORLD, ERR)
CALL MPI_BCAST(row, SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
!ALL MPI_BCAST(MBAND, SIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ERR)
!allocate(row(SIZE))

if (RANK .eq. 0) then  !0000000000
  NE_INT = row(1)
  allocate(iBAND(NE_INT,m)) !tut tozhe ne ponyatno nm NE_INT 1 ili 2 s por-m...
  do c_1 = 1, NE_INT
    do c_2 = 1, m
      iBAND(c_1,c_2) = tm(c_1,c_2)
    enddo
  enddo  !0000000000


  !do
  !enddo
else !---------------- !1231231313131313131
  NE_INT = row(RANK + 1)
  allocate(iBAND(NE_INT,m), t(NE_INT,m))
  CALL MPI_RECV(iBAND, NE_INT*m,MPI_DOUBLE_PRECISION,0,20+RANK,MPI_COMM_WORLD,ST,ERR) !т/т
  endif !1231231313131313131

  !http://www2.sscc.ru/Publikacii/Primery_Prll/1-4.htm
!SDVIG
  if (RANK .eq. 0) then !0000000000
    top = MPI_PROC_NULL                            !VERH
  else
    top = RANK - 1
  endif

  if (RANK .eq. SIZE - 1) then
    bottom = MPI_PROC_NULL                         !NIZ-KARNIZ
  else
    bottom = RANK + 1
  endif

  T = iBAND

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc                                                                        !ccc

do iter = 1, iterrations
err_0 = 0
error = 0

NE_INT = row(RANK + 1)

do c_1 = 2, (NE_INT - 1)
  do c_2 = 2, (m - 1)          !  1                    2                      3                    4
    T(c_1,c_2) = ((iBAND(c_1 - 1,c_2) + iBAND(c_1 + 1,c_2) + iBAND(c_1,c_2 - 1) + iBAND(c_1,c_2 + 1))/4)
    err_0 = err_0 + (t(c_1,c_2) - iBAND(c_1,c_2))**2
  enddo
enddo

iBAND = T

CALL MPI_SENDRECV(T(2,:), m, MPI_DOUBLE_PRECISION, top, RANK+(iter*(iter + 100)+iter), & !тег надо заменить, если лимит будет привышен
iBAND(NE_INT, :), m, MPI_DOUBLE_PRECISION, bottom, (RANK + 1)+(iter*(iter + 100)+iter), & !можно повысить степень ранка, а ост оставить
MPI_COMM_WORLD, ST, ERR)
CALL MPI_SENDRECV(T(NE_INT - 1,:), m, MPI_DOUBLE_PRECISION, bottom,  RANK+(iter*(iter+999)+iter), &
IBAND(1,:), m,MPI_DOUBLE_PRECISION,top,  (RANK - 1)+(iter*(iter + 999)+iter), &
MPI_COMM_WORLD, ST, ERR)

CALL MPI_ALLREDUCE(err_0, error, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ERR)
error=sqrt(error)

  if(error .lt. eps) goto 2000 !почему экзит не робит?
enddo !и что? 164 строка, но просит завершения
!ccc                                                                        !ccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
2000 continue
write(6,*)'Band on'
write(6,*) RANK,' = ', iBAND

if (RANK .eq. 0) then
  allocate(out(n,m))
  do c_1 = 1, NE_INT - 1 !perenos
    do c_2 = 1, m
      out(c_1,c_2) = iBAND (c_1,c_2)
    enddo
  enddo
deallocate(iBAND)
else
  call MPI_SEND(iBAND, NE_INT*m, MPI_DOUBLE_PRECISION, 0, 10+RANK, MPI_COMM_WORLD, ERR)
endif

if (RANK .eq. 0) then
  duck_1 = disp (1)
  do c_1 = 1,SIZE-1
    if (c_1 .eq. (SIZE - 1)) then
      mrbin = 0
    else
      mrbin = 1
    endif
    NE_INT = row (c_1 + 1)
    allocate(iBAND(NE_INT,m))
    CALL MPI_RECV(iBAND, NE_INT*m, MPI_DOUBLE_PRECISION, c_1, 10+c_1, MPI_COMM_WORLD, ST, ERR)
    do c_2 = 2, (NE_INT - mrbin) !!!!!!!!!!!
      do c_3 = 1, m
        out(c_2 + duck_1 - 1, c_3) = iBAND(c_2,c_3)
      enddo
    enddo
    duck_1 = duck_1 + disp(c_1 + 1)
    deallocate(iBAND)
  enddo



!if (RANK .eq. 0) then
write(6,*)'--------------------------------------------------------------------'
write(6,*)'ERROR = ', error
        !    write(6,*)'Max number of iter was reached = ', error
                do c_1 = 1, n
                   write(6, '(100f8.3)')(out(c_1,c_2), c_2=1,m)
                enddo
          endif
call MPI_FINALIZE(error)


end
