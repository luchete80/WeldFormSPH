module Neighbor
!use Domain

implicit none 

! Neighbor data
integer, dimension(3):: cellno
real, dimension(3):: blpf, trpr, cellsize
real :: rhomax, cellfac, hmax
integer, allocatable, dimension(:,:,:) :: HOC, pairs_t! Head of chain, pairs (proc, first, second)
integer, allocatable, dimension(:,:) :: cc !cell 
integer, allocatable, dimension(:) :: LL, pair_count, nbcount_t    !Linked list, paircount (per thread)
logical :: nballoc_pass !Pass only
integer, allocatable, dimension (:,:) :: Neib_t !Filled after pairs

contains

  !! ALLOCATE THINGS
  subroutine InitNb() !Kernel Function parameter to check cellfac
    use Domain
    implicit none
    allocate (cc(part_count,3))
    allocate (ll(part_count))
    allocate (pair_count(Nproc))
    allocate (nbcount_t(part_count))
    cellfac = 2.0
    hmax = pt%h(1)
    write(*,*) "hmax ", hmax
  end subroutine InitNb

  subroutine CellInitiate ()
    use Domain
    implicit none
    !integer, intent(in)::part_count
    !integer, dimension(3), intent(in)::cellsize
    integer:: i,p
    blpf(:) = pt%x(1,:)
    
    write (*,*) "Looking for particles ", part_count
    do p=1, part_count
      do i=1,Dim
        !write (*,*) "i ",i, " x(p,i) ",pt%x(p,i), "trpr ", trpr(i)
        if ( pt%x(p,i) > trpr(i) ) then 
          trpr(i) = pt%x(p,i)
        end if
        if ( pt%x(p,i) < blpf(i) ) then 
          blpf(i) = pt%x(p,i)
        end if
        if ( pt%rho(p) > rhomax ) then 
          rhomax = pt%rho(p)
        end if
      end do
    end do
    write (*,*) "trpr ", trpr(:)
    write (*,*) "blpf ", blpf(:)
    ! Override the calculated domain size
    
    write (*,*) "Done, allocating cellno "
    !Calculate Cell properties
    do i = 1, Dim
      if (real(ceiling((trpr(i)-blpf(i))/(cellfac*hmax)-(trpr(i)-blpf(i))/(cellfac*hmax)))<hmax/10.) then
        cellno(i) = (trpr(i)-blpf(i))/(cellfac*hmax)
      else 
        cellno(i) = floor((trpr(i)-blpf(i))/(cellfac*hmax))
      end if
    end do
    
    cellsize(:) = TRPR(:)-BLPF(:)/real(CellNo(:))

    write (*, *) "Cell No ", cellno(:)
    !write (*, *) "Cell Size ", cellsize(:)
    !Periodic BC correction

    !Initiate Head of Chain array for Linked-List
    allocate (HOC(cellno(1),cellno(2),cellno(3)))
    HOC(:,:,:) = -1
   
  end subroutine CellInitiate 

  subroutine ListGenerate()
    use Domain
    implicit none
  !  integer, intent(in)::part_count
  !  integer, dimension(3), intent(in)::cellsize
    integer, dimension(3)::ijk
    integer ::a, d, temp
    
    do a = 1, part_count 
      !write (*,*) "Part ",a
      do d = 1, Dim
        !write (*,*) "blpf(d) ",blpf(d), "cellsize ", cellsize(d)
        ijk = int(floor(pt%x(a,d)-blpf(d))/cellsize(d))
        if (ijk(d) < 0) then
          ijk(d) = 0
        else if (ijk(d) >= cellno(d)) then
          ijk(d)=CellNo(d)-1;
        end if  
      end do
      
      temp = HOC(ijk(1),ijk(2),ijk(3))
      HOC(ijk(1),ijk(2),ijk(3)) = a;
      ll(a) = temp
      cc(a,:) = ijk(:)
    end do 

  end subroutine ListGenerate

  subroutine CellReset()
    use Domain, only : part_count
    integer::a
    do a=1, part_count
    ll(a) = -1
    end do
    HOC(:,:,:) = -1 
    
  end subroutine CellReset

  subroutine YZPlaneCellsNeighbourSearch(q1)
    use omp_lib
    implicit none
    integer, intent(in)  :: q1
    integer :: q2, q3, temp1, temp2, t, i, j
    t = omp_get_thread_num()
    do q3 = 1, CellNo(3)
      do q2 = 1, CellNo(2)
        if (HOC(q1,q2,q3)== -1) then
        else
          temp1 = HOC(q1,q2,q3)
          do while (temp1 .ne. -1)
            temp2 = LL(temp1)
            do while (temp2 .ne. -1)
              !AllocateNbPair(temp1,temp2,t)
              temp2 = LL(temp2)
            end do 
          
            !(q1 + 1, q2 , q3)
            if (q1+1 < CellNo(1)) then !!Original: if (q1+1 < CellNo(0)) 
              temp2 = HOC(q1+1,q2,q3)
              do while (temp2 .ne. -1) 
                  !AllocateNbPair(temp1,temp2,T);
                temp2 = LL(temp2)
              end do !while temp2!=-1
            end if ! (q1 + 1, q2 , q3)   

            !(q1 + a, q2 + 1, q3) & a(-1,1)
            if (q2+1< CellNo(2)) then
              do i = q1 - 1, q1 + 1 !for (int i = q1-1; i <= q1+1; i++) 
                if (i<CellNo(1) .and. i>=1) then !ORIG if (i<CellNo(0) .and. i>=0) 
                  temp2 = HOC(i,q2+1,q3)
                  do while (temp2 .ne. -1)
                  
                    !AllocateNbPair(temp1,temp2,T);
                    temp2 = LL(temp2)
                  end do
                end if
              end do
            end if

            !(q1 + a, q2 + b, q3 + 1) & a,b(-1,1) => all 9 cells above the current cell
            if (q3+1< CellNo(2)) then
              do j = q2-1, q2 + 1 ! for (int j=q2-1; j<=q2+1; j++)
                do i = q1 - 1, q1 + 1 ! for (int i=q1-1; i<=q1+1; i++) {
                  if (i<CellNo(1) .and. i>=1 .and. j<CellNo(2) .and. j>=1) then ! if (i<CellNo(0) && i>=0 && j<CellNo(1) && j>=0) {
                    temp2 = HOC(i,j,q3+1)
                    do while (temp2 .ne. -1)
                      ! AllocateNbPair(temp1,temp2,T);
                      temp2 = LL(temp2)
                    end do
                  end if
                end do
              end do !j
            end if
            temp1 = LL(temp1)        
          end do !temp1 .ne. -1 
        end if
      end do
    end do
  end subroutine YZPlaneCellsNeighbourSearch

  logical function CheckRadius(temp1, temp2)
    use Domain
    implicit none
    integer, intent (in) :: temp1,temp2
    real, dimension(3) :: xij
    real :: h
    CheckRadius = .false.
    h = 0.5*(pt%h(temp1) + pt%h(temp2))
    !TODO: PERIODIC CORRECTION
    xij = pt%x(temp1,:) - pt%x(temp2,:) 
    if ((xij(1)*xij(1) + xij(2)*xij(2) + xij(3)*xij(3)) < cellfac*cellfac*h*h) then
      CheckRadius = .true.
    end if
  end function CheckRadius

  subroutine AllocateNbPair(temp1, temp2, t)
    integer ,intent(in):: temp1,temp2, t
    integer:: i,j
    if (CheckRadius(temp1,temp2))  then
      if (nballoc_pass) then 
        pairs_t(t,pair_count(t),1) = temp1    
        pairs_t(t,pair_count(t),2) = temp2
      endif
      pair_count(t) = pair_count(t) +1
    end if
  end subroutine AllocateNbPair

  subroutine MainNeighbourSearch()
  use omp_lib 
  PRINT *, "Hello from process: ", OMP_GET_THREAD_NUM()
  end subroutine MainNeighbourSearch

  subroutine ClearNbData()
    pair_count (:) = 0
    pairs_t(:,:,:) = 0
    nbcount_t(:) = 0
    call CellReset()
  end subroutine ClearNbData

end module Neighbor