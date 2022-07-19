module Neighbor
!use Domain

implicit none 

! Neighbor data
integer, dimension(3):: cellno
integer:: maxnbcount
real, dimension(3):: blpf, trpr, cellsize
real :: rhomax, cellfac, hmax
integer, allocatable, dimension(:,:,:) :: HOC, pairs_t! Head of chain, pairs (proc, first, second)
integer, allocatable, dimension(:,:) :: cc !cell 
integer, allocatable, dimension(:) :: LL, pair_count, nbcount_t    !Linked list, paircount (per thread)
! REDUCTION ARRAYS (Nishimura)
integer, allocatable, dimension(:,:) :: Anei_t, Aref_t
integer, allocatable, dimension(:) :: ipair_t, jpair_t, first_pair_perproc ! ipair(part), jpair (part), first_pair_perproc(t)
!-----

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
    nballoc_pass = .false.
    
    call ClearNbData()
    
    maxnbcount = 80 !PER PARTICLE
  end subroutine InitNb

  subroutine AllocateNbData()
    use Domain, only : nproc
    implicit none
    integer:: i, tot
    
    call CellInitiate()
    call ListGenerate()
    call MainNeighbourSearch();
    tot = 0
    do i = 1, nproc
      if (pair_count(i)>tot) then
      tot = pair_count(i)
      !tot = tot + pair_count(i)
      end if 
    end do
    tot = int (real(tot) * 1.25)
    print *, "Allocated ", tot, "pairs "
    allocate (pairs_t(nproc,tot,2))
    deallocate(HOC)
    call ClearNbData()
    nballoc_pass = .true.
  end subroutine AllocateNbData
  
  subroutine CellInitiate ()
    use Domain
    implicit none
    !integer, intent(in)::part_count
    !integer, dimension(3), intent(in)::cellsize
    integer:: i,p
    blpf(:) = pt%x(1,:)
    trpr(:) = pt%x(1,:)
    
    write (*,*) "Looking for particles ", part_count
    do p=1, part_count
      do i=1,Dim
        !write (*,*) "i ",i, " x(p,i) ",pt%x(p,i), "trpr ", trpr(i)
        if ( pt%x(p,i) > trpr(i) ) then 
          trpr(i) = pt%x(p,i)
        end if
        if ( pt%x(p,i) < blpf(i) ) then 
          !print *, "FOUND , part ", p, " dim ", i
          blpf(i) = pt%x(p,i)
        end if
        if ( pt%rho(p) > rhomax ) then 
          rhomax = pt%rho(p)
        end if
      end do
    end do

    ! Override the calculated domain size

	!!Because of Hexagonal close packing in x direction domain is modified
	!if (!BC.Periodic(0)) {TRPR(0) += hmax/2;	BLPF(0) -= hmax/2;}else{TRPR(0) += R; BLPF(0) -= R;}
    trpr(:) = trpr(:) + hmax/2.0; blpf(:) = blpf(:) - hmax/2.0
    write (*,*) "trpr ", trpr(:)
    write (*,*) "blpf ", blpf(:)
    
    write (*,*) "Done, allocating cellno "
    !Calculate Cell properties
    do i = 1, Dim
      if (real(ceiling(((trpr(i)-blpf(i))/(cellfac*hmax)))-((trpr(i)-blpf(i))/(cellfac*hmax))) < hmax/10.) then
        cellno(i) = ceiling((trpr(i)-blpf(i))/(cellfac*hmax))
      else 
        cellno(i) = floor(((trpr(i)-blpf(i))/(cellfac*hmax)))
      end if
    end do
    
    print *, "test ", real(CellNo(:))
    cellsize(:) = (trpr(:)-BLPF(:))/real(CellNo(:))

    write (*, *) "Cell No ", cellno(:)
    write (*, *) "Cell Size ", cellsize(:)
    !Periodic BC correction

    !Initiate Head of Chain array for Linked-List
    allocate (HOC(cellno(1),cellno(2),cellno(3)))
    HOC(:,:,:) = -1
   
  end subroutine CellInitiate 

  subroutine ListGenerate()
    use Domain
    implicit none
    integer, dimension(3)::ijk
    integer ::a, d, temp
    
    do a = 1, part_count 
      !write (*,*) "Part ",a
      !print *, "pt%x(a,d) ", pt%x(a,:), "blpf ", blpf, " cellsize", cellsize
      
      do d = 1, Dim
        !print *, "test ",real(floor(pt%x(a,d)-blpf(d)))
        !write (*,*) "blpf(d) ",blpf(d), "cellsize ", cellsize(d)        
        ijk(d) = floor((pt%x(a,d)-blpf(d))/cellsize(d))
        !print *, "particle ", a, " ijk ", ijk(:)
        if (ijk(d) < 0) then
          ijk(d) = 0
        else if (ijk(d) >= cellno(d)) then !Original  if (j>=CellNo(1))j=CellNo(1)-1;
          ijk(d)=CellNo(d)-1;
        end if  
      end do
      
      temp = HOC(ijk(1)+1,ijk(2)+1,ijk(3)+1)
      HOC(ijk(1)+1,ijk(2)+1,ijk(3)+1) = a;
      ll(a) = temp
      cc(a,:) = ijk(:) + 1
    end do 
  print *, "List generation done, HOC ", HOC(:,:,:) 

  end subroutine ListGenerate

  subroutine CellReset()
    use Domain, only : part_count
    integer::a
    !do a=1, part_count
    ll(:) = -1
    cc(:,:) = -1
    !end do
    !HOC(:,:,:) = -1 
    !deallocate(HOC)
  end subroutine CellReset

  subroutine YZPlaneCellsNeighbourSearch(q1)
    use omp_lib
    implicit none
    integer, intent(in)  :: q1
    integer :: q2, q3, temp1, temp2, t, i, j
    t = omp_get_thread_num() + 1
    print *,"******- thread ", t
    do q3 = 1, CellNo(3)
      do q2 = 1, CellNo(2)
        if (HOC(q1,q2,q3)== -1) then
        else
          temp1 = HOC(q1,q2,q3)
          do while (temp1 .ne. -1)
            temp2 = LL(temp1)
            do while (temp2 .ne. -1)
              call AllocateNbPair(temp1,temp2,t)
              temp2 = LL(temp2)
            end do 
          
            !(q1 + 1, q2 , q3)
            if (q1+1 < CellNo(1)+1) then !!Original: if (q1+1 < CellNo(0)) 
              temp2 = HOC(q1+1,q2,q3)
              do while (temp2 .ne. -1) 
                call AllocateNbPair(temp1,temp2,t)
                temp2 = LL(temp2)
              end do !while temp2!=-1
            end if ! (q1 + 1, q2 , q3)   

            !(q1 + a, q2 + 1, q3) & a(-1,1)
            if (q2+1< CellNo(2)+1) then !(q2+1< CellNo(1))
              do i = q1 - 1, q1 + 1 !for (int i = q1-1; i <= q1+1; i++) 
                if (i<CellNo(1)+1 .and. i>=1) then !ORIG if (i<CellNo(0) .and. i>=0) 
                  temp2 = HOC(i,q2+1,q3)
                  do while (temp2 .ne. -1)
                    call AllocateNbPair(temp1,temp2,t)
                    temp2 = LL(temp2)
                  end do
                end if
              end do
            end if

            !(q1 + a, q2 + b, q3 + 1) & a,b(-1,1) => all 9 cells above the current cell
            if (q3+1< CellNo(3)+1) then
              do j = q2-1, q2 + 1 ! for (int j=q2-1; j<=q2+1; j++)
                do i = q1 - 1, q1 + 1 ! for (int i=q1-1; i<=q1+1; i++) {
                  if (i<CellNo(1)+1 .and. i>=1 .and. j<CellNo(2)+1 .and. j>=1) then ! if (i<CellNo(0) && i>=0 && j<CellNo(1) && j>=0) {
                    temp2 = HOC(i,j,q3+1)
                    do while (temp2 .ne. -1)
                      call AllocateNbPair(temp1,temp2,t)
                      temp2 = LL(temp2)
                    end do
                  end if
                end do
              end do !j
            end if
            temp1 = LL(temp1)        
          end do !temp1 .ne. -1 
        end if !HOC .ne. -1
      end do !q2
    end do !q3
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
  implicit none
  integer :: q1
  PRINT *, "Hello from process: ", OMP_GET_THREAD_NUM()
  !TODO: MAKE PERIODIC THINGS
  !$omp parallel do 
  !schedule (dynamic) num_threads(Nproc)
  !print *, "Cell no", CellNo
  do q1=1, CellNo(1)
    call YZPlaneCellsNeighbourSearch(q1)
  end do
  !$omp end parallel do
  print *, "Nb Search done."
  print *, "Pair count ", pair_count(:)
  end subroutine MainNeighbourSearch

  subroutine ClearNbData()
    pair_count (:) = 0
    !pairs_t(:,:,:) = 0
    nbcount_t(:) = 0
    call CellReset()
  end subroutine ClearNbData

  subroutine InitRedArraysOnce()
    use Domain, only : part_count, nproc
    allocate(Anei_t(part_count,maxnbcount))
    allocate(Aref_t(part_count,maxnbcount))
    allocate (ipair_t(part_count))
    allocate (jpair_t(part_count))
    allocate (first_pair_perproc(nproc))
  end subroutine InitRedArraysOnce
  
  subroutine CalcPairPosList()
    use omp_lib
    use Domain, only:nproc
    implicit none 
    integer :: p, i, j, pp, pair_tot_count
    first_pair_perproc(:) = 0
    pair_tot_count = 0
    do p=1, nproc
      if ( p > 0 ) then 
        do j = 1, p
          first_pair_perproc(p) = first_pair_perproc(p) + pair_count(p)          
        end do
      end if
      pair_tot_count = pair_tot_count + pair_count(p)   
    end do
    
    !!$omp parallel do 
    Anei_t(:,:) = 0;Aref_t(:,:) = 0
    ipair_t(:)=0; jpair_t(:)=0
    !!$omp end parallel do     
    
    do p = 1, nproc
      do pp = 1, pair_count(pp)
        print *, "proc, pair, i, j ", p, ", ", pp, ", ", i, ", " ,j
        i = min(pairs_t(p,pp,1),pairs_t(p,pp,2)) 
        j = max(pairs_t(p,pp,1),pairs_t(p,pp,2)) 


        Anei_t(i,ipair_t(i))                   = j  !!Only stores j>i
        Anei_t(j,maxnbcount - jpair_t(j))  = i  !!Only stores j>i
        
        ! Aref_t(i,ipair_t(i)) = first_pair_perproc(p)+pp !!
        ! Aref_t(j,maxnbcount - jpair_t(j)) = first_pair_perproc(p) + pp
        
        ipair_t(i) = ipair_t(i) + 1             !!ngji in 
        jpair_t(j) = jpair_t(j) + 1             !!njli, pairs in which j has particles with index smaller than it        
      end do ! pairs
    end do !procs
    
    print *, "Nb part 1"
    do i=1,maxnbcount
      print *, Anei_t(1,i)
    end do 
  end subroutine CalcPairPosList
  
end module Neighbor