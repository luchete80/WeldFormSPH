module Neighbor
!use Domain

implicit none 

! Neighbor data
integer, dimension(3):: cellsize, cellno
real, dimension(3):: blpf, trpr
real :: rhomax, cellfac, hmax
integer, allocatable, dimension(:,:,:) :: HOC ! Head of chain
!integer, dimension (:,:) :: Neib_t !Filled after pairs

contains

subroutine CellInitiate ()
  use Domain
  implicit none
  !integer, intent(in)::part_count
  !integer, dimension(3), intent(in)::cellsize
  integer:: i,p
  blpf(:) = pt%x(1,:)
  
  do p=1, part_count
    do i=1,Dim
      if ( pt%x(p,i) > trpr(i) ) then 
        trpr(i) = pt%x(p,i)
      end if
      if ( pt%x(p,i) < blpf(i) ) then 
        trpr(i) = pt%x(p,i)
      end if
      if ( pt%rho(p) > rhomax ) then 
        rhomax = pt%rho(p)
      end if
    end do
  end do
  !Calculate Cell properties
  do i = 1, Dim
    if (real(ceiling((trpr(i)-blpf(i))/(cellfac*hmax)-(trpr(i)-blpf(i))/(cellfac*hmax)))<hmax/10.) then
      cellno(i) = (trpr(i)-blpf(i))/(cellfac*hmax)
    else 
      cellno(i) = floor((trpr(i)-blpf(i))/(cellfac*hmax))
    end if
  end do
  
  do i = 1, Dim
    cellno(i) = (TRPR(i)-BLPF(i))/CellNo(i)
  end do
  
  !Periodic BC correction

  !Initiate Head of Chain array for Linked-List
  allocate (HOC(cellno(1),cellno(2),cellno(3)))
  HOC(:,:,:) = -1
 
  !HOC = new int**[(int) CellNo[0]];
  ! for(int i =0; i<CellNo[0]; i++){
     ! HOC[i] = new int*[CellNo[1]];
     ! for(int j =0; j<CellNo[1]; j++){
         ! HOC[i][j] = new int[CellNo[2]];
         ! for(int k = 0; k<CellNo[2];k++){
            ! HOC[i][j][k] = -1;
         ! }
     ! }
  ! }  
end subroutine CellInitiate 

subroutine ListGenerate()
  use Domain
  implicit none
!  integer, intent(in)::part_count
!  integer, dimension(3), intent(in)::cellsize
  

  integer ::a
  do a = 1, part_count 
  end do 
end subroutine ListGenerate

subroutine MainNeighbourSearch()
use omp_lib 
PRINT *, "Hello from process: ", OMP_GET_THREAD_NUM()
end subroutine MainNeighbourSearch

end module Neighbor