module Neighbor
!use Domain

implicit none 

! Neighbor data
integer, dimension(3):: cellsize
real, dimension(3):: blpf, trpr
real :: rhomax, cellfac
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
        !rhomax = pt%rho(p,i)
      end if
    end do
  end do
  !Calculate Cell properties
  ! if (double (ceil(((TRPR(0)-BLPF(0))/(Cellfac*hmax)))-((TRPR(0)-BLPF(0))/(Cellfac*hmax)))<(hmax/10.0))
    ! CellNo[0] = int(ceil((TRPR(0)-BLPF(0))/(Cellfac*hmax)));
  ! else
    ! CellNo[0] = int(floor((TRPR(0)-BLPF(0))/(Cellfac*hmax)));

  ! if (double (ceil(((TRPR(1)-BLPF(1))/(Cellfac*hmax)))-((TRPR(1)-BLPF(1))/(Cellfac*hmax)))<(hmax/10.0))
    ! CellNo[1] = int(ceil((TRPR(1)-BLPF(1))/(Cellfac*hmax)));
  ! else
    ! CellNo[1] = int(floor((TRPR(1)-BLPF(1))/(Cellfac*hmax)));

  ! if (double (ceil(((TRPR(2)-BLPF(2))/(Cellfac*hmax)))-((TRPR(2)-BLPF(2))/(Cellfac*hmax)))<(hmax/10.0))
    ! CellNo[2] = int(ceil((TRPR(2)-BLPF(2))/(Cellfac*hmax)));
  ! else
    ! CellNo[2] = int(floor((TRPR(2)-BLPF(2))/(Cellfac*hmax)));

  ! CellSize  = Vec3_t ((TRPR(0)-BLPF(0))/CellNo[0],(TRPR(1)-BLPF(1))/CellNo[1],(TRPR(2)-BLPF(2))/CellNo[2]);  
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