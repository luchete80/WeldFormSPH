module Neighbor
!use Domain

implicit none 

! Neighbor data
integer, dimension(3):: cellsize
real, dimension(3):: blpf, trpr
!integer, dimension (:,:) :: Neib_t !Filled after pairs

contains

subroutine CellInitiate ()
  use Domain
  implicit none
  !integer, intent(in)::part_count
  !integer, dimension(3), intent(in)::cellsize
  integer:: i,p
  do p=1, part_count
    do i=1,Dim
      if(pt%x(i)>trpr(i)) 
        trpr(i) = pt%x(i)
      end if
    end do
  end do

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

end module Neighbor