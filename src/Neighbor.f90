module Neighbor
implicit none 
contains
subroutine ListGenerate(part_count, cellsize)
  implicit none
  integer, intent(in)::part_count
  integer, dimension(3), intent(in)::cellsize
  

  integer ::a
  do a = 1, part_count 
  end do 
end subroutine ListGenerate

end module Neighbor