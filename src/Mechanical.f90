module Mechanical

!Include mat_const.inc
real, allocatable, dimension(:) :: temp_pair
contains
  real function CalcAccIncNb(i, j)
    use omp_lib
    use Domain
    use Neighbor
    use Kernels
    implicit none

    real :: K, GK, h, xij(3), vij(3), cij, muij, rij
    integer,intent(in) :: i, j 

    xij(:) = pt%x(i,:) - pt%x(j,:)
    vij(:) = pt%v(i,:) - pt%v(j,:)
    h = 0.5 * (pt%h(i) + pt%h(j))
    
    rij = norm2(xij)
    GK = GradKernel (rij/h,h)
    K = GradKernel (rij/h,h)
    cij = 0.5*(pt%cs(i) + pt%cs(j))
    
    muij = h * dot_product(xij,vij) / (rij*rij+0.01*h*h)
    if (dot_product(xij,vij)<0.0) then
      
    end if
    !Assuming Gradient Type 0
    !Mult( GK*xij , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , temp);
    
  end function CalcAccIncNb

end module Mechanical