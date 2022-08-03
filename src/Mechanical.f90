module Mechanical
use ModPrecision, only : fp_kind
!Include mat_const.inc
real(fp_kind), allocatable, dimension(:) :: temp_pair
contains

function outer_product(A,B) result(AB)
  double precision, intent(in) :: A(:),B(:)
  double precision, allocatable :: AB(:,:)
  integer :: nA,nB
!  nA=size(A)
!  nB=size(B)
  allocate(AB(3,3))
  AB = spread(source = A, dim = 2, ncopies = 3) * &
       spread(source = B, dim = 1, ncopies = 3)
end function outer_product

  real(fp_kind) function CalcAccIncNb(i, j)
    use omp_lib
    use Domain
    use Neighbor
    use Kernels
    implicit none

    real(fp_kind) :: K, GK, h, xij(3), vij(3), cij, muij, rij, temp(3,3), PIij(3,3)
    integer,intent(in) :: i, j 

    xij(:) = pt%x(i,:) - pt%x(j,:)
    vij(:) = pt%v(i,:) - pt%v(j,:)
    h = 0.5 * (pt%h(i) + pt%h(j))
    
    rij = norm2(xij)
    !GK = GradKernel (rij/h,h)
    !K = Kernel (rij/h,h)
    cij = 0.5*(pt%cs(i) + pt%cs(j))
    
    muij = h * dot_product(xij,vij) / (rij*rij+0.01*h*h)
    if (dot_product(xij,vij)<0.0) then
      
    end if
    !Assuming Gradient Type 0
    !Mult( GK*xij , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , temp);
    !temp = matmul(GK*xij, PIij)
    
  end function CalcAccIncNb

end module Mechanical