module Mechanical
use ModPrecision, only : fp_kind
!Include mat_const.inc
real(fp_kind), allocatable, dimension(:,:) :: force_pair
contains

  ! function outer_product(A,B) result(AB)
    ! double precision, intent(in) :: A(:),B(:)
    ! double precision, allocatable :: AB(:,:)
    ! integer :: nA,nB
  ! !  nA=size(A)
  ! !  nB=size(B)
    ! allocate(AB(3,3))
    ! AB = spread(source = A, dim = 2, ncopies = 3) * &
         ! spread(source = B, dim = 1, ncopies = 3)
  ! end function outer_product
  ! !!!!!!------------------------------------
  ! !CHECK WICH TO USE, THIS IS NOT ALLOCATING
  ! !!----------------------------------------
  ! subroutine outer_product_s(A,B, AB) 
    ! double precision, intent(in) :: A(:),B(:)
    ! double precision, allocatable :: AB(:,:)
    ! integer :: nA,nB
    ! AB = spread(source = A, dim = 2, ncopies = 3) * &
         ! spread(source = B, dim = 1, ncopies = 3)
  ! end subroutine outer_product_s

  real(fp_kind) function CalcAccIncNb(i, j)
    use omp_lib
    use Domain
    use Neighbor
    use Kernels
    implicit none

    real(fp_kind) :: K, GK, h, xij(3), vij(3), cij, muij, rij, temp(3,1), PIij(3,3), xijm(3,1)
    integer,intent(in) :: i, j 

    xij(:) = pt%x(i,:) - pt%x(j,:)
    xijm(:,1) = pt%x(i,:) - pt%x(j,:)
    vij(:) = pt%v(i,:) - pt%v(j,:)
    h = 0.5 * (pt%h(i) + pt%h(j))
    
    rij = norm2(xij)
    !GK = GradKernel (rij/h,h)
    !K = Kernel (rij/h,h)
    cij = 0.5*(pt%cs(i) + pt%cs(j))
    
    muij = h * dot_product(xij,vij) / (rij*rij+0.01*h*h)
    if (dot_product(xij,vij)<0.0) then
      
    end if
    !call outer_product_s(GK*xij,PIij, temp)
    !Assuming Gradient Type 0
    !Mult( GK*xij , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , temp);
    temp = matmul(PIij, GK*xijm)
    !force_pair(
    
  end function CalcAccIncNb

end module Mechanical