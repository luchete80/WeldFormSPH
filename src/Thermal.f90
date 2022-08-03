module Thermal

use ModPrecision, only : fp_kind
!Include mat_const.inc
real(fp_kind), allocatable, dimension(:) :: temp_pair
contains
  real(fp_kind) function CalcTempIncNb(i, j)
    use omp_lib
    use Domain
    use Neighbor
    use Kernels
    implicit none
    
    real(fp_kind) :: GK, h, xij(3), nijinv,gkij(3)
    integer,intent(in) :: i, j 

    xij(:) = pt%x(i,:) - pt%x(j,:)
    h = 0.5 * (pt%h(i) + pt%h(j))
    GK = GradKernel (norm2 (xij)/h,h)
    !print *, "GK ", GK
    !print *, "norm xij " , norm2(xij)
    !print *, "rho ", pt%rho(i), ", ", pt%rho(j)
    CalcTempIncNb =  pt%m(j)/pt%rho(j) * 4. * ( pt%k_t(i) * pt%k_t(j)) / (pt%k_t(i) + pt%k_t(j)) * ( pt%t(i) - pt%t(j)) &
        * dot_product( xij , GK * xij )/ (dot_product(xij,xij)) 
  end function CalcTempIncNb
  
  !! TRADITIONAL FORM (SLOW)
  subroutine CalcTempIncPart(dTdt) !parallelized by particle
    use omp_lib
    use Domain
    use Neighbor
    use Kernels
    
    implicit none
    real(fp_kind), intent(out)::dTdt(part_count)
    real(fp_kind) :: m, GK, h, xij(3)
    integer :: i, j, k
    
    dTdt (:) = 0.
    
    !$omp parallel do num_threads(Nproc) private (i,j,k)
    do i = 1, part_count
      do k = 1, ipair_t(i)   
        j = Anei_t(i,k)
        dTdt(i) =  dTdt(i) + CalcTempIncNb(i, j)  
        !print *, "dTdt ", dTdt(i)
      end do

      do k = 1, jpair_t(i)   
        j = Anei_t(i,maxnbcount - k + 1)
        dTdt(i) = dTdt(i) + CalcTempIncNb(i, j)  
        !print *, "i, dTdt ", i, ", ", dTdt(i)
      end do
      !print *, "rho cp",  pt%
      dTdt(i) = dTdt(i)/(pt%rho(i)*pt%cp_t(i))
      !print *, "dTdt(i)", dTdt(i)
    end do
    !$omp end parallel do    

  end subroutine CalcTempIncPart

end module Thermal