module Thermal

!Include mat_const.inc

contains
  subroutine CalcTempIncNb(i, k)
    use omp_lib
    use Domain
    use Neighbor
    use Kernels
    implicit none
    
    real :: m, GK, h, xij(3), nijinv
    integer,intent(in) :: i, k 
    integer :: j    
    j = Anei_t(i,k)
    xij(:) = pt%x(i,:) - pt%x(j,:)
    h = 0.5 * (pt%h(i) + pt%h(j))
    GK = GradKernel (norm2 (xij)/h,h)
    !m =  pt%m(j)/pt%rho(j) * 4. * ( pt%t(i) * pt%t(j)) / (P1->k_T + P2->k_T) * ( P1->T - P2->T) * dot( xij , GK*xij )/ (norm2(xij)*norm(xij));  
  end subroutine CalcTempIncNb
  
  !! TRADITIONAL FORM (SLOW)
  subroutine CalcTempIncPart(dTdt) !parallelized by particle
    use omp_lib
    use Domain
    use Neighbor
    use Kernels
    implicit none
    real, intent(out)::dTdt(part_count)
    real :: m, GK, h, xij(3)
    integer :: i, j
    !$omp parallel do num_threads(Nproc) 
    do i = 1, part_count
      do j = 1, ipair_t(i)   
        call CalcTempIncNb(i, j)  
      end do
    end do
    !$omp end parallel do    
    !dTdt(i) = dTdt(i) + (1.0/Cp_t(i)) * (mass_t(j))
  end subroutine CalcTempIncPart

end module Thermal