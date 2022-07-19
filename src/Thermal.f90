module Thermal

!Include mat_const.inc

contains
  real function CalcTempIncNb(i, j)
    use omp_lib
    use Domain
    use Neighbor
    use Kernels
    implicit none
    
    real :: GK, h, xij(3), nijinv,gkij(3)
    integer,intent(in) :: i, j 

    xij(:) = pt%x(i,:) - pt%x(j,:)
    h = 0.5 * (pt%h(i) + pt%h(j))
    GK = GradKernel (norm2 (xij)/h,h)
    CalcTempIncNb =  pt%m(j)/pt%rho(j) * 4. * ( pt%t(i) * pt%t(j)) / (pt%t(i) + pt%t(j)) * ( pt%t(i) - pt%t(j)) &
        * dot_product( xij , GK * xij )/ (norm2(xij)*norm2(xij)) 
  end function CalcTempIncNb
  
  !! TRADITIONAL FORM (SLOW)
  subroutine CalcTempIncPart(dTdt) !parallelized by particle
    use omp_lib
    use Domain
    use Neighbor
    use Kernels
    
    implicit none
    real, intent(out)::dTdt(part_count)
    real :: m, GK, h, xij(3)
    integer :: i, j, k, const
    !$omp parallel do num_threads(Nproc) 
    do i = 1, part_count
      const = 1.0/(pt%rho(i)*pt%cp_t(i))
      do k = 1, ipair_t(i)   
        j = Anei_t(i,k)
        dTdt(i) =  dTdt(i) + const * CalcTempIncNb(i, j)  
      end do

      do k = 1, jpair_t(i)   
        j = Anei_t(i,maxnbcount - jpair_t(k)+1)
        dTdt(i) = dTdt(i) + const * CalcTempIncNb(i, j)  
      end do
    end do
    !$omp end parallel do    
    !dTdt(i) = dTdt(i) + (1.0/Cp_t(i)) * (mass_t(j))
  end subroutine CalcTempIncPart

end module Thermal