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
   !schedule (static)
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
  
  !Paralelization is done by table
  subroutine CalcTempInc()
    use omp_lib
    use Domain
    use Neighbor
    use Kernels
    implicit none
    
    real(fp_kind) :: GK, h, xij(3), nijinv,gkij(3)
    integer :: i, j, p, pp
    !!$omp parallel do num_threads(nproc) private (p,pp, i,j,h)
    !$omp parallel do num_threads(nproc) private (p,pp, i,j,h) schedule (static)
    do p = 1, nproc
      do pp = 1, pair_count(p)
        i = pairs_t(p,pp,1)
        j = pairs_t(p,pp,2)
        xij(:) = pt%x(i,:) - pt%x(j,:)
        h = 0.5 * (pt%h(i) + pt%h(j))
        GK = GradKernel (norm2 (xij)/h,h)
        !print *, "pair ", pp
        temp_pair(first_pair_perproc(p) + pp ) = pt%m(j)/pt%rho(j) * 4. * ( pt%k_t(i) * pt%k_t(j)) / (pt%k_t(i) + pt%k_t(j)) & 
           * ( pt%t(i) - pt%t(j)) * dot_product( xij , GK * xij )/ (dot_product(xij,xij)) 
      end do ! pairs
    end do !procs
    !$omp end parallel do    
  end subroutine CalcTempInc
  
  subroutine TempReduction (dTdt)
    !#pragma omp parallel for schedule (static) num_threads(Nproc)
    !for (int i=0; i<solid_part_count;i++)
    !  Particles[i]->a = 0.;
    !#pragma omp parallel for schedule (static) num_threads(Nproc)
    !Static: Split the loop variable evenly between threads beforehand (at compile-time);
    use Neighbor!, only : ipair_t
    use ParticleData
    use Domain 
    implicit none
    real(fp_kind), intent(out)::dTdt(part_count)
    integer :: i,n
    
    !$omp parallel do num_threads(nproc) private (i,n) schedule (static) 
    do i = 1, part_count
      do n = 1, ipair_t(i)  
        dTdt(i) = dTdt(i) + pt%rho(Anei_t(i,n)) * temp_pair(Anei_t(i,n))
      end do
      do n = 1, jpair_t(i)  
        dTdt(i) = dTdt(i) + pt%rho(Anei_t(i,maxnbcount - n)) * temp_pair(Anei_t(i,maxnbcount - n))
      end do
    end do
  end subroutine TempReduction

end module Thermal