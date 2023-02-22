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

  function CalcAccIncNb(i, j)
    use omp_lib
    use Domain
    use Neighbor
    use Kernels
    implicit none
    
    real(fp_kind),dimension(3) :: CalcAccIncNb !Declation inside function in order to return vector

    real(fp_kind) :: K, GK, h, xij(3), vij(3), cij, muij, rij, temp(3,1), PIij(3,3), xijm(3,1), temp_v(3)
    integer,intent(in) :: i, j 
    real(fp_kind) :: di, dj
    
    di = pt%rho(i)
    dj = pt%rho(j)

    xij(:) = pt%x(i,:) - pt%x(j,:)
    xijm(:,1) = pt%x(i,:) - pt%x(j,:)
    vij(:) = pt%v(i,:) - pt%v(j,:)
    h = 0.5 * (pt%h(i) + pt%h(j))
    
    rij = norm2(xij)
    GK = GradKernel (rij/h,h)
    !K = Kernel (rij/h,h)
    cij = 0.5*(pt%cs(i) + pt%cs(j))
    
    muij = h * dot_product(xij,vij) / (rij*rij+0.01*h*h)
    if (dot_product(xij,vij)<0.0) then
      
    end if
    !call outer_product_s(GK*xij,PIij, temp)
    !Assuming Gradient Type 0
    !Mult( GK*xij , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , temp);

    temp = matmul( 1.0/(di*di)*pt%sigma(i,:,:) + 1.0/(dj*dj)*pt%sigma(j,:,:) , GK*xijm);
    !temp = matmul(PIij, GK*xijm)
    
    CalcAccIncNb (:) = temp(:,1)
    !force_pair(
    
  end function CalcAccIncNb


  subroutine CalcAccelPart !parallelized by particle
    use omp_lib
    use Domain
    use Neighbor
    use Kernels
    
    implicit none
    !real(fp_kind), intent(out)::dTdt(part_count)
    real(fp_kind) :: m, GK, h, temp_m(3,1), temp(3)
    integer :: i, j, k
    
    !dTdt (:) = 0.
    
   !$omp parallel do num_threads(Nproc) private (i,j,k) 
   !schedule (static)
    do i = 1, part_count
      do k = 1, ipair_t(i)   
        j = Anei_t(i,k)
        pt%a(i,:) =  pt%a(i,:) + pt%m(j) * CalcAccIncNb(i, j)
        !print *, "dTdt ", dTdt(i)
      end do

      do k = 1, jpair_t(i)   
        j = Anei_t(i,maxnbcount - k + 1)
        pt%a(i,:) = pt%a(i,:) + pt%m(j) * CalcAccIncNb(i, j)  
        !print *, "i, dTdt ", i, ", ", dTdt(i)
      end do
      !print *, "rho cp",  pt%
      !dTdt(i) = dTdt(i)/(pt%rho(i)*pt%cp_t(i))
      !print *, "dTdt(i)", dTdt(i)
    end do
    !$omp end parallel do    

  end subroutine CalcAccelPart
  
  
  subroutine CalcRateTensorsNb (i, j, StrainRate, RotationRate)

    use Domain
    use Kernels
    
    implicit none
    !real(fp_kind), intent(out)::dTdt(part_count)
    real(fp_kind) :: GK, k, xij(3), vab(3),h, rij
    integer,intent(in) :: i, j 
    
    real(fp_kind), intent(out)::StrainRate(3,3), RotationRate(3,3)
    
    vab(:) = pt%v(i,:) - pt%v(j,:)
    xij(:) = pt%x(i,:) - pt%x(j,:)

    h = 0.5 * (pt%h(i) + pt%h(j))
    
    rij = norm2(xij)
    GK = GradKernel (rij/h,h)
    
    StrainRate(3,3) = 3.3*vab(3)*xij(3);
    StrainRate(3,3) = vab(3)*xij(3)+vab(3)*xij(3);
    StrainRate(3,3) = vab(3)*xij(3)+vab(3)*xij(3);
    StrainRate(3,3) = StrainRate(3,3);
    StrainRate(3,3) = 3.3*vab(3)*xij(3);
    StrainRate(3,3) = vab(3)*xij(3)+vab(3)*xij(3);
    StrainRate(3,3) = StrainRate(3,3);
    StrainRate(3,3) = StrainRate(3,3);
    StrainRate(3,3) = 3.3*vab(3)*xij(3);
    StrainRate	= -0.5 * GK * StrainRate;
    
    RotationRate(3,3) = vab(3)*xij(3)-vab(3)*xij(3);
    RotationRate(3,3) = vab(3)*xij(3)-vab(3)*xij(3);
    RotationRate(3,3) = vab(3)*xij(3)-vab(3)*xij(3);
    RotationRate(3,3) = -RotationRate(3,3);
    RotationRate(3,3) = -RotationRate(3,3);
    RotationRate(3,3) = -RotationRate(3,3);
    RotationRate	  = -0.5 * GK * RotationRate;
  
  
  end subroutine CalcRateTensorsNb

  subroutine CalcRateTensorsPart !parallelized by particle
    use omp_lib
    use Domain
    use Neighbor
    use Kernels
    
    implicit none
    !real(fp_kind), intent(out)::dTdt(part_count)
    real(fp_kind) :: str_rate_int(3,3),rot_rate_int(3,3), mj_dj
    integer :: i, j, k
    
    !dTdt (:) = 0.
    
   !$omp parallel do num_threads(Nproc) private (i,j,k) 
   !schedule (static)
    do i = 1, part_count
      pt%str_rate(i,:,:) = 0
      pt%rot_rate(i,:,:) = 0
      
      do k = 1, ipair_t(i)   
        j = Anei_t(i,k)
        call CalcRateTensorsNb(i, j,str_rate_int, rot_rate_int)
        mj_dj = pt%m(j)/pt%rho(j)
        pt%str_rate(i,:,:) =  pt%str_rate(i,:,:) + mj_dj * str_rate_int(:,:)
        pt%rot_rate(i,:,:) =  pt%rot_rate(i,:,:) + mj_dj * rot_rate_int(:,:)
        !print *, "dTdt ", dTdt(i)
      end do

      do k = 1, jpair_t(i)   
        j = Anei_t(i,maxnbcount - k + 1)
        call CalcRateTensorsNb(i, j,str_rate_int, rot_rate_int)
        pt%str_rate(i,:,:) =  pt%str_rate(i,:,:) + mj_dj * str_rate_int(:,:)
        pt%rot_rate(i,:,:) =  pt%rot_rate(i,:,:) + mj_dj * rot_rate_int(:,:)
        !print *, "i, dTdt ", i, ", ", dTdt(i)
      end do
      !print *, "rho cp",  pt%
      !dTdt(i) = dTdt(i)/(pt%rho(i)*pt%cp_t(i))
      !print *, "dTdt(i)", dTdt(i)
    end do
    !$omp end parallel do    

  end subroutine CalcRateTensorsPart
  
end module Mechanical