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
  
  subroutine StartVars ()
    use Domain
    pt%a(:,:) = 0.
    !pt%rho(:,:)
  end subroutine StartVars
  
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
   ! print *, "sigma", pt%sigma(i,:,:)
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
        if (i ==51) then 
        !print *, "part 51 position ", pt%x(52,:)
        !print *, "j " , j, "CalcAccIncNb ij ",CalcAccIncNb(i, j)
        !print *, "sigma j: ", j, ", " , pt%sigma(j,:,:) 
        end if
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
    
    StrainRate(1,1) = 2.0*vab(1)*xij(1);
    StrainRate(1,2) = vab(1)*xij(2)+vab(2)*xij(1);
    StrainRate(1,3) = vab(1)*xij(3)+vab(3)*xij(1);
    StrainRate(2,1) = StrainRate(1,2);
    StrainRate(2,2) = 2.0*vab(2)*xij(2);
    StrainRate(2,3) = vab(2)*xij(3)+vab(3)*xij(2);
    StrainRate(3,1) = StrainRate(1,3);
    StrainRate(3,2) = StrainRate(2,3);
    StrainRate(3,3) = 2.0*vab(3)*xij(3);
    StrainRate	= -0.5 * GK * StrainRate;
    
    RotationRate(1,1) = 0.0;     RotationRate(2,2) = 0.0;     RotationRate(3,3) = 0.0;
    RotationRate(1,2) = vab(1)*xij(2)-vab(2)*xij(1);
    RotationRate(1,3) = vab(1)*xij(3)-vab(3)*xij(1);
    RotationRate(2,3) = vab(2)*xij(3)-vab(3)*xij(2);
    RotationRate(2,1) = -RotationRate(1,2);
    RotationRate(3,1) = -RotationRate(1,3);
    RotationRate(3,2) = -RotationRate(2,3);
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
    real(fp_kind) :: GK, xij(3), vab(3),h, rij
    !dTdt (:) = 0.
    
   !$omp parallel do num_threads(Nproc) private (i,j,k) 
    do i = 1, part_count
      pt%str_rate(i,:,:) = 0.0
      pt%rot_rate(i,:,:) = 0.0
      
      do k = 1, ipair_t(i)   
        j = Anei_t(i,k)
        call CalcRateTensorsNb(i, j,str_rate_int, rot_rate_int)
        mj_dj = pt%m(j)/pt%rho(j)
        pt%str_rate(i,:,:) =  pt%str_rate(i,:,:) + mj_dj * str_rate_int(:,:)
        pt%rot_rate(i,:,:) =  pt%rot_rate(i,:,:) + mj_dj * rot_rate_int(:,:)
        if (i==52) then
            xij(:) = pt%x(i,:) - pt%x(j,:)

            h = 0.5 * (pt%h(i) + pt%h(j))
            
            rij = norm2(xij)
            GK = GradKernel (rij/h,h)
        print *, "nb k ", k, ", GK ", GK, "j particle ", j, "rot_rate inc",  rot_rate_int(:,:)
        end if
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
  
  real(fp_kind) function CalcDensIncNb(i, j)
    use omp_lib
    use Domain
    use Neighbor
    use Kernels
    
    implicit none
    real(fp_kind) :: GK, k, xij(3), h, rij, vij(3)
    integer,intent(in) :: i, j 
    
    xij(:) = pt%x(i,:) - pt%x(j,:)
    vij(:) = pt%v(i,:) - pt%v(j,:)
    h = 0.5 * (pt%h(i) + pt%h(j))
    
    rij = norm2(xij) 
    GK = GradKernel (rij/h,h)
    ! if (i==52 .and. j==675) then
    ! print *, "GK ", GK
    ! print *, "rij ", rij
    ! end if 
    CalcDensIncNb = dot_product(vij,GK*xij)
    
  end function CalcDensIncNb
  
  subroutine CalcDensIncPart
    use omp_lib
    use Domain
    use Neighbor
    use Kernels
    
    implicit none
    real(fp_kind) :: m
    integer :: i, j, k
    
    !dTdt (:) = 0.
    
   write (*,*) "ipair_t(i)   ", ipair_t(i)   
   print *, "Time ", time
   !$omp parallel do num_threads(Nproc) private (i,j,k) 
   !schedule (static)
    do i = 1, part_count
      do k = 1, ipair_t(i)   
        j = Anei_t(i,k)
        if (i==52) then
       ! print *," j ", j, "densij ",CalcDensIncNb(i,j)
        !print *, "vab ", pt%v(i,:) - pt%v(j,:)
        end if
        pt%rho(i) =  pt%rho(i) + pt%m(j)/pt%rho(j) * CalcDensIncNb(i, j)
        !print *, "dTdt ", dTdt(i)
      end do

      do k = 1, jpair_t(i)   
        j = Anei_t(i,maxnbcount - k + 1)
        ! if (i==1) then
        ! print *," j ", j
        ! end if
        pt%rho(i) = pt%rho(i) + pt%m(j)/pt%rho(j)  * CalcDensIncNb(i, j)  
        !print *, "i, dTdt ", i, ", ", dTdt(i)
      end do
      !print *, "rho cp",  pt%
      !dTdt(i) = dTdt(i)/(pt%rho(i)*pt%cp_t(i))
      !print *, "dTdt(i)", dTdt(i)
    end do
    !$omp end parallel do      
  end subroutine CalcDensIncPart
  
  
  
  subroutine CalcStressStrain (dt) 
  use Domain
  use Functions
	!Pressure = EOS(PresEq, Cs, P0,Density, RefDensity)

	!!!!Jaumann rate terms
	!Mat3_t RotationRateT, Stress,SRT,RS;
	!Trans(RotationRate,RotationRateT);
	!Mult(ShearStress,RotationRateT,SRT);
	!Mult(RotationRate,ShearStress,RS);
  implicit none
  real(fp_kind) :: RotationRateT(3,3), Stress(3,3), SRT(3,3), RS(3,3), ident(3,3)
  integer :: i
  real(fp_kind) ,intent(in):: dt
  
  real(fp_kind) :: p00
  
  p00 = 0.
  
  ident = 0.
  ident (1,1) = 1; ident (2,2) = 1.0; ident (3,3) = .1
  
  !$omp parallel do num_threads(Nproc) private (RotationRateT, Stress, SRT, RS)
  do i = 1, part_count
    !pt%pressure(i) = EOS(PresEq, Cs, P0,Density, RefDensity)
    pt%pressure(i) = EOS(0, pt%cs(i), p00,pt%rho(i), pt%rho_0(i))
    if (i==52) then
    !print *, "pt%pressure(i)", pt%pressure(i),", cs ", pt%cs(i), "p00", p00, ", rho", p00,pt%rho(i), ", rho 0", p00,pt%rho_0(i)
    end if
    RotationRateT = transpose (pt%rot_rate(i,:,:))
    SRT = MatMul(pt%shear_stress(i,:,:),RotationRateT)
    RS  = MatMul(pt%rot_rate(i,:,:), pt%shear_stress(i,:,:))
    
    !print *, "RS", RS
    pt%shear_stress(i,:,:)	= dt * (2.0 * mat_G *(pt%str_rate(i,:,:)-1.0/3.0 * &
                                 (pt%str_rate(i,1,1)+pt%str_rate(i,2,2)+pt%str_rate(i,3,3))*ident) &
                                 +SRT+RS) + pt%shear_stress(i,:,:)
    pt%sigma(i,:,:)			= -pt%pressure(i) * ident + pt%shear_stress(i,:,:)	!Fraser, eq 3.32
    !print *, "particle ", i, ", rot_rate ", pt%rot_rate(i,:,:)
    !pt%strain(i)			= dt*pt%str_rate(i + Strain;
  end do
  !$omp end parallel do    
	  !print *, "str_rate ", pt%str_rate(1,:,:)
	  !print *, "shear_stress ", pt%shear_stress(1,:,:)
    print *, "pressure ",pt%pressure(52)
	  print *, "sigma ", pt%sigma(52,:,:)
  !double dep = 0.;
  !double sig_trial = 0.;

	!ShearStress	= dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStress;

  !! FAIL : TODO

	!Sigma			= -Pressure * OrthoSys::I + ShearStress;	//Fraser, eq 3.32
	
  
	!!!!!if ( dep > 0.0 ) 
  !Strain	= dt*StrainRate + Strain;

  end subroutine CalcStressStrain

end module Mechanical