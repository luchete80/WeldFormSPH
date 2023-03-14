module Domain
use ParticleData
use ModPrecision, only : fp_kind
!use Neighbor 

implicit none 

integer :: Dim, part_count, Nproc
type(Particle)::pt
real(fp_kind), dimension(3):: dommax, dommin

real(fp_kind)::mat_G, mat_E!TODO: change to material
real(fp_kind):: time

contains 
  subroutine DomInit(proc)
    integer, intent(in)::proc
    nproc = proc
    Dim = 3
    DomMax (:) = -100000000000.0;
    DomMin (:) = 100000000000.0;
  end subroutine DomInit
  
  subroutine AllocateParticles(pt_count)
    integer, intent(in):: pt_count
    
    part_count = pt_count
    !!!GENERAL 
    allocate (pt%x(part_count,3))
    allocate (pt%rho(part_count))
    allocate (pt%drhodt(part_count))
    allocate (pt%h(part_count))
    allocate (pt%m(part_count))
    allocate (pt%id(part_count))
    
    !! THERMAL
    allocate (pt%cp_t(part_count))
    allocate (pt%t(part_count))
    allocate (pt%k_t(part_count))
    
    !MECHANICAL PROPERTIES
    !if (solver_type==1) then 
    
    allocate (pt%v(part_count,3))
    allocate (pt%a(part_count,3))
    allocate (pt%disp(part_count,3))
    
    allocate (pt%sigma(part_count,3,3))
    allocate (pt%str_rate(part_count,3,3))
    allocate (pt%rot_rate(part_count,3,3))
    allocate (pt%shear_stress(part_count,3,3))
    allocate (pt%strain(part_count))
    allocate (pt%pressure(part_count))
    allocate (pt%cs(part_count))
    
    allocate (pt%sigma_eq(part_count))
    
    allocate (pt%rho_0(part_count))
    
    !end if  
  end subroutine
  
  
  subroutine AddBoxLength(tag, V, Lx, Ly, Lz, r, Density,  h)			
    integer, intent(in):: tag
    !real(fp_kind), intent(in), allocatable :: V
    real(fp_kind), dimension(1:3), intent(in)  :: V ! input
    real(fp_kind), intent(in):: r, Lx, Ly, Lz, Density, h  
    real(fp_kind), dimension (1:3) :: Xp
    integer :: i, p

    Xp(3) = V(3) + r
    part_count = 0
    do while (Xp(3) <= (V(3)+Lz))
      Xp(2) = V(2) + r
      do while (Xp(2) <= (V(2)+Ly))
        Xp(1) = V(1) + r
        do while (Xp(1) <= (V(1)+Lx))
          part_count = part_count +1
          Xp(1) = Xp(1) + 2 * r
        end do
        Xp(2) = Xp(2) + 2 * r
      end do
      Xp(3) = Xp(3) + 2 * r
    end do
    
    call AllocateParticles(part_count)
    
    write (*,*) "Box particle count ", part_count
    
      Xp(:) = V(:) + r
    write (*,*) "xp ", Xp(:)    
    !write(*,*) "Box Particle Count is ", part_count
    p = 1
    do while (Xp(3) <= (V(3)+Lz))
      Xp(2) = V(2) + r
      do while (Xp(2) <= (V(2)+Ly))
        Xp(1) = V(1) + r
        do while (Xp(1) <= (V(1)+Lx))
            pt%x(p,:) = Xp(:)
          !print *,"particle ",p , "X: ",Xp(:)
          p = p + 1
          Xp(1) = Xp(1) + 2 * r

        end do
        Xp(2) = Xp(2) + 2 * r
      end do 
      Xp(3) = Xp(3) + 2 * r
    end do
  
    do p = 1, part_count
      pt%h(p) = h
    end do
  
    pt%m(:)   = Density * Lx * Ly * Lz / part_count
    pt%rho(:)   = Density
    pt%rho_0(:) = Density
    print *, "Particle mass ", pt%m(2)
    
    pt%id(:) = tag
    
    
  end subroutine AddBoxLength

!-------------------------------------------------------------------
!-  purpose:                                                       -
!-     calculate particle count in certain row                     -
!-     Return half (on the quadrant) particle count from a single position in an axis
!-     ri: particle separation 
!-     R: radius
!-------------------------------------------------------------------
  integer function CalcHalfPartCount(r, Rxy, xinc)
    implicit none 
    real(fp_kind), intent(in):: r, Rxy
    integer, intent (in):: xinc
    !decls
    real(fp_kind) :: xp, yp, rad
    integer:: ypartcount
    
    ypartcount = -1
    
    if ( xinc > 0 ) then
      ypartcount = 1
      yp = r
      xp = r + (xinc - 1 ) *2.*r; 
      rad = sqrt(yp*yp + xp*xp);
      do while( rad <= Rxy -r )
        yp = yp + 2.*r
        rad = sqrt(yp*yp + xp*xp);
        ypartcount = ypartcount + 1;
      end do 
      ypartcount = ypartcount - 1;
    end if
    
    CalcHalfPartCount =  ypartcount
  end function CalcHalfPartCount

  integer function AllocateCylXY (V, Rxy, Lz, r, mem_alloc) !if mem_alloc only calculate particle number
    implicit none 
    !real(fp_kind), intent(in), allocatable :: V
    real(fp_kind), dimension(1:), intent(in)  :: V ! input
    real(fp_kind), intent(in):: r, Lz, Rxy
    real(fp_kind), dimension (1:3) :: Xp
    real(fp_kind) :: numpartxy,numxpart,numypart, yinc_sign
    logical, intent(in):: mem_alloc
    !Function vars definitions
    integer :: i,j,k, part_per_row, xinc, yinc, id_part

    Xp(3) = V(3) + r
    numpartxy = calcHalfPartCount(r, Rxy, 1)
    write(*,*) "Particles at plane xy: ", numpartxy
    
    part_per_row = 0
    write(*,*) "Vector value is ", r
    !Calculate row count for non ghost particles
    k=0
    id_part = 1
    do while (Xp(3) <= (V(3)+Lz -r) )
      numypart = 2*numpartxy
      Xp(2) = V(2) - r - (2.*r*(numpartxy - 1) ) !First increment is radius, following ones are 2r
			yinc = numpartxy !particle row from the axis
			yinc_sign = -1
      !write(*,*) "y part ", numypart
      do j=1,numypart
        !write(*,*) "j ", j
				numxpart = calcHalfPartCount(r, Rxy, yinc)
        Xp(1) = V(1) - r - (2.*r*(numxpart - 1) ) !First increment is radius, following ones are 2r
				do i=1,2*numxpart
          !write(*,*) "i ", i
          if (mem_alloc .eqv. .TRUE.) then 
            pt%x(id_part,:) = Xp(:)
            id_part = id_part + 1
          end if
					if (Xp(3) == V(3)+r) then
						part_per_row = part_per_row  + 1            
          end if
          Xp(1) = Xp(1) + 2.0 * r
        end do !i
				yinc = yinc + yinc_sign
				Xp(2) = Xp(2) + 2.0 * r
        if (yinc<1) then 
					yinc = 1
					yinc_sign = 1
				end if
      end do !j 
      !Xp(2) = V(2) - r - (2.*r*(numpartxy - 1) ); //First increment is radius, following ones are 2r     
      !write(*,*) "xp3", Xp(3)      
      k = k + 1 
      Xp(3) = Xp(3) + 2.* r
    end do
    write(*,*) "Cylinder row count is ", k
    
    write(*,*) "Particle count is ", part_per_row * k  
    
    AllocateCylXY = 2.0 * part_per_row * k 
  end function AllocateCylXY
  
  subroutine AddCylinderLength(tag, V, Rxy, Lz, r, Density, h)
    implicit none 
    integer, intent(in):: tag
    !real(fp_kind), intent(in), allocatable :: V
    real(fp_kind), dimension(1:), intent(in)  :: V ! input
    real(fp_kind), intent(in):: r, Lz, Rxy, Density, h
    real(fp_kind), dimension (1:3) :: Xp
    real(fp_kind) :: numpartxy,numxpart,numypart, yinc_sign
    !Function vars definitions
    integer :: i,j,k, part_per_row, xinc, yinc

    part_count = AllocateCylXY(V, Rxy, Lz, r, .false.) !First allocate
    
    call AllocateParticles(part_count)
    part_count = AllocateCylXY(V, Rxy, Lz, r, .true.) !First allocate

    pt%m(:)   = 3.1415926 * Density * Rxy * Rxy * Lz / part_count
    pt%rho(:) = Density
    print *, "Particle mass ", pt%m(1)
    
    pt%h(:) = h
    pt%id(:) = tag
    
  end subroutine AddCylinderLength

end module Domain