module Domain
use ParticleData
use ModPrecision, only : fp_kind
!use Neighbor 

implicit none 

integer :: Dim, part_count, Nproc
type(Particle)::pt
real(fp_kind), dimension(3):: dommax, dommin

contains 
  subroutine DomInit(proc)
    integer, intent(in)::proc
    nproc = proc
    Dim = 3
    DomMax (:) = -100000000000.0;
    DomMin (:) = 100000000000.0;
  end subroutine DomInit
  
  
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
    allocate (pt%x(part_count,3))
    allocate (pt%rho(part_count))
    allocate (pt%h(part_count))
    allocate (pt%m(part_count))
    !! THERMAL
    allocate (pt%cp_t(part_count))
    allocate (pt%t(part_count))
    allocate (pt%k_t(part_count))
    
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
    pt%rho(:) = Density
    print *, "Particle mass ", pt%m(2)
    
    !MECHANICAL PROPERTIES
    !if (solver_type==1) then 
    allocate (pt%v(part_count,3))
    allocate (pt%a(part_count,3))
    allocate (pt%sigma(part_count,3,3))
    !end if
    
  end subroutine AddBoxLength

!-------------------------------------------------------------------
!-  purpose:                                                       -
!-     calculate particle count in certain row                     -
!-     Return half (on the quadrant) particle count from a single position in an axis
!-     ri: particle separation 
!-     R: radius
!-------------------------------------------------------------------
  integer function CalcHalfPartCount(ri, R, xinc)
    implicit none 
    real(fp_kind), intent(in):: ri, R
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
      do while( rad <= R -r )
        yp = yp + 2.*r;
        rad = sqrt(yp*yp + xp*xp);
        ypartcount = ypartcount + 1;
      end do 
      ypartcount = ypartcount - 1;
    end if
    
    CalcHalfPartCount =  ypartcount
  end function CalcHalfPartCount
  
  subroutine AddCylinderLength(tag, V, Rxy, Lz, r)
    implicit none 
    integer, intent(in):: tag
    !real(fp_kind), intent(in), allocatable :: V
    real(fp_kind), dimension(1:), intent(in)  :: V ! input
    real(fp_kind), intent(in):: r, Lz, Rxy
    real(fp_kind), dimension (1:3) :: Xp
    real(fp_kind) :: numpartxy,numxpart,numypart, yinc_sign
    !Function vars definitions
    integer :: i,j,k, part_per_row, xinc, yinc

    Xp(3) = V(3) + r
    numpartxy = calcHalfPartCount(r, Rxy, 1)
    write(*,*) "Particles at plane xy: ", numpartxy
    
    part_per_row = 0
    write(*,*) "Vector value is ", r
    !Calculate row count for non ghost particles
    k=0
    do while (Xp(3) <= (V(3)+Lz -r) )
      k = k + 1 
      Xp(3) = Xp(3) + 2.* r
      numypart = 2*numpartxy
      do j=1,numypart
				numxpart = calcHalfPartCount(r, Rxy, yinc)
				do i=1,numxpart
					if (Xp(3) == V(2)+r) then
						part_per_row = part_per_row  + 1
            
          end if
        end do !i
				Xp(1) = V(1) - r - (2.*r*(numxpart - 1) ) !First increment is radius, following ones are 2r
      end do !j
      part_per_row = part_per_row + 1      
      !Xp(2) = V(2) - r - (2.*r*(numpartxy - 1) ); //First increment is radius, following ones are 2r     
      write(*,*) "xp3", Xp(3)      
    end do
    write(*,*) "Cylinder row count is ", k
    
    write(*,*) "Particle count is ", part_per_row
    
  end subroutine AddCylinderLength

end module Domain