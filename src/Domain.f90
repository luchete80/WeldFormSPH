module Domain
use ParticleData
!use Neighbor 

implicit none 

integer :: Dim, part_count, Nproc
type(Particle)::pt
real, dimension(3):: dommax, dommin

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
    !real, intent(in), allocatable :: V
    real, dimension(1:3), intent(in)  :: V ! input
    real, intent(in):: r, Lx, Ly, Lz, Density, h  
    real, dimension (1:3) :: Xp
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
  
    pt%m(:) = Density * Lx * Ly * Lz / part_count
    print *, "Particle mass ", pt%m(2)
    
  end subroutine AddBoxLength

!-------------------------------------------------------------------
!-  purpose:                                                       -
!-     calculate particle count in certain row                     -
!-     Return half (on the quadrant) particle count from a single position in an axis
!-     ri: particle separation 
!-     R: radius
!-------------------------------------------------------------------
  ! integer function CalcHalfPartCount(ri, R, xinc)
    ! implicit none 
    ! real, intent(in):: ri, R
    ! integer intent (in):: xinc
    ! !decls
    ! real :: xp, yp, rad
    ! integer:: ypartcount(-1)
    ! ! if ( xinc > 0 )
      ! ! ypartcount = 1
      ! ! yp = r
      ! ! xp = r + (xinc - 1 ) *2.*r; 
      ! ! rad = sqrt(yp*yp + xp*xp);
      ! ! do while( rad <= R -r )
        ! ! yp += 2.*r;
        ! ! rad = sqrt(yp*yp + xp*xp);
        ! ! ypartcount++;
      ! ! end 
      ! ! ypartcount-=1;
    
    ! CalcHalfPartCount =  ypartcount
  ! end function CalcHalfPartCount
  
  subroutine AddCylinderLength(tag, V, Rxy, Lz, r)
    implicit none 
    integer, intent(in):: tag
    !real, intent(in), allocatable :: V
    real, dimension(1:), intent(in)  :: V ! input
    real, intent(in):: r, Lz, Rxy
    !Function vars definitions
    integer :: k
    real :: zp 
    zp = V(2) + r
    
    write(*,*) "Vector value is ", r
    !Calculate row count for non ghost particles
    k=0
    do while (zp <= (V(3)+Lz -r) )
      k = k + 1 
      zp = zp + 2.* r      
    end do
    write(*,*) "Cylinder row count is ", k
    
  end subroutine AddCylinderLength

end module Domain