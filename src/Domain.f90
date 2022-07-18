module Domain
use ParticleData
!use Neighbor 

implicit none 

integer :: Dim, part_count
type(Particle)::pt

contains 
  subroutine Init()
  
  end subroutine Init
  
  
  subroutine AddBoxLength(tag, V, Lx, Ly, Lz, r, Density,  h)			
    integer, intent(in):: tag
    !real, intent(in), allocatable :: V
    real, dimension(1:3), intent(in)  :: V ! input
    real, intent(in):: r, Lx, Ly, Lz, Density, h  
    real, dimension (1:3) :: Xp
    integer :: i, p
    Dim = 3
    do i=1, Dim
      Xp(i) = V(i) + r
    end do 
    
    part_count = 1
    do while (Xp(3) <= (V(3)+Lz-r))
      do while (Xp(2) <= (V(2)+Ly-r))
        do while (Xp(1) <= (V(1)+Lx-r))
          part_count = part_count +1
          Xp(1) = Xp(1) + 2 * r
        end do
        Xp(2) = Xp(2) + 2 * r
      end do
      Xp(3) = Xp(3) + 2 * r
    end do
    allocate (pt%x(part_count+10,3))

    do i=1, Dim
      Xp(i) = V(i) + r
    end do    
    write (*,*) "xp ", V(1)    
    !write(*,*) "Box Particle Count is ", part_count
    p = 1
    do while (Xp(3) <= (V(3)+Lz-r))
      do while (Xp(2) <= (V(2)+Ly-r))
        do while (Xp(1) <= (V(1)+Lx-r))
          do i=1, Dim
            pt%x(p,i) = Xp(i)
            !write (*,*) "particle " ,p, " xp ", Xp(i)
          end do 
          p = p + 1
          Xp(1) = Xp(1) + 2 * r
        end do
        Xp(2) = Xp(2) + 2 * r
      end do 
      Xp(3) = Xp(3) + 2 * r
    end do
    
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