program WeldFormSPH
use Integrate
use ParticleData
use Domain 
use Neighbor 

implicit none


!Type(Particle), POINTER :: pt

  real, dimension(1:3) :: V
  real :: dx, r, Rxy, Lz, h 
  integer:: i
  
  dx    = 0.005
  Rxy  = 0.15
  Lz = 0.56
  r = dx / 2.0
  h = dx * 1.2

 V(1) = 0.;V(2) = 0.;V(3) = 0.
 !AddBoxLength(tag, V, Lx, Ly, Lz, r, Density,  h)			
 call AddBoxLength(0, V, 1.0,1.0,1.0,r, 2700.0, h)
 
 do i = 1, part_count
 !write (*,*) "Particle", i ," position is ", pt%x(i,1), pt%x(i,1), pt%x(i,3)
 end do 
 !call AddCylinderLength(0, V, Rxy, Lz, r)
 
 call MainNeighbourSearch()

end program WeldFormSPH

