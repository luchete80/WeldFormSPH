program WeldFormSPH
use Integrate
use ParticleData
use Domain 

implicit none


!Type(Particle), POINTER :: pt

  real, dimension(1:3) :: V
  real :: dx, r, Rxy, Lz, h 
  
  dx    = 0.005
  Rxy  = 0.15
  Lz = 0.56
  r = dx / 2.0
  h = dx * 1.2

 V(1) = 0.;V(2) = 0.;V(3) = 0.
 !AddBoxLength(tag, V, Lx, Ly, Lz, r, Density,  h)			
 call AddBoxLength(0, V, 1.0,1.0,1.0,r, 2700.0, h)
 
 !call AddCylinderLength(0, V, Rxy, Lz, r)

end program WeldFormSPH

