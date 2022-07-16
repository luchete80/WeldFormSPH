program WeldFormSPH
use Integrate
use ParticleData
use Geometry 

implicit none

Type(Particle):: pt
!Type(Particle), POINTER :: pt

 allocate (pt%x(100))
 !pt => x(10) = 0.
 pt%x(1) = 0. 

end program WeldFormSPH

