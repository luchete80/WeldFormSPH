program MechTest

use Integrate
use ParticleData
use Domain 
use Neighbor 
use Thermal
use omp_lib
use Thermal
use Mechanical
use ModPrecision, only : fp_kind

implicit none

  real(fp_kind), dimension(1:3) :: V
  real(fp_kind):: r, Lz, Rxy, dx, h 
  
  Lz = 0.56
  Rxy = 0.15
  
  dx = 0.015
  h	= dx*1.2 
  r = dx/2.
  
  call AddCylinderLength(0, V, Rxy, Lz + 2.0 * Lz/10., r) !(tag, V, Rxy, Lz, r)

  print *, "Program End."

end program MechTest