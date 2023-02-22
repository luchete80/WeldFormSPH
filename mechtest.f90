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
  real(fp_kind):: r, Lz, Rxy
  
  call AddCylinderLength(0, V, Rxy, Lz, r) !(tag, V, Rxy, Lz, r)

end program MechTest