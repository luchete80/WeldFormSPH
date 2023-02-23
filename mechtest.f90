program MechTest

use Integrate
use ParticleData
use Domain 
use Neighbor 
use Thermal
use omp_lib
use Thermal
use Mechanical
use Solver
use ModPrecision, only : fp_kind

implicit none

  real(fp_kind), dimension(1:3) :: V
  real(fp_kind):: r, Lz, Rxy, dx, h 
  
  real(fp_kind):: dt, t
  Lz = 0.56
  Rxy = 0.15
  
  dx = 0.015
  h	= dx*1.2 
  r = dx/2.
  
  !! MATERIAL
  real(fp_kind)::rho
  
  rho = 1.
  
  !call AddCylinderLength(0, V, Rxy, Lz + 2.0 * Lz/10., r) !(tag, V, Rxy, Lz, r)
  
  call AddBoxLength(0, V, Lz, Lz, Lz, r, rho, h)


  dt = 1.e-4
  t = 0.01
  call SolveDiffUpdateFraser(t,dt)
  
  
  print *, "Program End."

end program MechTest