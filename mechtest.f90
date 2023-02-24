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
  real(fp_kind):: r, Lx, Ly, Lz, Rxy, dx, h 
  
  real(fp_kind):: dt, t
  !! MATERIAL
  real(fp_kind)::rho
  integer :: i
  
  Lz = 0.56
  Rxy = 0.15
  
  rho	= 2700.0;
  dx = 0.002;
  h	= dx*1.2 
  r = dx/2.
  
  
  ! CALLING NBS TO THIS DOMAIN  STALLS THE PROGRAM
  !call AddCylinderLength(0, V, Rxy, Lz + 2.0 * Lz/10., r, rho, h) !(tag, V, Rxy, Lz, r)
  Lx = 0.1
  Ly = 0.024
  Lz = 0.012	  
  call AddBoxLength(0, V, Lx, Ly, Lz, r, rho, h)


  dt = 1.e-4
  t = 0.01
  call SolveDiffUpdateFraser(t,dt)
  
  open (1,file='test.csv')!, position='APPEND')  
  write (1,*) "X, Y, Z, temp"

  do i=1,part_count  
    write (1,*) pt%x(i,1), ", ", pt%x(i,2), ", " ,pt%x(i,3), ", " ,pt%t(i) 
  end do
  close(1)
  
  print *, "Program End."

end program MechTest