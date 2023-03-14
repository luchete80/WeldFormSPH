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
  real(fp_kind):: r, Lx, Ly, Lz, Rxy, dx, h, e, nu, k_i, g_i
  
  real(fp_kind):: dt, t
  !! MATERIAL
  real(fp_kind)::rho
  integer :: i
  
  
  
  nproc = 1
  
  Lz = 0.56
  Rxy = 0.15
  
  rho	= 2700.0;
  dx = 0.002;
  h	= dx*1.2 
  r = dx/2.

  e = 70.0e9
  nu = 0.330
  
  k_i= E / ( 3.*(1.-2*nu) );
  g_i= E / (2.* (1.+nu));  
  
  mat_G = g_i
  
  ! CALLING NBS TO THIS DOMAIN  STALLS THE PROGRAM
  !call AddCylinderLength(0, V, Rxy, Lz + 2.0 * Lz/10., r, rho, h) !(tag, V, Rxy, Lz, r)
  Lx = 0.1
  Ly = 0.024
  Lz = 0.012	  
  
  call AddBoxLength(0, V, Lx + Lx/20., Ly, Lz, r, rho, h)
  
  pt%cs(:)  = sqrt(K_i/rho)
  
  !Boundary zones
  do i=1,part_count
    if (pt%x(i,1) < Lx/40.) then
      pt%id(i) = 2
    end if
    if (pt%x(i,1) > Lx+Lx/40.) then
      pt%id(i) = 3
    end if
  end do

  dt = 0.4 * h / pt%cs(1)
  print *, "h ", h, "Cs ", pt%cs(1), "Time step: ", dt
  
  !t = 1.0e-4
  t = 2.*dt
  call SolveDiffUpdateFraser(t,dt)
  
  call CalcEquivalentStress()
    
  open (1,file='test.csv')!, position='APPEND')  
  write (1,*) "X, Y, Z, id, rho, Nb, Vx, Vy, Vz, Ax, Ay, Az, &
              &Ux, Uy, Uz, sig_eq"
  
  do i=1,part_count  
    write (1,*) pt%x(i,1), ", ", pt%x(i,2), ", " ,pt%x(i,3), ", " ,pt%id(i), ", ", pt%rho(i), ", ", ipair_t(i)+jpair_t(i), ", " ,&
                &pt%v(i,1), ", ", pt%v(i,2), ", ", pt%v(i,3), ", ",&
                &pt%a(i,1), ", ", pt%a(i,2), ", ", pt%a(i,3), ", ",&
                &pt%disp(i,1), ", ", pt%disp(i,2), ", ", pt%disp(i,3), ", ",&
                &pt%sigma_eq(i)
  end do
  close(1)
  
  print *, "Program End."

end program MechTest