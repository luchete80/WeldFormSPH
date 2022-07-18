program WeldFormSPH
use Integrate
use ParticleData
use Domain 
use Neighbor 

use omp_lib

implicit none


!Type(Particle), POINTER :: pt

  real, dimension(1:3) :: V
  real :: dx, r, Rxy, Lz, h 
  integer:: i, tnr, maxt
  
  !$ call omp_set_num_threads(4);
  
  maxt = omp_get_max_threads()
  write( *, * ) 'Max threads ', maxt
  
  dx    = 0.1
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
  !$omp parallel private(tnr) 
  tnr = omp_get_thread_num()
  do i = 1, 10000
     !write( *, * ) 'Thread', tnr, ':'!,  i
     print *, "Hello thread ",tnr  
  end do
  !$omp end parallel   

  ! !$omp parallel 
  ! !$omp do private(tnr)
  ! do i = 1, 20
     ! tnr = omp_get_thread_num()
     ! write( *, * ) 'Thread', tnr, ':',  i
  ! end do
  ! !$omp end do
  ! !$omp end parallel

  ! !$omp parallel
  ! !$omp do private(tnr)
     ! tnr = omp_get_thread_num()
     ! write( *, * ) 'Thread', tnr, ':',  i

  ! !$omp end parallel
  call DomInit(4)
  call InitNb()
  call CellInitiate()
  call ListGenerate()

!open (12,file='temp.log', status='old', position='APPEND')  
  print *, "Program End."
end program WeldFormSPH

