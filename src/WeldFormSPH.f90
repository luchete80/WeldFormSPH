program WeldFormSPH
use Integrate
use ParticleData
use Domain 
use Neighbor 
use Thermal
use omp_lib
use Thermal

implicit none


!Type(Particle), POINTER :: pt

  real, dimension(1:3) :: V
  real :: dx, r, Rxy, Lz, h 
  integer:: i, tnr, maxt
  real,allocatable, dimension(:):: dTdt
  real :: t_, deltat
  
  call omp_set_num_threads(4);
  
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

  call DomInit(4)
  call InitNb()
  call AllocateNbData()
  
  call CellInitiate()
  call ListGenerate()
  
  call MainNeighbourSearch()
  call InitRedArraysOnce()
  call CalcPairPosList()
  
  allocate (dTdt(part_count))
  
  pt%t(:)     = 20.
  pt%cp_t(:)  = 960.
  pt%k_t(:)   = 120.  
  pt%t(900:1000) = 100.
  
  deltat = 0.01
  t_ = 0.
  do while (t_ < 1.)
    call CalcTempIncPart(dTdt)
    print *, dTdt(900:1000)
    pt%t(:) = pt%t(:) + dTdt(:) * deltat
    
    t_ = t_ + deltat
  end do
  !print *, "Temperatures "
  !do p = 0, part_count
   ! print *, pt%t(:) 
  !end do
!open (12,file='temp.log', status='old', position='APPEND')  
  print *, "Program End."
end program WeldFormSPH

