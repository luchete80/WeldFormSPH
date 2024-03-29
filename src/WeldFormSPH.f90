program WeldFormSPH
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

!Type(Particle), POINTER :: pt

  real(fp_kind), dimension(1:3) :: V
  real(fp_kind) :: dx, r, Rxy, Lz, h 
  integer:: i, tnr, maxt ,nn
  real(fp_kind),allocatable, dimension(:):: dTdt
  real(fp_kind) :: t_, deltat
  real(fp_kind) :: start, finish
  real(fp_kind) :: L, rho 

   ! Variables for clock
   integer count_0, count_1
   integer count_rate, count_max
   double precision time_init, time_final, elapsed_time

    
  call omp_set_num_threads(12);
  
  maxt = omp_get_max_threads()
  write( *, * ) 'Max threads ', maxt
  
  dx    = 0.02
  Rxy  = 0.15
  Lz = 0.56
  r = dx / 2.0
  h = dx * 1.2

  V(1) = 0.;V(2) = 0.;V(3) = 0.
  !AddBoxLength(tag, V, Lx, Ly, Lz, r, Density,  h)		
  L = 1.	
  rho = 1000.
  call AddBoxLength(0, V, L, L, L, r, rho, h)

  do i = 1, part_count
  !write (*,*) "Particle", i ," position is ", pt%x(i,1), pt%x(i,1), pt%x(i,3)
  end do 
  !call AddCylinderLength(0, V, Rxy, Lz, r)

  call DomInit(12)
  call InitNb()
  call AllocateNbData()
  
  call CellInitiate()
  call ListGenerate()
  
  call MainNeighbourSearch()
  call InitRedArraysOnce()
  call CalcPairPosList()
  
  allocate (dTdt(part_count))
  allocate (temp_pair(pair_tot_count))
  pt%t(:)     = 20.
  pt%cp_t(:)  = 1.
  pt%k_t(:)   = 3000.  
  
  print *, "Size of floating point: ", sizeof(pt%t(1))

  !call cpu_time(start)
  !deltat = 0.00036
  deltat = 0.3*h*h*rho*pt%cp_t(1)/pt%k_t(1)	
  
  print *,'Reduction by particle... '    
  t_ = 0.

   ! Starting time
   call system_clock(count_0, count_rate, count_max)
   time_init = count_0*1.0/count_rate

  
  nn=(L/dx)*(L/dx)
    pt%t(1:nn) = 500.
  !!! FASTEST ALGORITHM
  do while (t_ <= 1000.0*deltat)

    call CalcTempIncPart(dTdt) !THIS IS THE FASTEST BY NOW
    !print *, "dTdt 0",  dTdt(1)
    pt%t(:) = pt%t(:) + dTdt(:) * deltat
    
    t_ = t_ + deltat
  end do

  ! Ending time
  call system_clock(count_1, count_rate, count_max)
  time_final = count_1*1.0/count_rate
  ! Elapsed time
  elapsed_time = time_final - time_init

  ! Write elapsed time
  write(*,1003) int(elapsed_time),elapsed_time-int(elapsed_time)
  !write(*,*) "Elapsed time: ", int(elapsed_time),elapsed_time-int(elapsed_time)
   
  !call cpu_time(finish)
  print *,'Time: ', t_ 
  print '("CPU Time = ",f6.3," seconds.")',finish-start
  print *, "Program End."

  
  open (1,file='temp.csv')!, position='APPEND')  
  write (1,*) "X, Y, Z, temp"

  do i=1,part_count  
    write (1,*) pt%x(i,1), ", ", pt%x(i,2), ", " ,pt%x(i,3), ", " ,pt%t(i) 
  end do
  close(1) 


  ! print *,'Reduction by pair... '
  ! call cpu_time(start)
  ! t_ = 0.
  ! do while (t_ <= 2.0)
    ! pt%t(1:400) = 500.
    ! call CalcTempInc()
    ! call TempReduction(dTdt)
    ! !print *, "dTdt 0",  dTdt(1)
    ! pt%t(:) = pt%t(:) + dTdt(:) * deltat
    
    ! t_ = t_ + deltat
  ! end do
  ! call cpu_time(finish)
  
  ! open (1,file='temp.csv')!, position='APPEND')  
  ! write (1,*) "X, Y, Z, temp"
  ! print *,'Time: ', t_ 
  ! print '("CPU Time = ",f6.3," seconds.")',finish-start
  ! print *, "Program End."
  ! do i=1,part_count  
    ! write (1,*) pt%x(i,1), ", ", pt%x(i,2), ", " ,pt%x(i,3), ", " ,pt%t(i) 
  ! end do
  ! close(1) 

  1003 format('  Elapsed Time  = ',i0,f0.9)
  
end program WeldFormSPH

