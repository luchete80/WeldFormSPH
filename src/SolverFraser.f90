module Solver
use ModPrecision, only : fp_kind

contains 

subroutine SolveDiffUpdateFraser (tf, dt)
  
  use omp_lib
  use Domain
  use Mechanical !Tensors, accel
  use Neighbor
  use Kernels
  implicit none
 
 real(fp_kind)::tf, dt
 
  real(fp_kind),dimension(3)::du
  integer :: i

  call omp_set_num_threads(12); 
  call DomInit(12)
  
  call InitNb()
  call AllocateNbData()
  
  call CellInitiate()
  call ListGenerate()
  
  call MainNeighbourSearch()
  call InitRedArraysOnce()
  call CalcPairPosList()
  
  do while (tf <= 100.0*dt)
    call CalcDensIncPart
    call CalcRateTensorsPart
    call CalcAccelPart
    
    !$omp parallel do num_threads(Nproc) private (du)
    do i = 1, part_count
    du = pt%v(i,:) * dt * + 0.5 * pt%a(i,:) * dt * dt 
    pt%x(i,:)     = pt%x(i,:)   + du
    pt%disp(i,:)  = pt%disp(i,:) + du
    end do
    !$omp end parallel do  

    !$omp parallel do num_threads(Nproc)  
    do i = 1, part_count
      pt%v(i,:) = pt%a(i,:) * dt 
    end do
    !$omp end parallel do  

  end do 
  

end subroutine SolveDiffUpdateFraser

end module Solver