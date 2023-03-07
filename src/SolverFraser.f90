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
 
  real(fp_kind),intent(in)::tf, dt
 
  real(fp_kind),dimension(3)::du
  
  real(fp_kind),dimension(3)::dumax
  
  integer :: i, j

  call omp_set_num_threads(12); 
  call DomInit(12)
  
  call InitNb()
  call AllocateNbData()
  
  call CellInitiate()
  call ListGenerate()
  
  call MainNeighbourSearch()
  call InitRedArraysOnce()
  call CalcPairPosList()
 
  
  time = 0.
  do while (time < tf)
    call StartVars !Similar to start acceleration in weldform
    !APPLY BC!
    do i=1,part_count
      if (pt%id(i) == 2 ) then
        pt%v(i,:) = 0.0
      end if
      if (pt%id(i) == 3) then
        !print *, "id ", i 
        pt%v(i,1) = 0.0
        pt%v(i,2) = -0.48
        pt%v(i,3) = 0.0
      end if
    end do 
    
    call CalcDensIncPart
    !$omp parallel do num_threads(Nproc) private (du)
    do i = 1, part_count
      pt%rho(i) = pt%rho(i) + pt%drhodt(i) * dt
    end do
    !$omp end parallel do  
    
    print *, "dens 52",  pt%rho(52)
    !print *, "dens 674",  pt%rho(674)
  ! do i = 1, part_count
    ! print *, "dens ",  pt%rho(i)
  ! end do
    
    call CalcRateTensorsPart
    print *, "str rate 52",  pt%str_rate(52,:,:)
    call CalcStressStrain(dt)

    print *, "Sigma 103 " , pt%sigma(1298,:,:)
    
    call CalcAccelPart
    print *, "Accel 52" , pt%a(52,:)
    ! !REINFORCE bc vel
    do i=1,part_count
       if (pt%id(i) == 2 .or. pt%id(i) == 3) then
        pt%a(i,:) = 0.
       end if
    end do
    
    !$omp parallel do num_threads(Nproc) private (du)
    do i = 1, part_count
    !print *, "acc ",  pt%a(i,1), ", ",  pt%a(i,2), ", ", pt%a(i,3)
    !    print *, "StrainRate ",  pt%str_rate(i,:,:)
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

    !!REINFORCE bc vel again
    do i=1,part_count
      if (pt%id(i) == 2 ) then
        pt%v(i,:) = 0.
      end if
      if (pt%id(i) == 3) then
        !print *, "id ", i 
        pt%v(i,1) = 0.
        pt%v(i,2) = -0.48
        pt%v(i,3) = 0.
      end if
    end do    

      print *, "dens ",  pt%rho(1)
    
    time = time + dt
  end do 
  
  dumax = 0.0
  do i=1,part_count  
    do j=1,3
      if (abs(pt%disp(i,j)) > dumax(j) ) then
        dumax(j) = pt%disp(i,j)
      end if
    end do
  end do 
  print *, "Max Displacements ", dumax
  
end subroutine SolveDiffUpdateFraser

end module Solver