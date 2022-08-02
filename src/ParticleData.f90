
! TYPE Struct Of Arrays (SOA), for better GPU performance
Module ParticleData 

Type Particle
  Integer, Dimension(:), Allocatable :: ID
  real, dimension(:,:), Allocatable :: x
  !GENERAL
  real, dimension(:), Allocatable :: h, t, cs, rho, m !influence radius, temp
  !THERMAL
  real, dimension(:), allocatable :: cp_t, k_t
  !Mechanical
  !real, dimension(:), allocatable :: cs
  real, dimension(:,:), allocatable :: v
  
End Type

End Module ParticleData