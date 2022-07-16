
! TYPE Struct Of Arrays (SOA), for better GPU performance
Module ParticleData 

Type Particle
  Integer, Dimension(:), Allocatable :: ID
  real, dimension(:), Allocatable :: x
  real, dimension(:), Allocatable :: y
  real, dimension(:), Allocatable :: z
End Type

End Module ParticleData