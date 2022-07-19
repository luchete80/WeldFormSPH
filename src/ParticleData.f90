
! TYPE Struct Of Arrays (SOA), for better GPU performance
Module ParticleData 

Type Particle
  Integer, Dimension(:), Allocatable :: ID
  real, dimension(:,:), Allocatable :: x
  real, dimension(:), Allocatable :: h, t, cs, rho, m !influence radius, temp
  
End Type

End Module ParticleData