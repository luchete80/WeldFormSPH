
! TYPE Struct Of Arrays (SOA), for better GPU performance
Module ParticleData 
  use ModPrecision, only : fp_kind

  public :: bc_type, bc_vel, bc_disp
  enum, bind(c)
    enumerator :: bc_type = 0
    enumerator :: bc_vel = 157839
    enumerator :: bc_disp = 23097
  end enum
  !Usage 
  ! integer(kind(bc_type)) :: bc
  ! bc = bc_disp  

Type Particle
  Integer, Dimension(:), Allocatable :: ID
  real(fp_kind), dimension(:,:), Allocatable :: x
  !GENERAL
  real(fp_kind), dimension(:), Allocatable :: h, t, cs, rho, m, rho_0 !influence radius, temp

  !THERMAL
  real(fp_kind), dimension(:), allocatable :: cp_t, k_t
  !Mechanical
  !real(fp_kind), dimension(:), allocatable :: cs
  real(fp_kind), dimension(:,:), allocatable :: v, a, disp
  real(fp_kind), dimension(:,:,:), allocatable :: sigma, str_rate, rot_rate, shear_stress
  real(fp_kind), dimension(:), allocatable:: pressure, strain, mat_g
  
  Integer solver_type
 
  
End Type

End Module ParticleData