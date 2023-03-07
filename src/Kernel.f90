module Kernels

use ModPrecision, only : fp_kind

real(fp_kind), parameter :: m_pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406
contains

  real(fp_kind) function Kernel(q, h)
    real(fp_kind), intent(in) :: q,h
		real(fp_kind) :: C

    !Qubic Spline
    !Dim == 2 ? C = 10.0/(7.0*h*h*M_PI) : C = 1.0/(h*h*h*M_PI);
    C = 1.0/(h*h*h*m_pi)
    if 		(q < 1.0)	then 
    Kernel = C*(1.0-(3.0/2.0)*q*q+(3.0/4.0)*q*q*q)
    else if (q < 2.0)	then
    Kernel = C*((1.0/4.0)*(2.0-q)*(2.0-q)*(2.0-q))
    else						
    Kernel = 0.0
    end if
  end function Kernel

	real(fp_kind) function GradKernel(q, h)
    real(fp_kind), intent(in) :: q,h
		real(fp_kind)::C
    C = 1.0/(h*h*h*h*m_pi)
				!Dim ==2 ? C = 10.0/(7.0*h*h*h*M_PI) : C = 1.0/(h*h*h*h*M_PI);

    if 		  (q==0.0)	then
      GradKernel = C/h    *(-3.0+(9.0/2.0)*q)
    else if (q<1.0)		then 
      GradKernel = C/(q*h)*(-3.0*q+(9.0/4.0)*q*q)
    else if (q<2.0)		then 
      GradKernel = C/(q*h)*((-3.0/4.0)*(2.0-q)*(2.0-q))
    else							
      GradKernel = 0.0  
    end if 
  end function GradKernel
  
end module Kernels