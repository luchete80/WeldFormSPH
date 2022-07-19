module Kernels

real, parameter :: m_pi = 3.1415926
contains

  real function Kernel(q, h)
  
		real :: C

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

	real function GradKernel(q, h)
		real::C
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