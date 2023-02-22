module Functions
  use ModPrecision, only : fp_kind
	!real(fp_kind),
  contains

  real(fp_kind) function SoundSpeed(EQ, Cs0, Density, Density0)
		!switch (EQ)
		!{
		!	case 0:
				SoundSpeed = Cs0;
				! break;

			! case 1:
				! return sqrt((Cs0*Cs0)*pow(Density/Density0,6.0));
				! break;

			! case 2:
				! return Cs0;
				! break;

			! default:
				! std::cout << "Please correct Pressure Equation No and run again" << std::endl;
				! std::cout << "0 => P0+(Cs*Cs)*(Density-Density0)" << std::endl;
				! std::cout << "1 => P0+(Density0*Cs*Cs/7)*(pow(Density/Density0,7)-1)" << std::endl;
				! std::cout << "2 => (Cs*Cs)*Density" << std::endl;
				! abort();
				! break;
  
  end function SoundSpeed


end module Functions