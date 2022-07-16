module Geometry
use ParticleData

implicit none 
contains 

subroutine AddCylinderLength(tag, V, pt)

		zp = V(2)+r;
		//Calculate row count for non ghost particles
		do while (zp <= (V(2)+Lz -r)){
			k++; 
      zp += 2.*r;      
		}
    
end subroutine AddCylinderLength

end module Geometry