module Thermal
use Neighbor
use Domain

contains

  subroutine CalcTempInc(dTdt)
    real, intent(out)::dTdt(part_count)
  
  end subroutine CalcTempInc

end module Thermal