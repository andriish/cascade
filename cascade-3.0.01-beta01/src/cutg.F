      SUBROUTINE CUTG(PT2,NACC)
      IMPLICIT NONE
	Double Precision PT2
	INTEGER NACC
#include "casparas.inc"
      NACC=1
      IF(PT2.LT.PT2CUT(IPRO)) NACC=0
      RETURN
      END
