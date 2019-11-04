      FUNCTION DOT1(I,J)
	IMPLICIT NONE
#include "caslujets.inc"
      Double Precision DOT1
	INTEGER I,J
C+++++++++++
C    DOT PRODUCT OF FOUR VECTOR IN MINKOWSKI METRIK
C    WITH VECTORS FROM LUJETS
C+++++++++++
      DOT1= P(I,4)*P(J,4)-(P(I,1)*P(J,1))-(P(I,2)*P(J,2))-
     .     (P(I,3)*P(J,3))
      RETURN
      END
