      FUNCTION AALAM(A,B,C)
      IMPLICIT NONE
	double precision A,B,C,AALAM
      AALAM = A**2 + B**2 + C**2 - 2.D0*A*B - 2.D0*B*C - 2.D0*C*A
      RETURN
      END
