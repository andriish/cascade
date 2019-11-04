      SUBROUTINE CAPRIKT(S,PT,PHI)
	IMPLICIT NONE
      Double Precision dcasrn
	Double Precision S,ST,PT,PHI
      EXTERNAL dcasrn
C...Size (PT) and azimuthal angle (PHI) of primordial kt according
C...to a Gaussian distribution.
c	S=1.4d0
      ST = S
      PT=ST*SQRT(-DLOG(dcasrn()))
      PHI=6.2832d0*dcasrn()
      RETURN
      END
