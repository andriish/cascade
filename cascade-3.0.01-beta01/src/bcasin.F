      SUBROUTINE BCASIN

      IMPLICIT NONE
      CHARACTER*23  FLAG,CHECK
#include "caspar.inc"

C...Give sensible default values to switches and parameters.
      INTEGER I
      DOUBLE PRECISION EPARL
      Double Precision X,Y,Q2
      DOUBLE PRECISION PARL
      Integer LST
      COMMON /CASTOU/PARL(30),X,Y,Q2,LST(40)

      INTEGER KLST,KIPY
      Double Precision EPYPAR,EPYVAR
      DIMENSION KLST(40),EPARL(30)
      DIMENSION KIPY(80),EPYPAR(80),EPYVAR(80)
C...CASTOU: Cuts, basic switches and parameters.
      DATA KLST/0,1,5,1,3,1,1,12,5,1,0,4,5,1,1,1,0,2,3,21*0/
**      DATA EPARL/1.,1.,0.44,0.75,0.226,0.,0.,0.015,2.,0.,0.01,4.,
**     +0.1,0.44,0.01,7.29735E-03,1.16637E-05,0.044,0.03,1.,10*0./
      DATA EPARL/1.,0.7,0.44,0.75,0.226,0.,0.,0.015,2.,0.,0.01,4.,
     +0.1,0.44,0.01,7.29735E-03,1.16637E-05,0.044,0.03,1.,10*0./
      DATA KIPY/
     + 0,     0,     2,     2,     6,     1,     1,     4,     3,     1,
     + 3,     1,     1,     2,     1,     1,     4,     1,     1,     1,
     + 0,     1,     1,     1,     1,     1,     1,     0,     0,     0,
     + 1,     2,     1,     1,    30,    33,     1,     1,     7,     0,
     + 0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     + 0,     0,     0,     1,   100,     0,     0,     0,     0,     0,
     + 0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     + 0,     0,     0,     0,     0,     0,     0,     0,     0,     0/
      DATA (EPYPAR(I),I=1,40)/
     +   7.299E-03,   2.290E-01,   2.000E-01,   2.500E-01,   4.000E+00,
     +   1.000E+00,   4.400E-01,   4.400E-01,   7.500E-02,   0.000E+00,
     +   2.000E+00,   2.000E+00,   1.000E+00,   0.000E+00,   3.000E+00,
     +   1.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   1.000E+00,
     +   2.500E-01,   1.000E+00,   2.000E+00,   1.000E-03,   1.000E+00,
     +   1.000E+00,   1.000E+00,  -2.000E-02,  -1.000E-02,   0.000E+00,
     +   0.000E+00,   1.600E+00,   0.500E+00,   0.200E+00,   3.894E-01,
     +   1.000E+00,   3.300E-01,   6.600E-01,   0.000E+00,   1.000E+00/
      DATA (EPYPAR(I),I=41,80)/
     +   2.260E+00,   1.000E+04,   1.000E-04,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     +   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00/
      DATA EPYVAR/80*0./
      DATA CHECK /'CASCADE 2.0 INITIALISED'/
      DATA FLAG  /'                      '/

      IF (FLAG .EQ. CHECK) THEN
         WRITE(6,*) ' BLOCK DATA for CASCADE 2.0 REINITIALISED.'
      ENDIF

      CALL UCOPY(KLST,LST,40)
      CALL UCOPY(EPARL,PARL,30)
      CALL UCOPY(KIPY,IPY,80)
      CALL UCOPY(EPYPAR,PYPAR,80)
      CALL UCOPY(EPYVAR,PYVAR,80)
      FLAG = CHECK
      WRITE(6,*) ' BLOCK DATA for CASCADE  initialised'

      RETURN
      END
