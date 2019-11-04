      FUNCTION ALPHAS(RQ)
      IMPLICIT NONE
#include "caspara.inc"
#include "casludat1.inc"
#include "casvalues.inc"
#include "casscal.inc"
	Double Precision ALPHAS,RQ,RQ1,PYALPS,ALAM
      Double Precision QG0
	Integer ikincut
      COMMON /GLUDAT/QG0,ikincut
      Double Precision SCALE,Qgcut
	Double Precision Qscal
	Integer Iqqbar
	Common/GLUDAT2/Qscal,Iqqbar
      Logical First/.True./
      Double Precision scale_mz
      RQ1 = RQ
      IF(IRUNA.EQ.0) THEN
         ALPHAS=ALPHS
      ELSE

         ALAM=PARU(112)
	   QGcut=max(ALAM,Qg0)
         SCALE = RQ1*RQ1
	   SCALE = SCALE * SCALFA * Qscal
         IF(SCALE.LT.Qgcut**2) SCALE=Qgcut**2
c         if(sqrt(scale).lt.1.5)  	
c     &      write(6,*) 'lambda ',alam,qgcut,sqrt(scale),PYALPS(SCALE)
         ALPHAS =PYALPS(SCALE)
         
      If(First) then
         First=.False.
         scale_mz = 91.8*91.8
         
         write(6,*) ' alphas order: ', MSTU(111),PARU(112), '  alphas(M_Z) = ', PYALPS(scale_mz)
      
      endif
	
      ENDIF
      RETURN
      END
