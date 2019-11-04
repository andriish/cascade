      Subroutine Xv1(X,XSEC)
      IMPLICIT None
      Double Precision X(20)
#include "caslujets.inc"
#include "casparas.inc"
#include "caspart.inc"
#include "casluco.inc"
	Integer NCALL,I,J
	Double Precision GEV2NB,XSEC,WEIG,WPART,WT1
      DATA NCALL/0/
      DATA GEV2NB/.3893857E+6/
      DO 10  I=1,20
         DO 10 J=1,5
            K(I,J) = 0
   10 P(I,J) = 0.0D0
      N=0
      XSEC=0.D0
      WEIG=0.D0
      WPART=0
c	write(6,*) ' in xvq '
      CALL P_SEMIH(X,WPART)
      IF(WPART.GT.0.) THEN
         CALL mev(WT1)
c         write(6,*) ' after mevq ',WPART,WT1,IPRO
      ENDIF	
c      write(6,*) WPART,WT1,IPRO
      WEIG=WPART*WT1*GEV2NB
C BOOST BACK TO OVERALL CMS
      IF(WEIG.EQ.0.0) Return
      ncall = ncall + 1
c      call PYLIST(1)
      XSEC=WEIG
      RETURN
      END
