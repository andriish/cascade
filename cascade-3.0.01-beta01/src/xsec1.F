      Subroutine xsec1(X,XSEC)
      IMPLICIT None
	Double Precision X(20)
#include "caslujets.inc"
#include "casparas.inc"
#include "caspart.inc"
#include "casluco.inc"
      Integer IDEBUG
      COMMON/ INTERN/IDEBUG
#include "caspara5.inc"
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
	
c	IDEBUG = 1
c      write(6,*) KBEAM(1,2),KBEAM(2,2)
c      write(6,*) ' calling meoffsh in xsec1 '

      IF((KBEAM(1,2).eq.22.or.Iabs(KBEAM(1,2)).eq.11).and.
     & 	(KBEAM(2,2).eq.22.or.Iabs(KBEAM(2,2)).eq.11)) Then
         CALL P_gamgam(X,WPART)
c	   write(6,*) ' gamgam called'
	Else
         CALL P_SEMIH(X,WPART)
c	   write(6,*) ' standard called', wpart
	Endif
      IF(WPART.GT.0.) THEN
         CALL MEOFFSH(WT1)
c         write(6,*) 'after meoffsh ',wt1
      ENDIF
      WEIG=WPART*WT1*GEV2NB
C BOOST BACK TO OVERALL CMS
      IF(WEIG.EQ.0.0) Return
      ncall = ncall + 1
c      call PYLIST(1)
      XSEC=WEIG
      RETURN
      END
