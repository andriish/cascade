*CMZ :  1.00/08 14/10/2001  16.41.42  by  Hannes Jung
*CMZ :  1.00/07 28/09/2001  20.49.12  by  Hannes Jung
*CMZ :  1.00/00 12/06/2001  12.05.04  by  Hannes Jung
*-- Author :
      FUNCTION XGLU(XX,QQ2)
      Implicit none
#include "casglu.inc"
      DOUBLE PRECISION XX,QQ2,XGLU,QQ22,PP2
      DOUBLE PRECISION XXX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU,
     +VAL(20)
      LOGICAL PDFFIRST
      COMMON /W50516/PDFFIRST
      CHARACTER*20 PARM(20)
      Integer NPDF,IP2
      DATA NPDF/0/,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU/9*0.D0/
      QQ=SQRT(MAX(0.D0,QQ2))
      XXX = XX

ccc for LHAPDF 5.7            
      PARM(1) = 'DEFAULT'
      VAL(1) = icolglu
CNEW
      NPDF=NPDF+1
      PDFFIRST=.FALSE.
      IF(NPDF.LE.1) THEN
         PDFFIRST=.TRUE.
      ENDIF
c call PDFSET each time
      CALL PDFSET(PARM,VAL)
	IF(VAL(1).EQ.3) THEN
         QQ22=QQ**2
         PP2 =0d0
	   IP2=0
         CALL STRUCTP(XXX,QQ22,PP2,IP2,UPV,DNV,USEA,DSEA,STR,CHM,
     &	   BOT,TOP,GLU)
	ELSE
         CALL STRUCTM(XXX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
	ENDIF
      XGLU = GLU
      RETURN
      END
