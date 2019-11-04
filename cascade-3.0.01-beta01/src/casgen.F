      SUBROUTINE casgen(NDIMEN,XGR)
      IMPLICIT NONE
      Double Precision ACC1,ACC2
	Integer IINT,NCB
      COMMON /INTEGR/ ACC1,ACC2,IINT,NCB
	Double Precision XGR
	Double Precision XGS,FXNB
	Integer NDIMEN,MXTRY,I
      DIMENSION XGR(20)
      COMMON/XFXNB/XGS(20)
      EXTERNAL FXNB
      IF(iint.eq.0) then
         MXTRY = 500
         CALL SPRING( FXNB, MXTRY )
      elseif(iint.eq.1) then
         MXTRY = 10000
         CALL MCGEN( FXNB, MXTRY )
      else
      write(6,*)' CASCADE: IINT = ',IINT,' not implemented: STOP '
      endif
      DO I=1,NDIMEN
         XGR(I) =XGS(I)
      ENDDO
      RETURN
      END

