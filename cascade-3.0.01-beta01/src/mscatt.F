      FUNCTION MSCATT(X,WEIGHT,ierr)

      IMPLICIT none
      Double Precision X(20)
#include "caslujets.inc"
#include "casparas.inc"
#include "caspart.inc"
#include "casluco.inc"
#include "casmult.inc"
#include "caspar.inc"

	Double Precision xg1,xg1_bar,yx1,yx1_bar 
	Double Precision xg2,xg2_bar,yx2,yx2_bar 
	Common /resolv/xg1,xg1_bar,yx1,yx1_bar,xg2,xg2_bar,yx2,yx2_bar
      DOUBLE PRECISION WEIGL,WEIGC,WEIGB,WEIGT
      COMMON /weiglh/weigl,weigc,weigb,weigt
      Integer IGENFL,IMEPHI
      COMMON/GENWEI/IGENFL,IMEPHI
      Integer IPROO
      Double Precision WEIG,MSCATT,XSEC,WEIGHT,pt2cut10,pt2cut11
	Integer NCALL,I,J
	Double Precision GEV2NB,WPART,WT1
      Integer Nhs
      Double Precision DOT,DOT1
	Double Precision CMS(4)
      Double Precision  XGR(20)
	Integer Len
      DATA NCALL/0/
      DATA GEV2NB/.3893857E+6/
      Integer IDEBUG
      COMMON/ INTERN/IDEBUG
	Integer IDEBUGO
	Integer Ierr
      IMEPHI = 0
      WEIG=0.D0
      MSCATT=0.D0
      IPROO = IPRO
	IDEBUGO = IDEBUG
	IDEBUG = 0
      pt2cut10 = pt2cut(10)
      pt2cut11 = pt2cut(11)
	pt2cut(10) = 0.4
c        pt2cut(10) = 1.
c	pt2cut(11) = 0.
c	write(6,*) ' mscatt pt2cut ',pt2cut(10),pt2cut(11)
      Weig = 0.D0
      Weigl = 0.D0
      Weigc = 0.D0
      Weigb = 0.D0
      Weigt = 0.D0
      Xsec = 0.D0
	do I=1,100
	   do J=1,5
	      K(I,J) = 0
	      P(I,J) = 0
	   enddo
	enddo
c	call pyedit(12)
c give temporay beam vectors for reduced beam energy
c      write(6,*) ' mscatt ',pyvar(31),pyvar(32),plepin,ppin
c      write(6,*) ' mscatt p1:',(1-pyvar(31))*plepin
c      write(6,*) ' mscatt p2:',(1-pyvar(32))*ppin
      DO J=1,5
      K(1,J)=KBEAM(1,J)
      K(2,J)=KBEAM(2,J)
      Enddo
      P(1,1) = 0.0D0
      P(1,2) = 0.0D0
      P(1,3) = (1d0-pyvar(31))*PLEPIN
      P(1,5) = PYMASS(KE)
      P(1,4) = DSQRT(P(1,1)**2+P(1,2)**2+P(1,3)**2+P(1,5)**2)
C...  GIVE PROTON FOUR VECTOR
      P(2,1)= 0.0D0
      P(2,2)= 0.0D0
      P(2,3)=(1d0-pyvar(32))*PPIN
      P(2,5)=PYMASS(KP)
      P(2,4)=DSQRT(P(2,1)**2+P(2,2)**2+P(2,3)**2+P(2,5)**2)
	
      DO 20 I=1,4
         CMS(I) =P(1,I)+P(2,I)
   20 CONTINUE
      N = 2
C BOOST TO EP CMS
c	write(6,*) 'mscatt: before secondary interaction N= ',N
c	call pylist(1)
      CALL PYROBO(0,N,0.D0,0.D0,-CMS(1)/CMS(4),-CMS(2)/CMS(4),
     &	-CMS(3)/CMS(4))
      SSS = DOT(CMS,CMS)
      if(sss.lt.10.) then
	  ierr = 1
	  goto 345
	endif
c	write(6,*) 'mscatt: 1st after boost back to cms ' 
c	call pylist(1)

	IPRO = 10
	
      CALL p_mscatt(X,WPART)
      WT1 = 0
      IF(WPART.GT.0.) THEN
         CALL MEOFFSH(WT1)
c	   WT1 = 1
      ENDIF
      Weig=WPART*WT1*GEV2NB

c	write(6,*) ' after secondary interaction  N= ',N
c	write(6,*) ' mscatt: Weig = ',weig,wpart,wt1,(X(i),I=1,8)
c	call pylist(1)
cc	Xsec=Weig
c      Weigl = Xsec
c      Weig = Xsec
      ncall = ncall + 1
	goto 345

345   Continue	
      IMEPHI = 0
      IPRO = IPROO
	IDEBUG = IDEBUGO
	pt2cut(10)=pt2cut10
	pt2cut(11)=pt2cut11
      mscatt=WEIG
c	write(6,*) 'mscatt: end before boost back to cms ' 
c	call pylist(1)
      CALL PYROBO(0,N,0.D0,0.D0,CMS(1)/CMS(4),CMS(2)/CMS(4),
     &	CMS(3)/CMS(4))
      WEIGHT = WEIG
c	write(6,*) 'mscatt: after boost back to cms ' 
c	call pylist(1)
      RETURN
      END
