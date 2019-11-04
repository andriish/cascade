      SUBROUTINE cauniglu(KF,XI,KT2I,TI,XPQ)
      Implicit None
      Integer kf,i,Iccfmset
      Double Precision  XI,KT2I,TI
      Double Precision  XPQ(-6:6)
      Double Precision XPQS(-6:6),XPQQ(-6:6)
      Double Precision X,KT2,T,scale
      Double Precision DGRV_NLN
      DOUBLE PRECISION DGRV_NL,BLUEML,d_xgx,XGLU
      Double Precision SIGM0,LAM,X0,R02,alphas,glu
      Integer IGLU
      Common/CAGLUON/Iglu
      Integer ISEL
      Double Precision UPV,DNV,USEA,DSEA
      Common /collpdf/UPV,DNV,USEA,DSEA 
	Double Precision xminpdf,xmaxpdf,q2minpdf,q2maxpdf
	common/W50513/xminpdf,xmaxpdf,q2minpdf,q2maxpdf
      EXTERNAL DGRV_NLN,DGRV_NL,BLUEML,D_XGX,XGLU
      Double Precision d_dnv,d_upv
      External d_dnv,d_upv
      LOGICAL FIRST,PDFFIRST
      Data PDFFIRST/.TRUE./
      Data FIRST/.TRUE./
      Logical firstTMD
      Data firstTMD/.true./
      Integer Iquark
      Data Iquark/0/
C TMDlib variables   
      Integer isel_save   
      Integer IsetTMDlib
      Character*132 name, namet
      Double Precision up,ubar,dn,dbar,strange,sbar,charm,cbar,bottom,bbar,xbar,eps
      Double Precision kt
      Data eps/1.e-4/
      SAVE ISEL  
      Double Precision kt0,x00
c ISEL = 1 CCFM
c      = 2 derivative of collinear gluon
c      = 3 Bluemlein
c      = 4 KMS
c      = 5 Golec-Biernat Wuesthoff
c      = 6 KMR
      IF(FIRST) THEN
	   If(IGLU.GT.1000.and.IGLU.lt.10000) then
	   ISEL = Int(IGLU/1000)
         else
         ISEL = IGLU
	   ENDIF
c	   write(6,*) ' isel ',ISEL,int(Iglu/1000),int(Iglu/10000),iglu
         FIRST=.FALSE.
         Write(6,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
         Write(6,*) '++++  TMD density selection    ++++++++++++++'
         Write(6,*) '++++  ISEL= ',ISEL,
     +      '                  ++++'
         If(Isel.eq.1) then 
	     if(iglu.lt.1000) then
           Write(6,*) '++++  ISEL=1 --> CCFM   JS2001             ++++'
	     elseif(iglu.gt.1000) then
	        Iccfmset =  MOD(IGLU,1000)
     	        write(6,*) '++++  Iccfmset = ',Iccfmset
	        if(iccfmset.eq.1) then
           Write(6,*) '++++  CCFM  J2003 set 1 (Ipgg0ns-1)        ++++'
	        elseif(iccfmset.eq.2) then
           Write(6,*) '++++  CCFM J2003 set 2 (Ipgg1ns1)          ++++'
	        elseif(iccfmset.eq.3) then
           Write(6,*) '++++  CCFM J2003 set 3 (Ipgg2ns2)          ++++'
	        elseif(iccfmset.eq.10) then 
           Write(6,*) '++++  CCFM set A0                          ++++'
	        elseif(iccfmset.eq.11) then
           Write(6,*) '++++  CCFM set A0+                         ++++'
	        elseif(iccfmset.eq.12) then
           Write(6,*) '++++  CCFM set A0-                         ++++'
	        elseif(iccfmset.eq.13) then
           Write(6,*) '++++  CCFM set A1                          ++++'
	        elseif(iccfmset.eq.20) then
           Write(6,*) '++++  CCFM set B0                          ++++'
	        elseif(iccfmset.eq.21) then
           Write(6,*) '++++  CCFM set B0+                         ++++'
	        elseif(iccfmset.eq.22) then
           Write(6,*) '++++  CCFM set B0-                         ++++'
	        elseif(iccfmset.eq.23) then
           Write(6,*) '++++  CCFM set B1                          ++++'
	        elseif(iccfmset.eq.101) then
           Write(6,*) '++++  CCFM set C                           ++++'
	        elseif(iccfmset.eq.111) then
           Write(6,*) '++++  CCFM test set                        ++++'
	        elseif(iccfmset.ge.200.and.iccfmset.le.300) then
           Write(6,*) '++++  CCFM JH 2013 set1                    ++++'
	        elseif(iccfmset.ge.300.and.iccfmset.le.400) then
           Write(6,*) '++++  CCFM JH 2013 set2                    ++++'
	        else
	     Write(6,*) '++++ no valid set specified:',iglu,'++++'
	        endif
	     endif
	     elseif(isel.eq.2) then
         Write(6,*) '++++  ISEL=2 --> derivative of coll. gluon ++++'
	     elseif(isel.eq.3) then
         Write(6,*) '++++  ISEL=3 --> Bluemlein                 ++++'
	     elseif(isel.eq.4) then
         Write(6,*) '++++  ISEL=4 --> Kwiecinski,Martin,Stasto  ++++'
	     elseif(isel.eq.5) then
         Write(6,*) '++++  ISEL=5 --> Wuesthoff Golec-Biernat   ++++'
	     elseif(isel.eq.6) then
         Write(6,*) '++++  ISEL=6 --> Kimber,Martin,Ryskin      ++++'
	     elseif(isel.eq.7) then
         Write(6,*) '++++  ISEL=7 --> Ryskin,Shabelski          ++++'
	   endif
         Write(6,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
      ENDIF

      X = XI
      KT2 = KT2I
c for new p_max
      T=max(0.0001,TI)
      DO 10 I=-6,6
         XPQ(i)=0.d0
   10 XPQS(I)=0.0
      SCALE = KT2
      If(ISEL.EQ.1) THEN
	   if(iabs(kf).eq.2212) then 
            call ccfmstfu(kf,x,scale,t,xpqs)
	   elseif(kf.eq.22.or.iabs(kf).eq.11) then
            call ccfmstfg(kf,x,scale,t,xpqs)
	   else
            write(6,*) ' cauniglu: no valid kf = ',kf	      
	   endif 
         if(iquark.eq.0.and.iccfmset.lt.111) then
            call ccfmstfq(kf,x,scale,t,xpqq)
            if(kf.gt.0) then 
              do I=1,6
                xpqs(i) = xpqq(i)
                xpqs(-i) = xpqq(-i)
              end do 
              else
              do I=1,6
                xpqs(i) = xpqq(-i)
                xpqs(-i) = xpqq(i)
              end do
              endif
         elseif(iquark.eq.0.and.iccfmset.ge.111) then
c           write(6,*)  ' kinematic ', x,scale,t
c           write(6,*)  ' ccfmstu ',xpqs(0),xpqs(1),xpqs(2)
         else
c	   XPQS(1) = dnv/q2minpdf
c	   XPQS(2) = upv/q2minpdf
c	   XPQS(1) = dnv/scale 
c	   XPQS(2) = upv/scale
         SCALE = MAX(q2minpdf,KT2)  
         call xpdf(x,scale)
c         call xpdf(dble(x),dble(t**2))
	   if(kt2.le.q2minpdf) then
  	      XPQS(1) = (dnv)/q2minpdf 
	      XPQS(2) = (upv)/q2minpdf
 	   elseif(kt2.le.100.) then
            XPQS(1) = d_dnv(X,SCALE)
            XPQS(1) = max(0.,xpqs(1))
            XPQS(2) = d_upv(X,SCALE)
            XPQS(2) = max(0.,xpqs(2))            
         endif
c  	      XPQS(1) = (dnv)/scale
c	      XPQS(2) = (dnv)/scale
         endif
c         write(6,*) ' cauniglu kf,x,scale,t ',kf,x,scale,t
c         write(6,*) ' cauniglu ',xpqs(0),xpqs(1),xpqs(2)
      ELSEIF(ISEL.EQ.2) THEN
	   If(pdffirst) then
	      XPQS(0) = XGLU(X,100d0)
		pdffirst=.false.
		endif
         SCALE = MAX(q2minpdf,KT2)
	   if(kt2.le.q2minpdf) then
	      XPQS(0) = XGLU(X,SCALE)/q2minpdf
	   else
c analytical derivative of grv (fast)
cc         SCALE = MAX(0.35,KT2)
cc         XPQS(0) = Real(dgrv_nln(X,SCALE))
c numerical derivative of gluon (slow)
            XPQS(0) = d_xgx(X,SCALE)
	   endif
c old version
cc         write(6,*) ' old version used'
cc         SCALE = MAX(0.35,KT2)
cc         XPQS(0) = d_xgx(X,SCALE)
cend old version	   
         IF(XPQS(0).LE.0.) XPQS(0)=0.
      ELSEIF(ISEL.EQ.3) THEN
	   if(kf.eq.22.or.iabs(kf).eq.11) then
	   write(6,*) ' cauniglu : '
	   write(6,*) ' requested set',isel,' not available for kf = ',kf
         Stop
	   endif 
c select Bluemlein parametrisation
         XPQS(0) = SNGL(BLUEML(X,SCALE,T))
      ELSEIF(ISEL.EQ.4) THEN
	   if(kf.eq.22.or.iabs(kf).eq.11) then
	   write(6,*) ' cauniglu : '
	   write(6,*) ' requested set',isel,' not available for kf = ',kf
         Stop
	   endif 
         SCALE = KT2
         if(kt2.lt.1.) then
            scale = 1.01
            call ktstfu(kf,x,scale,xpqs)
         else
            if(kt2.gt.1000.) scale = 1000.
            call ktstfu(kf,x,scale,xpqs)
         endif
            call ccfmstfq(kf,x,scale,t,xpqq)
            if(kf.gt.0) then 
              do I=1,6
                xpqs(i) = xpqq(i)
                xpqs(-i) = xpqq(-i)
              end do
              else
              do I=1,6
                xpqs(i) = xpqq(-i)
                xpqs(-i) = xpqq(i)
              end do
              endif
      ELSEIF(ISEL.EQ.5) THEN
	   if(kf.eq.22.or.iabs(kf).eq.11) then
	   write(6,*) ' cauniglu : '
	   write(6,*) ' requested set',isel,' not available for kf = ',kf
         Stop
	   endif 
         SCALE = KT2
cccc this is from SATRAP
c ori      SIGM0 = 23.03*2.568
c          mb -> 1/GeV^2 = 1/0.389379 = 2.568
c         SIGM0 = 23.03d0*2.568d0 ! w/o charm
*! with charm
         SIGM0 = 29.12d0*2.568d0
c ori      LAM  = 0.288  ! w/o charm
*! with charm
         LAM = 0.277d0
c ori      X0 =  3.04/10000. ! w/o charm
*! with charm
         X0 = 0.41d0/10000.d0
         ALPHAS = 0.2d0
         R02 = (dble(X)/X0)**LAM
         glu = 3.d0*SIGM0*R02*SCALE*EXP(-R02*SCALE)
c pi**2 = 9.8696
         XPQS(0) = glu/4.d0/9.8696d0/ALPHAS
      ELSEIF(ISEL.EQ.6) THEN
	   if(kf.eq.22.or.iabs(kf).eq.11) then
	   write(6,*) ' cauniglu : '
	   write(6,*) ' requested set',isel,' not available for kf = ',kf
         Stop
	   endif 
         call ccfmstfu(kf,x,scale,t,xpqs)
      ELSEIF(ISEL.EQ.7) THEN
	   if(kf.eq.22.or.iabs(kf).eq.11) then
	   write(6,*) ' cauniglu : '
	   write(6,*) ' requested set',isel,' not available for kf = ',kf
         Stop
	   endif 
         call ryssha(kf,x,scale,xpqs)
      Elseif(ISEL.GE.10000) then
c first call initial set for setting up min and max values
         isetTMDlib=isel 
         xbar=0
c         write(6,*) ' isel ',isel
         if(firstTMD) then 
            write(6,*) ' cauniglu: initialize TMDs'
cccc            iccfmset =1201
cccc            isel_save = isel
ccc            iglu = 1201
ccc            call ccfmstfu(kf,x,scale,t,xpqs)
cccc            isel = isel_save
            call TMDinit(isetTMDlib)  
            firstTMD=.false.
         endif   
         kt = dsqrt(scale)
c         write(6,*) ' before calling TMDpdf '
         call TMDpdf(kf,x,xbar,kt,t,up ,ubar,dn,dbar,strange,sbar,charm,cbar,bottom,bbar,glu)
c           write(6,*)  ' kinematic ', x,kt,t
c         write(6,*)  ' cauniglu: TMDpdf ',kf,x,xbar,kt,t,up ,ubar,dn,dbar,strange,sbar,charm,cbar,bottom,bbar,glu
         xpqs(0) = glu
         xpqs(1) = dn
         xpqs(-1) = dbar
         xpqs(2) = up
         xpqs(-2) = ubar
         xpqs(3) = strange
         xpqs(-3) = sbar
         xpqs(4) = charm
         xpqs(-4) = cbar
         xpqs(5) = bottom 
         xpqs(-5) = bbar 
      ELSE
         write(6,*) ' ISEL = ',ISEL,' not implemented '
         STOP
      ENDIF
c test for saturation
      kt0 = 2.
      x00 = 0.1      
      DO I=-6,6
         XPQ(I)=XPQS(I)
c	 if(kt2.le.kt0) then
c         XPQ(I)=XPQS(I)*kt2**2
c	 endif
ccccccccc        XPQ(I)=XPQS(I)*(scale)/(scale + kt0**2)
ccc         XPQ(I)=XPQS(I)*(x00/(x + x00))      
	 ENDDO
      RETURN
      END
      Subroutine XPDF(XX,QQ2)
      Implicit none
#include "casglu.inc"
      Double Precision PARP,PARI
      INTEGER MSTP,MSTI
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)

      DOUBLE PRECISION XX,QQ2,QQ22,PP2,XPQC(-6:6)
      DOUBLE PRECISION XXX,QQ,STR,CHM,BOT,TOP,GLU,
     +VAL(20)
      Double Precision UPV,DNV,USEA,DSEA
      Common /collpdf/UPV,DNV,USEA,DSEA 
      LOGICAL PDFFIRST
      COMMON /W50516/PDFFIRST
      CHARACTER*20 PARM(20)
      Integer NPDF,IP2
      DATA NPDF/0/,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU/9*0.D0/
      CHARACTER*20 LHAPARM(20)
      DOUBLE PRECISION LHAVALUE(20)
      COMMON/LHACONTROL/LHAPARM,LHAVALUE
      SAVE/LHACONTROL/

      QQ=SQRT(MAX(2.0D0,QQ2))
      XXX = Max(1e-5,XX)
            

ccc for LHAPDF 5.7            
      PARM(1) = 'DEFAULT'
      VAL(1) = ICOLGLU
cccc      write(6,*) ' XPDF', Val(3)
ccc      write(6,*) ' xpdf ',xx,qq2  
CNEW
      NPDF=NPDF+1
      PDFFIRST=.FALSE.
      IF(NPDF.LE.1) THEN
         PDFFIRST=.TRUE.
      ENDIF
c call PDFSET each time
      CALL PDFSET(PARM,VAL)
      CALL STRUCTM(XXX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
      RETURN
      END
      
      
      FUNCTION D_DNV(XX,QQ2)
      Implicit None
      DOUBLE PRECISION XDNV,X,XX,QQ2,DELTA,DXQX,RERR,D_DNV
      DOUBLE PRECISION DQ2
      COMMON/gluon1/X
      EXTERNAL XDNV
      DELTA = 1.D-1*QQ2
      RERR = -9999.D0
      DXQX = -99999.D0
      DQ2 = QQ2
      X = XX
      CALL DFRIDR(XDNV,DQ2,DELTA,DXQX,RERR)
      IF(RERR.EQ.-9999.D0) THEN
         write(6,*) ' D_DNV: error in calculating derivative of xq(x,'
     +   //'q2)'
         write(6,*) ' error = ',RERR
         DXQX = -99999.D0
      ENDIF
      IF(RERR/dxqx.GT.5.) THEN
         write(6,*) ' D_DNV: large error for derivative of xq(x,'
     +   //'q2)'
         write(6,*) ' error = ',RERR/dxqx,' larger than 5. '
         write(6,*) ' with RERR = ',RERR,' and dxqx = ',dxqx
	   write(6,*) ' results are not reliable'
	   write(6,*) ' check selected set of pdf"s'
	   write(6,*) ' are they valid for Q2 = ',DQ2,' and x = ',XX
	   write(6,*) ' stop the program; derivative set to 0.0000 '
         DXQX = 0.D0
      ENDIF
c      write(6,*) ' D_DNV: DFRIDR DXQX,DELTA,RERR ',DXQX,DELTA,RERR
      D_DNV = DXQX
      RETURN
      END
      FUNCTION D_UPV(XX,QQ2)
      Implicit None
      DOUBLE PRECISION XUPV,X,XX,QQ2,DELTA,DXQX,RERR,D_UPV
      DOUBLE PRECISION DQ2
      COMMON/gluon1/X
      EXTERNAL XUPV
      DELTA = 1.D-1*QQ2
      RERR = -9999.D0
      DXQX = -99999.D0
      DQ2 = QQ2
      X = XX
      CALL DFRIDR(XUPV,DQ2,DELTA,DXQX,RERR)
      IF(RERR.EQ.-9999.D0) THEN
         write(6,*) ' D_UPV: error in calculating derivative of xq(x,'
     +   //'q2)'
         write(6,*) ' error = ',RERR
         DXQX = -99999.D0
      ENDIF
      IF(RERR/dxqx.GT.5.) THEN
         write(6,*) ' D_UPV: large error for derivative of xq(x,'
     +   //'q2)'
         write(6,*) ' error = ',RERR/dxqx,' larger than 5. '
         write(6,*) ' with RERR = ',RERR,' and dxqx = ',dxqx
	   write(6,*) ' results are not reliable'
	   write(6,*) ' check selected set of pdf"s'
	   write(6,*) ' are they valid for Q2 = ',DQ2,' and x = ',XX
	   write(6,*) ' stop the program; derivative set to 0.0000 '
         DXQX = 0.D0
      ENDIF
c      write(6,*) ' D_DNV: DFRIDR DXQX,DELTA,RERR ',DXQX,DELTA,RERR
      D_UPV = DXQX
      RETURN
      END
      FUNCTION XUPV(QQ2)
	Implicit none
      DOUBLE PRECISION X,QQ2
      COMMON/gluon1/X
      DOUBLE PRECISION XXX,Q2,XUPV
      Double Precision UPV,DNV,USEA,DSEA
      Common /collpdf/UPV,DNV,USEA,DSEA 
      Q2=QQ2
      XXX = X
      call XPDF(XXX,Q2)
c      XUPV = upv+usea
      XUPV = upv 
      RETURN
      END
      FUNCTION XDNV(QQ2)
	Implicit none
      DOUBLE PRECISION X,QQ2
      COMMON/gluon1/X
      DOUBLE PRECISION XXX,Q2,XDNV
      Double Precision UPV,DNV,USEA,DSEA
      Common /collpdf/UPV,DNV,USEA,DSEA 
      Q2=QQ2
      XXX = X
      call XPDF(XXX,Q2)
c      XDNV = dnv+dsea
      XDNV = dnv 
      RETURN
      END
