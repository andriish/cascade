      SUBROUTINE RYSSHA(KF,X,Q2,XPQ)
      IMPLICIT None
#include "casglu.inc"
      Double Precision XPQ(-6:6),X,Q2
      Integer KF
      DOUBLE PRECISION XXX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU,
     +VAL(20)
      LOGICAL PDFFIRST
      COMMON /W50516/PDFFIRST
      CHARACTER*20 PARM(20)
      Integer MST51,NPDF
      DATA NPDF/0/,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU/9*0.D0/
      Double Precision  phi,Cq,k0x,q02
      Double Precision a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,alpha,beta,eps
      Data a1/163.2/,a2/-0.137/,a3/0.888/,a4/-0.423/
      Data a5/0.000187/,a6/1.0324/,a7/1.1038/,a8/0.000141/
      Data a9/-4.344/,a10/1.3846/,a11/-1.1098/
      Data alpha/11.94/,beta/1.789/,eps/0.00001/
      Double Precision pi
      Data pi/3.1415927/
      Logical First
      Data First/.TRUE./

ccc for LHAPDF 5.7            
      PARM(1) = 'DEFAULT'
      VAL(1) = ICOLGLU
      IF(FIRST) THEN
         CALL PDFSET(PARM,VAL)
         FIRST=.FALSE.
      ENDIF

*! for GRV 1990
      q02=4.
      k0x=2.+0.053*0.053*exp(3.56*sqrt(max(0.,log(1./(3.*x)))))
      if(q2.lt.k0x) then
         Cq=1.
      else
         Cq=k0x/q2
      endif

      phi=a1/(a2+a3+a4)*(a2+a3*q02/q2+a4*(q02/q2)**2 +
     + alpha*x+beta/(eps+log(1./x)))*Cq*sqrt(a5/(a5+x)) *
     + (1.-a6*x**a7 *log(q2/a8))*(1.+a11*x)*(1.-x)**(a9+a10*log(q2/a8))
      if(q2.gt.q02) then
         xpq(0) = phi/4./sqrt(2.)/pi**3
c     write(6,*) ' phi = ',phi
      else
         QQ=sqrt(q02)
         XXX = X
         CALL STRUCTM(XXX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
c     write(6,*) ' glu = ',glu
         xpq(0) = glu/q02
      endif
      return
      END
