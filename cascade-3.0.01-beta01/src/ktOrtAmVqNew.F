      function DOT4(VV,U)
      implicit none
#include "caslujets.inc"
#include "caspart.inc"


      double precision DOT4
      double precision VV(4),U(4)

      DOT4=VV(4)*U(4)-VV(1)*U(1)-VV(2)*U(2)-VV(3)*U(3)

      return
      end

      function DOTP4(I,U)
      implicit none
#include "caslujets.inc"
#include "caspart.inc"


      double precision DOTP4
      double precision U(4)
      integer I

      DOTP4=P(I,4)*U(4)-P(I,1)*U(1)-P(I,2)*U(2)-P(I,3)*U(3)

      return
      end

      function DOTE(I,J)
      implicit none
#include "caslujets.inc"
#include "caspart.inc"


      double precision DOTE
      integer I,J

      DOTE=P(I,1)*P(J,1)+P(I,2)*P(J,2)

      return
      end

      function ME(NI1,NI2,N1,N2,kpz,ud)
      implicit none

#include "caslujets.inc"
#include "caspart.inc"
      DOUBLE PRECISION SIN2W,XMW2
      COMMON/ELWEAK/SIN2W,XMW2
	DOUBLE PRECISION XMW


      integer NI1,NI2,N1,N2,a,b,c,i,j,kk,vb
      double precision ME
      double precision sh,s1h,s2h,t1h,t2h,u1h,u2h,t1,t2
      double precision m1,m2,mw,m1t,m2t,mwt,mqt,s
      double precision er1,er2,xw,dot1,dot,dote,dotp4,dot4
      double precision q1tk1t,q1tk2t,q1tkwt,q2tk1t,q2tk2t,q2tkwt
      double precision q1tq2t,k1tk2t,k1tkwt,k2tkwt,vq,aq
      double precision P1,P2,P3,P4,P5,P6,P7,P8,q1k1,q2k1,aa2,bb2,cc1,
     &ac,bc,i1a,i2b,i1b,i2a,i1c,i2c,iin,aa,bb,cc,ab,cc2,
     &k1k2,k1kw,k2kw,q1k2,q2k2,q1kw,q2kw,rw,q1t(4),q2t(4)
ccccccccccccccccccccc
c first index projection, second index component
ccccccccccccccccccccc
      double precision L(4),U(4),A8(4,4),AS1(4,4),AS2(4,4),
     & B8(4),BS1(4),BS2(4),B58(4),B51(4),B52(4)
      double precision A58(4,4),A51(4,4),A52(4,4),C1(4),C2(4),
     & D1(4),D2(4)
      double precision A1,A2,B1,B2,A15,A25,B15,B25,UN,LN,qn1,qn2
ccccccccccccccccccccc
      double precision pAe1,pBe1,k1e1,k2e1,kwe1,q1e1,q2e1,e1e2,
     & pAe2,pBe2,k1e2,k2e2,kwe2,q1e2,q2e2,Le1,Ue1,Le2,Ue2,e1(4),e2(4)
      double precision N0(4),N5(4),q1U,q2U,q1L,q2L,pAL,q2tL,q2tU,
     &                 q1tL,q1tU,w,q1q2,kwpA,k1pA,k2pA,k1pB,k2pB,
     &                 kwpB,q1pB,q2pA,
     &                 pAU,pBL,pBU,M(4,4,4),M5(4,4,4),d(4,4),dd(4,4)
      double precision q1(4),q2(4),k1(4),k2(4),kw(4),pA(4),pB(4)
      
      Integer kpz,ud

      a=N1
      b=N2
      c=N1+1

      m1=sqrt(DOT1(a,a))
      m2=sqrt(DOT1(b,b))
      mw=sqrt(DOT1(c,c))

c      m1=0
c      m2=0

      s=SSS
      w=sqrt(s)/2.
      t1=DOT1(NI1,NI1)
      qn1=sqrt(dabs(t1))
      t2=DOT1(NI2,NI2)
      qn2=sqrt(dabs(t2))
      sh=t1+t2+2*DOT1(NI1,NI2)
      s1h=m1**2+mw**2+2*DOT1(a,c)
      s2h=m2**2+mw**2+2*DOT1(b,c)
      t1h=m1**2+t1-2*DOT1(NI1,a)
      t2h=m2**2+t2-2*DOT1(NI2,b)
      u1h=t1+m2**2-2*DOT1(NI1,b)
      u2h=t2+m1**2-2*DOT1(NI2,a)
      m1t=sqrt(m1**2+DOTE(a,a))
      m2t=sqrt(m2**2+DOTE(b,b))
      mwt=sqrt(mw**2+DOTE(c,c))
      mqt=sqrt(sh-t1-t2+2*DOTE(NI1,NI2))
c      er1=4*aa2*cc1/(mwt*m1t*s)
c      er2=4*bb2*cc1/(mwt*m2t*s)
      rw=DOT1(c,2)/mwt/w
      q1tk1t=DOTE(NI1,a)
      q1tk2t=DOTE(NI1,b)
      q1tkwt=DOTE(NI1,c)
      q2tk1t=DOTE(NI2,a)
      q2tk2t=DOTE(NI2,b)
      q2tkwt=DoTE(NI2,c)
      q1tq2t=DOTE(NI1,NI2)
      k1tk2t=DOTE(a,b)
      k1tkwt=DOTE(a,c)
      k2tkwt=DOTE(b,c)

      q1q2=DOT1(NI1,NI2)
      q1pB=DOT1(NI1,2)
      q2pA=DOT1(NI2,1)
      k1pA=DOT1(a,1)
      k2pA=DOT1(b,1)
      kwpA=DOT1(c,1)
      k1pB=DOT1(a,2)
      k2pB=DOT1(b,2)
      kwpB=DOT1(c,2)
      k1k2=DOT1(a,b)
      k1kw=DOT1(a,c)
      k2kw=DOT1(b,c)
      q1k1=DOT1(a,NI1)
      q2k1=DOT1(a,NI2)
      q1k2=DOT1(b,NI1)
      q2k2=DOT1(b,NI2)
      q1kw=DOT1(c,NI1)
      q2kw=DOT1(c,NI2)

      P1=1./(sh*(s2h-m1**2))
      P2=1./(sh*(s1h-m2**2))
      P3=1./((t1h-m1**2)*(s2h-m1**2))
      P4=1./((s1h-m2**2)*(t2h-m2**2))
      P5=1./((s2h-m1**2)*(u2h-m1**2))
      P6=1./((s1h-m2**2)*(u1h-m2**2))
      P7=1./((t2h-m2**2)*(t1h-m1**2))
      P8=1./((u1h-m2**2)*(u2h-m1**2))

      do i=1,4
        q1(i)=P(NI1,i)
        q2(i)=P(NI2,i)
        k1(i)=P(a,i)
        k2(i)=P(b,i)
        kw(i)=P(c,i)
        pA(i)=P(1,i)
        pB(i)=P(2,i)
      enddo

c       write(6,*) '*'
c       write(6,*) q1
c       write(6,*) pA
c       write(6,*) q2
c       write(6,*) pB

      q1t(1)=q1(1)
      q1t(2)=q1(2)
      q1t(3)=0.
      q1t(4)=0.

      q2t(1)=q2(1)
      q2t(2)=q2(2)
      q2t(3)=0.
      q2t(4)=0.

      e1(1)=q1(1)/sqrt(q1(1)**2+q1(2)**2)
      e1(2)=q1(2)/sqrt(q1(1)**2+q1(2)**2)
      e1(3)=0.
      e1(4)=0.

      e2(1)=q2(1)/sqrt(q2(1)**2+q2(2)**2)
      e2(2)=q2(2)/sqrt(q2(1)**2+q2(2)**2)
      e2(3)=0.
      e2(4)=0.

      U(1)=(P(a,2)*P(b,3)*P(c,4)+P(a,3)*P(b,4)*P(c,2)
     &    +P(a,4)*P(b,2)*P(c,3)-P(a,2)*P(b,4)*P(c,3)
     &    -P(a,3)*P(b,2)*P(c,4)-P(a,4)*P(b,3)*P(c,2))
     &/sqrt(p(a,4)*p(b,4)*p(c,4))
      U(2)=-(P(a,1)*P(b,3)*P(c,4)+P(a,3)*P(b,4)*P(c,1)
     &    +P(a,4)*P(b,1)*P(c,3)-P(a,1)*P(b,4)*P(c,3)
     &    -P(a,3)*P(b,1)*P(c,4)-P(a,4)*P(b,3)*P(c,1))
     &/sqrt(p(a,4)*p(b,4)*p(c,4))
      U(3)=-(P(a,2)*P(b,1)*P(c,4)+P(a,1)*P(b,4)*P(c,2)
     &    +P(a,4)*P(b,2)*P(c,1)-P(a,2)*P(b,4)*P(c,1)
     &    -P(a,1)*P(b,2)*P(c,4)-P(a,4)*P(b,1)*P(c,2))
     &/sqrt(p(a,4)*p(b,4)*p(c,4))
      U(4)=-(P(a,1)*P(b,3)*P(c,2)+P(a,2)*P(b,1)*P(c,3)
     &    +P(a,3)*P(b,2)*P(c,1)-P(a,1)*P(b,2)*P(c,3)
     &    -P(a,3)*P(b,1)*P(c,2)-P(a,2)*P(b,3)*P(c,1))
     &/sqrt(p(a,4)*p(b,4)*p(c,4))

      L(1)=(P(a,2)*P(b,3)*U(4)+P(a,3)*P(b,4)*U(2)
     &    +P(a,4)*P(b,2)*U(3)-P(a,2)*P(b,4)*U(3)
     &    -P(a,3)*P(b,2)*U(4)-P(a,4)*P(b,3)*U(2))
     &/(p(a,4)*p(b,4)*abs(U(4)))
      L(2)=-(P(a,1)*P(b,3)*U(4)+P(a,3)*P(b,4)*U(1)
     &    +P(a,4)*P(b,1)*U(3)-P(a,1)*P(b,4)*U(3)
     &    -P(a,3)*P(b,1)*U(4)-P(a,4)*P(b,3)*U(1))
     &/(p(a,4)*p(b,4)*abs(U(4)))
      L(3)=-(P(a,2)*P(b,1)*U(4)+P(a,1)*P(b,4)*U(2)
     &    +P(a,4)*P(b,2)*U(1)-P(a,2)*P(b,4)*U(1)
     &    -P(a,1)*P(b,2)*U(4)-P(a,4)*P(b,1)*U(2))
     &/(p(a,4)*p(b,4)*abs(U(4)))
      L(4)=-(P(a,1)*P(b,3)*U(2)+P(a,2)*P(b,1)*U(3)
     &    +P(a,3)*P(b,2)*U(1)-P(a,1)*P(b,2)*U(3)
     &    -P(a,3)*P(b,1)*U(2)-P(a,2)*P(b,3)*U(1))
     &/(p(a,4)*p(b,4)*abs(U(4)))

c      write(6,*) '1.',DOT4(U,U),DOT4(L,L)

      UN=DOT4(U,U)
      LN=DOT4(L,L)

      do i=1,4
        L(i)=L(i)/sqrt(-LN)
        U(i)=U(i)/sqrt(-UN)
      enddo

      ab=DOT1(a,b)

      UN=DOT4(U,U)
      LN=DOT4(L,L)

c      write(6,*) UN,LN

      N0(1)=-m1*m2+4*ab
      N0(2)=-4*ab*UN-4*m1*m2*UN
      N0(3)=-4*ab*LN-4*m1*m2*LN
      N0(4)=-4*ab*UN*LN+4*m1*m2*UN*LN

      N5(1)=-m1*m2-4*ab
      N5(2)=-4*ab*UN+4*m1*m2*UN
      N5(3)=-4*ab*LN+4*m1*m2*LN
      N5(4)=4*ab*UN*LN+4*m1*m2*UN*LN

      q1U=DOTP4(NI1,U)
      q2U=DOTP4(NI2,U)
      q1tU=-(P(NI1,1)*U(1)+P(NI1,2)*U(2))
      q2tU=-(P(NI2,1)*U(1)+P(NI2,2)*U(2))
      q1tL=-(P(NI1,1)*L(1)+P(NI1,2)*L(2))
      q2tL=-(P(NI2,1)*L(1)+P(NI2,2)*L(2))
      q1L=DOTP4(NI1,L)
      q2L=DOTP4(NI2,L)
      pAU=DOTP4(1,U)
      pBU=DOTP4(2,U)
      pAL=DOTP4(1,L)
      pBL=DOTP4(2,L)

      e1e2=DOT4(e1,e2)
      pAe1=DOTP4(1,e1)
      pBe1=DOTP4(2,e1)
      k1e1=DOTP4(a,e1)
      k2e1=DOTP4(b,e1)
      kwe1=DOTP4(c,e1)
      q1e1=DOTP4(NI1,e1)
      q2e1=DOTP4(NI2,e1)
      Ue1=DOT4(U,e1)
      Le1=DOT4(L,e1)

      pAe2=DOTP4(1,e2)
      pBe2=DOTP4(2,e2)
      k1e2=DOTP4(a,e2)
      k2e2=DOTP4(b,e2)
      kwe2=DOTP4(c,e2)
      q1e2=DOTP4(NI1,e2)
      q2e2=DOTP4(NI2,e2)
      Ue2=DOT4(U,e2)
      Le2=DOT4(L,e2)

      
      if(iabs(kpz).eq.24) then
cccc W
         vq=1./(2.*sqrt(2.))
         aq=1./(2.*sqrt(2.))  
      elseif(iabs(kpz).eq.23.and.ud.eq.0) then
cccc Z u-type
         vq=(1./2.-4./3.*SIN2W)/(2*sqrt(1.-SIN2W))
         aq=1./2./(2*sqrt(1.-SIN2W))
      elseif(iabs(kpz).eq.23.and.ud.eq.1) then
cccc Z d-type
         vq=(-1./2.+2./3.*SIN2W)/(2*sqrt(1.-SIN2W))
         aq=-1./2./(2*sqrt(1.-SIN2W))
      else
         Write(6,*) ' ktOrtAmVqNew: invalid combination of z and ud ',
     &     kpz,ud
         write(6,*) ' program stopped '
         stop
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccc

      do vb=1,4

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                kt1kt2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      A8(1,vb) =
     &  - 8*k1(vb)*q1tq2t*q1tk1t*P1*m2 - 8*k1(vb)*q1tq2t*q1tk1t*P2*
     & m2 - 8*k1(vb)*q1tq2t*q1tk2t*P1*m1 + 8*k1(vb)*q1tq2t*q1tk2t*
     & P2*m1 + 16*k1(vb)*q1tq2t*q1tk2t*P2*m2 + 8*k1(vb)*q1tq2t*
     & q2tk1t*P1*m2 + 8*k1(vb)*q1tq2t*q2tk1t*P2*m2 + 8*k1(vb)*
     & q1tq2t*q2tk2t*P1*m1 - 8*k1(vb)*q1tq2t*q2tk2t*P2*m1 - 16*
     & k1(vb)*q1tq2t*q2tk2t*P2*m2 - 4*k1(vb)*q1tq2t*P1*m1*q1k2 + 4*
     & k1(vb)*q1tq2t*P1*m1*q2k2 - 4*k1(vb)*q1tq2t*P1*m2*q1k1 + 4*
     & k1(vb)*q1tq2t*P1*m2*q2k1 + 4*k1(vb)*q1tq2t*P2*m1*q1k2 - 4*
     & k1(vb)*q1tq2t*P2*m1*q2k2 - 4*k1(vb)*q1tq2t*P2*m2*q1k1 + 8*
     & k1(vb)*q1tq2t*P2*m2*q1k2 + 4*k1(vb)*q1tq2t*P2*m2*q2k1 - 8*
     & k1(vb)*q1tq2t*P2*m2*q2k2 + 4*k1(vb)*q1tk1t*t2*P1*m2 + 4*k1(vb)
     & *q1tk1t*t2*P2*m2 + 4*k1(vb)*q1tk2t*t2*P1*m1 - 4*k1(vb)*q1tk2t
     & *t2*P2*m1 - 8*k1(vb)*q1tk2t*t2*P2*m2 - 4*k1(vb)*q2tk1t*t1*P1*
     & m2 - 4*k1(vb)*q2tk1t*t1*P2*m2 - 4*k1(vb)*q2tk2t*t1*P1*m1 + 4*
     & k1(vb)*q2tk2t*t1*P2*m1
      A8(1,vb) = A8(1,vb) + 8*k1(vb)*q2tk2t*t1*P2*m2 + 16*k2(vb)*
     & q1tq2t*q1tk1t*P1*m1 + 8*k2(vb)*q1tq2t*q1tk1t*P1*m2 - 8*
     & k2(vb)*q1tq2t*q1tk1t*P2*m2 - 8*k2(vb)*q1tq2t*q1tk2t*P1*m1 - 
     & 8*k2(vb)*q1tq2t*q1tk2t*P2*m1 - 16*k2(vb)*q1tq2t*q2tk1t*P1*m1
     &  - 8*k2(vb)*q1tq2t*q2tk1t*P1*m2 + 8*k2(vb)*q1tq2t*q2tk1t*P2*
     & m2 + 8*k2(vb)*q1tq2t*q2tk2t*P1*m1 + 8*k2(vb)*q1tq2t*q2tk2t*
     & P2*m1 + 8*k2(vb)*q1tq2t*P1*m1*q1k1 - 4*k2(vb)*q1tq2t*P1*m1*
     & q1k2 - 8*k2(vb)*q1tq2t*P1*m1*q2k1 + 4*k2(vb)*q1tq2t*P1*m1*q2k2
     &  + 4*k2(vb)*q1tq2t*P1*m2*q1k1 - 4*k2(vb)*q1tq2t*P1*m2*q2k1 - 4
     & *k2(vb)*q1tq2t*P2*m1*q1k2 + 4*k2(vb)*q1tq2t*P2*m1*q2k2 - 4*
     & k2(vb)*q1tq2t*P2*m2*q1k1 + 4*k2(vb)*q1tq2t*P2*m2*q2k1 - 8*
     & k2(vb)*q1tk1t*t2*P1*m1 - 4*k2(vb)*q1tk1t*t2*P1*m2 + 4*k2(vb)*
     & q1tk1t*t2*P2*m2 + 4*k2(vb)*q1tk2t*t2*P1*m1 + 4*k2(vb)*q1tk2t*
     & t2*P2*m1 + 8*k2(vb)*q2tk1t*t1*P1*m1 + 4*k2(vb)*q2tk1t*t1*P1*m2
     &  - 4*k2(vb)*q2tk1t*t1*P2*m2 - 4*k2(vb)*q2tk2t*t1*P1*m1 - 4*
     & k2(vb)*q2tk2t*t1*P2*m1
      A8(1,vb) = A8(1,vb) + 8*kw(vb)*q1tq2t*q1tk1t*P1*m2 - 8*kw(vb)
     & *q1tq2t*q1tk1t*P2*m2 - 8*kw(vb)*q1tq2t*q1tk2t*P1*m1 + 8*
     & kw(vb)*q1tq2t*q1tk2t*P2*m1 - 8*kw(vb)*q1tq2t*q2tk1t*P1*m2 + 
     & 8*kw(vb)*q1tq2t*q2tk1t*P2*m2 + 8*kw(vb)*q1tq2t*q2tk2t*P1*m1
     &  - 8*kw(vb)*q1tq2t*q2tk2t*P2*m1 - 4*kw(vb)*q1tq2t*P1*m1*q1k2
     &  + 4*kw(vb)*q1tq2t*P1*m1*q2k2 + 4*kw(vb)*q1tq2t*P1*m2*q1k1 - 4
     & *kw(vb)*q1tq2t*P1*m2*q2k1 + 4*kw(vb)*q1tq2t*P2*m1*q1k2 - 4*
     & kw(vb)*q1tq2t*P2*m1*q2k2 - 4*kw(vb)*q1tq2t*P2*m2*q1k1 + 4*
     & kw(vb)*q1tq2t*P2*m2*q2k1 - 4*kw(vb)*q1tk1t*t2*P1*m2 + 4*kw(vb)
     & *q1tk1t*t2*P2*m2 + 4*kw(vb)*q1tk2t*t2*P1*m1 - 4*kw(vb)*q1tk2t
     & *t2*P2*m1 + 4*kw(vb)*q2tk1t*t1*P1*m2 - 4*kw(vb)*q2tk1t*t1*P2*
     & m2 - 4*kw(vb)*q2tk2t*t1*P1*m1 + 4*kw(vb)*q2tk2t*t1*P2*m1 + 4*
     & q1(vb)*q1tq2t*P1*m1*q1k2 + 4*q1(vb)*q1tq2t*P1*m1*q2k2 - 4*
     & q1(vb)*q1tq2t*P1*m2*q1k1 - 4*q1(vb)*q1tq2t*P1*m2*q2k1 - 4*
     & q1(vb)*q1tq2t*P2*m1*q1k2 - 4*q1(vb)*q1tq2t*P2*m1*q2k2 + 4*
     & q1(vb)*q1tq2t*P2*m2*q1k1
      A8(1,vb) = A8(1,vb) + 4*q1(vb)*q1tq2t*P2*m2*q2k1 - 4*q2(vb)*
     & q1tq2t*P1*m1*q1k2 - 4*q2(vb)*q1tq2t*P1*m1*q2k2 + 4*q2(vb)*
     & q1tq2t*P1*m2*q1k1 + 4*q2(vb)*q1tq2t*P1*m2*q2k1 + 4*q2(vb)*
     & q1tq2t*P2*m1*q1k2 + 4*q2(vb)*q1tq2t*P2*m1*q2k2 - 4*q2(vb)*
     & q1tq2t*P2*m2*q1k1 - 4*q2(vb)*q1tq2t*P2*m2*q2k1 + 4*q1t(vb)*t2*
     & P1*m1*q1k2 + 4*q1t(vb)*t2*P1*m1*q2k2 - 4*q1t(vb)*t2*P1*m2*q1k1
     &  - 4*q1t(vb)*t2*P1*m2*q2k1 - 4*q1t(vb)*t2*P2*m1*q1k2 - 4*q1t(vb)
     & *t2*P2*m1*q2k2 + 4*q1t(vb)*t2*P2*m2*q1k1 + 4*q1t(vb)*t2*P2*m2*
     & q2k1 - 8*q1t(vb)*q1tq2t*P1*m1*q1k2 - 8*q1t(vb)*q1tq2t*P1*m1*
     & q2k2 + 8*q1t(vb)*q1tq2t*P1*m2*q1k1 + 8*q1t(vb)*q1tq2t*P1*m2*
     & q2k1 + 8*q1t(vb)*q1tq2t*P2*m1*q1k2 + 8*q1t(vb)*q1tq2t*P2*m1*
     & q2k2 - 8*q1t(vb)*q1tq2t*P2*m2*q1k1 - 8*q1t(vb)*q1tq2t*P2*m2*
     & q2k1 - 4*q2t(vb)*t1*P1*m1*q1k2 - 4*q2t(vb)*t1*P1*m1*q2k2 + 4*
     & q2t(vb)*t1*P1*m2*q1k1 + 4*q2t(vb)*t1*P1*m2*q2k1 + 4*q2t(vb)*t1*
     & P2*m1*q1k2 + 4*q2t(vb)*t1*P2*m1*q2k2 - 4*q2t(vb)*t1*P2*m2*q1k1
     &  - 4*q2t(vb)*t1*P2*m2*q2k1
      A8(1,vb) = A8(1,vb) + 8*q2t(vb)*q1tq2t*P1*m1*q1k2 + 8*q2t(vb)*
     & q1tq2t*P1*m1*q2k2 - 8*q2t(vb)*q1tq2t*P1*m2*q1k1 - 8*q2t(vb)*
     & q1tq2t*P1*m2*q2k1 - 8*q2t(vb)*q1tq2t*P2*m1*q1k2 - 8*q2t(vb)*
     & q1tq2t*P2*m1*q2k2 + 8*q2t(vb)*q1tq2t*P2*m2*q1k1 + 8*q2t(vb)*
     & q1tq2t*P2*m2*q2k1

      A8(2,vb) =
     & 8*k1(vb)*q2U*q1tq2t*P1*m1*m2 - 8*k1(vb)*q2U*q1tq2t*P1*q1k2
     &  - 8*k1(vb)*q2U*q1tq2t*P1*q2k2 + 8*k1(vb)*q2U*q1tq2t*P1*k1k2
     &  - 8*k1(vb)*q2U*q1tq2t*P2*m1*m2 - 8*k1(vb)*q2U*q1tq2t*P2*
     & q1k2 - 8*k1(vb)*q2U*q1tq2t*P2*q2k2 - 8*k1(vb)*q2U*q1tq2t*P2*
     & k1k2 + 8*k1(vb)*q1tq2t*q1tU*P1*m1*m2 - 8*k1(vb)*q1tq2t*q1tU*
     & P1*q1k2 - 8*k1(vb)*q1tq2t*q1tU*P1*q2k2 + 8*k1(vb)*q1tq2t*
     & q1tU*P1*k1k2 - 8*k1(vb)*q1tq2t*q1tU*P2*m1*m2 - 8*k1(vb)*
     & q1tq2t*q1tU*P2*q1k2 - 8*k1(vb)*q1tq2t*q1tU*P2*q2k2 - 8*
     & k1(vb)*q1tq2t*q1tU*P2*k1k2 - 8*k1(vb)*q1tq2t*q2tU*P1*m1*m2
     &  + 8*k1(vb)*q1tq2t*q2tU*P1*q1k2 + 8*k1(vb)*q1tq2t*q2tU*P1*
     & q2k2 - 8*k1(vb)*q1tq2t*q2tU*P1*k1k2 + 8*k1(vb)*q1tq2t*q2tU*
     & P2*m1*m2 + 8*k1(vb)*q1tq2t*q2tU*P2*q1k2 + 8*k1(vb)*q1tq2t*
     & q2tU*P2*q2k2 + 8*k1(vb)*q1tq2t*q2tU*P2*k1k2 - 4*k1(vb)*q1tU*
     & t2*P1*m1*m2 + 4*k1(vb)*q1tU*t2*P1*q1k2 + 4*k1(vb)*q1tU*t2*P1*
     & q2k2
      A8(2,vb) = A8(2,vb) - 4*k1(vb)*q1tU*t2*P1*k1k2 + 4*k1(vb)*
     & q1tU*t2*P2*m1*m2 + 4*k1(vb)*q1tU*t2*P2*q1k2 + 4*k1(vb)*q1tU*
     & t2*P2*q2k2 + 4*k1(vb)*q1tU*t2*P2*k1k2 + 4*k1(vb)*q2tU*t1*P1*m1
     & *m2 - 4*k1(vb)*q2tU*t1*P1*q1k2 - 4*k1(vb)*q2tU*t1*P1*q2k2 + 4*
     & k1(vb)*q2tU*t1*P1*k1k2 - 4*k1(vb)*q2tU*t1*P2*m1*m2 - 4*k1(vb)*
     & q2tU*t1*P2*q1k2 - 4*k1(vb)*q2tU*t1*P2*q2k2 - 4*k1(vb)*q2tU*t1
     & *P2*k1k2 + 8*k2(vb)*q2U*q1tq2t*P1*m1*m2 + 8*k2(vb)*q2U*
     & q1tq2t*P1*q1k1 + 8*k2(vb)*q2U*q1tq2t*P1*q2k1 + 8*k2(vb)*q2U*
     & q1tq2t*P1*k1k2 - 8*k2(vb)*q2U*q1tq2t*P2*m1*m2 + 8*k2(vb)*q2U
     & *q1tq2t*P2*q1k1 + 8*k2(vb)*q2U*q1tq2t*P2*q2k1 - 8*k2(vb)*q2U
     & *q1tq2t*P2*k1k2 + 8*k2(vb)*q1tq2t*q1tU*P1*m1*m2 + 8*k2(vb)*
     & q1tq2t*q1tU*P1*q1k1 + 8*k2(vb)*q1tq2t*q1tU*P1*q2k1 + 8*
     & k2(vb)*q1tq2t*q1tU*P1*k1k2 - 8*k2(vb)*q1tq2t*q1tU*P2*m1*m2
     &  + 8*k2(vb)*q1tq2t*q1tU*P2*q1k1 + 8*k2(vb)*q1tq2t*q1tU*P2*
     & q2k1 - 8*k2(vb)*q1tq2t*q1tU*P2*k1k2 - 8*k2(vb)*q1tq2t*q2tU*
     & P1*m1*m2
      A8(2,vb) = A8(2,vb) - 8*k2(vb)*q1tq2t*q2tU*P1*q1k1 - 8*k2(vb)
     & *q1tq2t*q2tU*P1*q2k1 - 8*k2(vb)*q1tq2t*q2tU*P1*k1k2 + 8*
     & k2(vb)*q1tq2t*q2tU*P2*m1*m2 - 8*k2(vb)*q1tq2t*q2tU*P2*q1k1
     &  - 8*k2(vb)*q1tq2t*q2tU*P2*q2k1 + 8*k2(vb)*q1tq2t*q2tU*P2*
     & k1k2 - 4*k2(vb)*q1tU*t2*P1*m1*m2 - 4*k2(vb)*q1tU*t2*P1*q1k1 - 
     & 4*k2(vb)*q1tU*t2*P1*q2k1 - 4*k2(vb)*q1tU*t2*P1*k1k2 + 4*k2(vb)
     & *q1tU*t2*P2*m1*m2 - 4*k2(vb)*q1tU*t2*P2*q1k1 - 4*k2(vb)*q1tU*
     & t2*P2*q2k1 + 4*k2(vb)*q1tU*t2*P2*k1k2 + 4*k2(vb)*q2tU*t1*P1*m1
     & *m2 + 4*k2(vb)*q2tU*t1*P1*q1k1 + 4*k2(vb)*q2tU*t1*P1*q2k1 + 4*
     & k2(vb)*q2tU*t1*P1*k1k2 - 4*k2(vb)*q2tU*t1*P2*m1*m2 + 4*k2(vb)*
     & q2tU*t1*P2*q1k1 + 4*k2(vb)*q2tU*t1*P2*q2k1 - 4*k2(vb)*q2tU*t1
     & *P2*k1k2 + 8*kw(vb)*q2U*q1tq2t*P1*m1*m2 + 8*kw(vb)*q2U*
     & q1tq2t*P1*k1k2 - 8*kw(vb)*q2U*q1tq2t*P2*m1*m2 - 8*kw(vb)*q2U
     & *q1tq2t*P2*k1k2 + 8*kw(vb)*q1tq2t*q1tU*P1*m1*m2 + 8*kw(vb)*
     & q1tq2t*q1tU*P1*k1k2 - 8*kw(vb)*q1tq2t*q1tU*P2*m1*m2 - 8*
     & kw(vb)*q1tq2t*q1tU*P2*k1k2
      A8(2,vb) = A8(2,vb) - 8*kw(vb)*q1tq2t*q2tU*P1*m1*m2 - 8*
     & kw(vb)*q1tq2t*q2tU*P1*k1k2 + 8*kw(vb)*q1tq2t*q2tU*P2*m1*m2
     &  + 8*kw(vb)*q1tq2t*q2tU*P2*k1k2 - 4*kw(vb)*q1tU*t2*P1*m1*m2
     &  - 4*kw(vb)*q1tU*t2*P1*k1k2 + 4*kw(vb)*q1tU*t2*P2*m1*m2 + 4*
     & kw(vb)*q1tU*t2*P2*k1k2 + 4*kw(vb)*q2tU*t1*P1*m1*m2 + 4*kw(vb)*
     & q2tU*t1*P1*k1k2 - 4*kw(vb)*q2tU*t1*P2*m1*m2 - 4*kw(vb)*q2tU*
     & t1*P2*k1k2 + 16*U(vb)*q1tq2t*q1tk1t*P1*m1*m2 - 8*U(vb)*q1tq2t
     & *q1tk1t*P1*q1k2 - 8*U(vb)*q1tq2t*q1tk1t*P1*q2k2 + 16*U(vb)*
     & q1tq2t*q1tk1t*P1*k1k2 - 8*U(vb)*q1tq2t*q1tk1t*P2*q1k2 - 8*
     & U(vb)*q1tq2t*q1tk1t*P2*q2k2 + 8*U(vb)*q1tq2t*q1tk2t*P1*q1k1
     &  + 8*U(vb)*q1tq2t*q1tk2t*P1*q2k1 - 16*U(vb)*q1tq2t*q1tk2t*P2
     & *m1*m2 + 8*U(vb)*q1tq2t*q1tk2t*P2*q1k1 + 8*U(vb)*q1tq2t*
     & q1tk2t*P2*q2k1 - 16*U(vb)*q1tq2t*q1tk2t*P2*k1k2 - 16*U(vb)*
     & q1tq2t*q2tk1t*P1*m1*m2 + 8*U(vb)*q1tq2t*q2tk1t*P1*q1k2 + 8*
     & U(vb)*q1tq2t*q2tk1t*P1*q2k2 - 16*U(vb)*q1tq2t*q2tk1t*P1*k1k2
     &  + 8*U(vb)*q1tq2t*q2tk1t*P2*q1k2
      A8(2,vb) = A8(2,vb) + 8*U(vb)*q1tq2t*q2tk1t*P2*q2k2 - 8*U(vb)
     & *q1tq2t*q2tk2t*P1*q1k1 - 8*U(vb)*q1tq2t*q2tk2t*P1*q2k1 + 16*
     & U(vb)*q1tq2t*q2tk2t*P2*m1*m2 - 8*U(vb)*q1tq2t*q2tk2t*P2*q1k1
     &  - 8*U(vb)*q1tq2t*q2tk2t*P2*q2k1 + 16*U(vb)*q1tq2t*q2tk2t*P2
     & *k1k2 + 8*U(vb)*q1tq2t*P1*m1*m2*q1k1 - 8*U(vb)*q1tq2t*P1*m1*m2
     & *q2k1 - 8*U(vb)*q1tq2t*P1*q1k1*q2k2 + 8*U(vb)*q1tq2t*P1*q1k1*
     & k1k2 + 8*U(vb)*q1tq2t*P1*q1k2*q2k1 - 8*U(vb)*q1tq2t*P1*q2k1*
     & k1k2 - 8*U(vb)*q1tq2t*P2*m1*m2*q1k2 + 8*U(vb)*q1tq2t*P2*m1*m2*
     & q2k2 - 8*U(vb)*q1tq2t*P2*q1k1*q2k2 + 8*U(vb)*q1tq2t*P2*q1k2*
     & q2k1 - 8*U(vb)*q1tq2t*P2*q1k2*k1k2 + 8*U(vb)*q1tq2t*P2*q2k2*
     & k1k2 - 8*U(vb)*q1tk1t*t2*P1*m1*m2 + 4*U(vb)*q1tk1t*t2*P1*q1k2
     &  + 4*U(vb)*q1tk1t*t2*P1*q2k2 - 8*U(vb)*q1tk1t*t2*P1*k1k2 + 4*
     & U(vb)*q1tk1t*t2*P2*q1k2 + 4*U(vb)*q1tk1t*t2*P2*q2k2 - 4*U(vb)*
     & q1tk2t*t2*P1*q1k1 - 4*U(vb)*q1tk2t*t2*P1*q2k1 + 8*U(vb)*
     & q1tk2t*t2*P2*m1*m2 - 4*U(vb)*q1tk2t*t2*P2*q1k1 - 4*U(vb)*
     & q1tk2t*t2*P2*q2k1
      A8(2,vb) = A8(2,vb) + 8*U(vb)*q1tk2t*t2*P2*k1k2 + 8*U(vb)*
     & q2tk1t*t1*P1*m1*m2 - 4*U(vb)*q2tk1t*t1*P1*q1k2 - 4*U(vb)*
     & q2tk1t*t1*P1*q2k2 + 8*U(vb)*q2tk1t*t1*P1*k1k2 - 4*U(vb)*
     & q2tk1t*t1*P2*q1k2 - 4*U(vb)*q2tk1t*t1*P2*q2k2 + 4*U(vb)*
     & q2tk2t*t1*P1*q1k1 + 4*U(vb)*q2tk2t*t1*P1*q2k1 - 8*U(vb)*
     & q2tk2t*t1*P2*m1*m2 + 4*U(vb)*q2tk2t*t1*P2*q1k1 + 4*U(vb)*
     & q2tk2t*t1*P2*q2k1 - 8*U(vb)*q2tk2t*t1*P2*k1k2

      A8(3,vb) =
     &  - 8*k1(vb)*q1L*q1tq2t*q1tk2t*P1 - 8*k1(vb)*q1L*q1tq2t*
     & q1tk2t*P2 + 8*k1(vb)*q1L*q1tq2t*q2tk2t*P1 + 8*k1(vb)*q1L*
     & q1tq2t*q2tk2t*P2 - 4*k1(vb)*q1L*q1tq2t*P1*m1*m2 + 8*k1(vb)*
     & q1L*q1tq2t*P1*q2k2 - 4*k1(vb)*q1L*q1tq2t*P1*k1k2 + 4*k1(vb)*
     & q1L*q1tq2t*P2*m1*m2 + 8*k1(vb)*q1L*q1tq2t*P2*q2k2 + 4*k1(vb)
     & *q1L*q1tq2t*P2*k1k2 + 4*k1(vb)*q1L*q1tk2t*t2*P1 + 4*k1(vb)*
     & q1L*q1tk2t*t2*P2 - 4*k1(vb)*q1L*q2tk2t*t1*P1 - 4*k1(vb)*q1L
     & *q2tk2t*t1*P2 - 8*k1(vb)*q2L*q1tq2t*q1tk2t*P1 - 8*k1(vb)*
     & q2L*q1tq2t*q1tk2t*P2 + 8*k1(vb)*q2L*q1tq2t*q2tk2t*P1 + 8*
     & k1(vb)*q2L*q1tq2t*q2tk2t*P2 + 4*k1(vb)*q2L*q1tq2t*P1*m1*m2
     &  - 8*k1(vb)*q2L*q1tq2t*P1*q1k2 + 4*k1(vb)*q2L*q1tq2t*P1*k1k2
     &  - 4*k1(vb)*q2L*q1tq2t*P2*m1*m2 - 8*k1(vb)*q2L*q1tq2t*P2*
     & q1k2 - 4*k1(vb)*q2L*q1tq2t*P2*k1k2 + 4*k1(vb)*q2L*q1tk2t*t2*
     & P1 + 4*k1(vb)*q2L*q1tk2t*t2*P2 - 4*k1(vb)*q2L*q2tk2t*t1*P1
     &  - 4*k1(vb)*q2L*q2tk2t*t1*P2
      A8(3,vb) = A8(3,vb) + 8*k1(vb)*q1tq2t*q1tL*P1*m1*m2 - 8*
     & k1(vb)*q1tq2t*q1tL*P1*q1k2 - 8*k1(vb)*q1tq2t*q1tL*P1*q2k2 + 
     & 8*k1(vb)*q1tq2t*q1tL*P1*k1k2 - 8*k1(vb)*q1tq2t*q1tL*P2*m1*m2
     &  - 8*k1(vb)*q1tq2t*q1tL*P2*q1k2 - 8*k1(vb)*q1tq2t*q1tL*P2*
     & q2k2 - 8*k1(vb)*q1tq2t*q1tL*P2*k1k2 - 8*k1(vb)*q1tq2t*q2tL*
     & P1*m1*m2 + 8*k1(vb)*q1tq2t*q2tL*P1*q1k2 + 8*k1(vb)*q1tq2t*
     & q2tL*P1*q2k2 - 8*k1(vb)*q1tq2t*q2tL*P1*k1k2 + 8*k1(vb)*
     & q1tq2t*q2tL*P2*m1*m2 + 8*k1(vb)*q1tq2t*q2tL*P2*q1k2 + 8*
     & k1(vb)*q1tq2t*q2tL*P2*q2k2 + 8*k1(vb)*q1tq2t*q2tL*P2*k1k2 - 
     & 4*k1(vb)*q1tL*t2*P1*m1*m2 + 4*k1(vb)*q1tL*t2*P1*q1k2 + 4*
     & k1(vb)*q1tL*t2*P1*q2k2 - 4*k1(vb)*q1tL*t2*P1*k1k2 + 4*k1(vb)*
     & q1tL*t2*P2*m1*m2 + 4*k1(vb)*q1tL*t2*P2*q1k2 + 4*k1(vb)*q1tL*
     & t2*P2*q2k2 + 4*k1(vb)*q1tL*t2*P2*k1k2 + 4*k1(vb)*q2tL*t1*P1*m1
     & *m2 - 4*k1(vb)*q2tL*t1*P1*q1k2 - 4*k1(vb)*q2tL*t1*P1*q2k2 + 4*
     & k1(vb)*q2tL*t1*P1*k1k2 - 4*k1(vb)*q2tL*t1*P2*m1*m2 - 4*k1(vb)*
     & q2tL*t1*P2*q1k2
      A8(3,vb) = A8(3,vb) - 4*k1(vb)*q2tL*t1*P2*q2k2 - 4*k1(vb)*
     & q2tL*t1*P2*k1k2 + 8*k2(vb)*q1L*q1tq2t*q1tk1t*P1 + 8*k2(vb)*
     & q1L*q1tq2t*q1tk1t*P2 - 8*k2(vb)*q1L*q1tq2t*q2tk1t*P1 - 8*
     & k2(vb)*q1L*q1tq2t*q2tk1t*P2 - 4*k2(vb)*q1L*q1tq2t*P1*m1*m2
     &  - 8*k2(vb)*q1L*q1tq2t*P1*q2k1 - 4*k2(vb)*q1L*q1tq2t*P1*k1k2
     &  + 4*k2(vb)*q1L*q1tq2t*P2*m1*m2 - 8*k2(vb)*q1L*q1tq2t*P2*
     & q2k1 + 4*k2(vb)*q1L*q1tq2t*P2*k1k2 - 4*k2(vb)*q1L*q1tk1t*t2*
     & P1 - 4*k2(vb)*q1L*q1tk1t*t2*P2 + 4*k2(vb)*q1L*q2tk1t*t1*P1
     &  + 4*k2(vb)*q1L*q2tk1t*t1*P2 + 8*k2(vb)*q2L*q1tq2t*q1tk1t*
     & P1 + 8*k2(vb)*q2L*q1tq2t*q1tk1t*P2 - 8*k2(vb)*q2L*q1tq2t*
     & q2tk1t*P1 - 8*k2(vb)*q2L*q1tq2t*q2tk1t*P2 + 4*k2(vb)*q2L*
     & q1tq2t*P1*m1*m2 + 8*k2(vb)*q2L*q1tq2t*P1*q1k1 + 4*k2(vb)*q2L
     & *q1tq2t*P1*k1k2 - 4*k2(vb)*q2L*q1tq2t*P2*m1*m2 + 8*k2(vb)*
     & q2L*q1tq2t*P2*q1k1 - 4*k2(vb)*q2L*q1tq2t*P2*k1k2 - 4*k2(vb)*
     & q2L*q1tk1t*t2*P1 - 4*k2(vb)*q2L*q1tk1t*t2*P2 + 4*k2(vb)*q2L
     & *q2tk1t*t1*P1
      A8(3,vb) = A8(3,vb) + 4*k2(vb)*q2L*q2tk1t*t1*P2 + 8*k2(vb)*
     & q1tq2t*q1tL*P1*m1*m2 + 8*k2(vb)*q1tq2t*q1tL*P1*q1k1 + 8*
     & k2(vb)*q1tq2t*q1tL*P1*q2k1 + 8*k2(vb)*q1tq2t*q1tL*P1*k1k2 - 
     & 8*k2(vb)*q1tq2t*q1tL*P2*m1*m2 + 8*k2(vb)*q1tq2t*q1tL*P2*q1k1
     &  + 8*k2(vb)*q1tq2t*q1tL*P2*q2k1 - 8*k2(vb)*q1tq2t*q1tL*P2*
     & k1k2 - 8*k2(vb)*q1tq2t*q2tL*P1*m1*m2 - 8*k2(vb)*q1tq2t*q2tL*
     & P1*q1k1 - 8*k2(vb)*q1tq2t*q2tL*P1*q2k1 - 8*k2(vb)*q1tq2t*
     & q2tL*P1*k1k2 + 8*k2(vb)*q1tq2t*q2tL*P2*m1*m2 - 8*k2(vb)*
     & q1tq2t*q2tL*P2*q1k1 - 8*k2(vb)*q1tq2t*q2tL*P2*q2k1 + 8*
     & k2(vb)*q1tq2t*q2tL*P2*k1k2 - 4*k2(vb)*q1tL*t2*P1*m1*m2 - 4*
     & k2(vb)*q1tL*t2*P1*q1k1 - 4*k2(vb)*q1tL*t2*P1*q2k1 - 4*k2(vb)*
     & q1tL*t2*P1*k1k2 + 4*k2(vb)*q1tL*t2*P2*m1*m2 - 4*k2(vb)*q1tL*
     & t2*P2*q1k1 - 4*k2(vb)*q1tL*t2*P2*q2k1 + 4*k2(vb)*q1tL*t2*P2*
     & k1k2 + 4*k2(vb)*q2tL*t1*P1*m1*m2 + 4*k2(vb)*q2tL*t1*P1*q1k1 + 
     & 4*k2(vb)*q2tL*t1*P1*q2k1 + 4*k2(vb)*q2tL*t1*P1*k1k2 - 4*k2(vb)
     & *q2tL*t1*P2*m1*m2
      A8(3,vb) = A8(3,vb) + 4*k2(vb)*q2tL*t1*P2*q1k1 + 4*k2(vb)*
     & q2tL*t1*P2*q2k1 - 4*k2(vb)*q2tL*t1*P2*k1k2 - 4*kw(vb)*q1L*
     & q1tq2t*P1*m1*m2 - 4*kw(vb)*q1L*q1tq2t*P1*k1k2 + 4*kw(vb)*q1L
     & *q1tq2t*P2*m1*m2 + 4*kw(vb)*q1L*q1tq2t*P2*k1k2 + 4*kw(vb)*
     & q2L*q1tq2t*P1*m1*m2 + 4*kw(vb)*q2L*q1tq2t*P1*k1k2 - 4*kw(vb)
     & *q2L*q1tq2t*P2*m1*m2 - 4*kw(vb)*q2L*q1tq2t*P2*k1k2 + 8*
     & kw(vb)*q1tq2t*q1tL*P1*m1*m2 + 8*kw(vb)*q1tq2t*q1tL*P1*k1k2
     &  - 8*kw(vb)*q1tq2t*q1tL*P2*m1*m2 - 8*kw(vb)*q1tq2t*q1tL*P2*
     & k1k2 - 8*kw(vb)*q1tq2t*q2tL*P1*m1*m2 - 8*kw(vb)*q1tq2t*q2tL*
     & P1*k1k2 + 8*kw(vb)*q1tq2t*q2tL*P2*m1*m2 + 8*kw(vb)*q1tq2t*
     & q2tL*P2*k1k2 - 4*kw(vb)*q1tL*t2*P1*m1*m2 - 4*kw(vb)*q1tL*t2*
     & P1*k1k2 + 4*kw(vb)*q1tL*t2*P2*m1*m2 + 4*kw(vb)*q1tL*t2*P2*k1k2
     &  + 4*kw(vb)*q2tL*t1*P1*m1*m2 + 4*kw(vb)*q2tL*t1*P1*k1k2 - 4*
     & kw(vb)*q2tL*t1*P2*m1*m2 - 4*kw(vb)*q2tL*t1*P2*k1k2 + 4*q1(vb)*
     & q1L*q1tq2t*P1*m1*m2 + 4*q1(vb)*q1L*q1tq2t*P1*k1k2 - 4*q1(vb)
     & *q1L*q1tq2t*P2*m1*m2
      A8(3,vb) = A8(3,vb) - 4*q1(vb)*q1L*q1tq2t*P2*k1k2 + 4*q1(vb)*
     & q2L*q1tq2t*P1*m1*m2 + 4*q1(vb)*q2L*q1tq2t*P1*k1k2 - 4*q1(vb)
     & *q2L*q1tq2t*P2*m1*m2 - 4*q1(vb)*q2L*q1tq2t*P2*k1k2 - 4*
     & q2(vb)*q1L*q1tq2t*P1*m1*m2 - 4*q2(vb)*q1L*q1tq2t*P1*k1k2 + 4
     & *q2(vb)*q1L*q1tq2t*P2*m1*m2 + 4*q2(vb)*q1L*q1tq2t*P2*k1k2 - 
     & 4*q2(vb)*q2L*q1tq2t*P1*m1*m2 - 4*q2(vb)*q2L*q1tq2t*P1*k1k2
     &  + 4*q2(vb)*q2L*q1tq2t*P2*m1*m2 + 4*q2(vb)*q2L*q1tq2t*P2*
     & k1k2 - 8*q1t(vb)*q1L*q1tq2t*P1*m1*m2 - 8*q1t(vb)*q1L*q1tq2t*
     & P1*k1k2 + 8*q1t(vb)*q1L*q1tq2t*P2*m1*m2 + 8*q1t(vb)*q1L*
     & q1tq2t*P2*k1k2 + 4*q1t(vb)*q1L*t2*P1*m1*m2 + 4*q1t(vb)*q1L*t2
     & *P1*k1k2 - 4*q1t(vb)*q1L*t2*P2*m1*m2 - 4*q1t(vb)*q1L*t2*P2*
     & k1k2 - 8*q1t(vb)*q2L*q1tq2t*P1*m1*m2 - 8*q1t(vb)*q2L*q1tq2t*
     & P1*k1k2 + 8*q1t(vb)*q2L*q1tq2t*P2*m1*m2 + 8*q1t(vb)*q2L*
     & q1tq2t*P2*k1k2 + 4*q1t(vb)*q2L*t2*P1*m1*m2 + 4*q1t(vb)*q2L*t2
     & *P1*k1k2 - 4*q1t(vb)*q2L*t2*P2*m1*m2 - 4*q1t(vb)*q2L*t2*P2*
     & k1k2
      A8(3,vb) = A8(3,vb) + 8*q2t(vb)*q1L*q1tq2t*P1*m1*m2 + 8*
     & q2t(vb)*q1L*q1tq2t*P1*k1k2 - 8*q2t(vb)*q1L*q1tq2t*P2*m1*m2
     &  - 8*q2t(vb)*q1L*q1tq2t*P2*k1k2 - 4*q2t(vb)*q1L*t1*P1*m1*m2
     &  - 4*q2t(vb)*q1L*t1*P1*k1k2 + 4*q2t(vb)*q1L*t1*P2*m1*m2 + 4*
     & q2t(vb)*q1L*t1*P2*k1k2 + 8*q2t(vb)*q2L*q1tq2t*P1*m1*m2 + 8*
     & q2t(vb)*q2L*q1tq2t*P1*k1k2 - 8*q2t(vb)*q2L*q1tq2t*P2*m1*m2
     &  - 8*q2t(vb)*q2L*q1tq2t*P2*k1k2 - 4*q2t(vb)*q2L*t1*P1*m1*m2
     &  - 4*q2t(vb)*q2L*t1*P1*k1k2 + 4*q2t(vb)*q2L*t1*P2*m1*m2 + 4*
     & q2t(vb)*q2L*t1*P2*k1k2 + 16*L(vb)*q1tq2t*q1tk1t*P1*m1*m2 - 8*
     & L(vb)*q1tq2t*q1tk1t*P1*q1k2 - 8*L(vb)*q1tq2t*q1tk1t*P1*q2k2
     &  + 16*L(vb)*q1tq2t*q1tk1t*P1*k1k2 - 8*L(vb)*q1tq2t*q1tk1t*P2
     & *q1k2 - 8*L(vb)*q1tq2t*q1tk1t*P2*q2k2 + 8*L(vb)*q1tq2t*
     & q1tk2t*P1*q1k1 + 8*L(vb)*q1tq2t*q1tk2t*P1*q2k1 - 16*L(vb)*
     & q1tq2t*q1tk2t*P2*m1*m2 + 8*L(vb)*q1tq2t*q1tk2t*P2*q1k1 + 8*
     & L(vb)*q1tq2t*q1tk2t*P2*q2k1 - 16*L(vb)*q1tq2t*q1tk2t*P2*k1k2
     &  - 16*L(vb)*q1tq2t*q2tk1t*P1*m1*m2
      A8(3,vb) = A8(3,vb) + 8*L(vb)*q1tq2t*q2tk1t*P1*q1k2 + 8*L(vb)
     & *q1tq2t*q2tk1t*P1*q2k2 - 16*L(vb)*q1tq2t*q2tk1t*P1*k1k2 + 8*
     & L(vb)*q1tq2t*q2tk1t*P2*q1k2 + 8*L(vb)*q1tq2t*q2tk1t*P2*q2k2
     &  - 8*L(vb)*q1tq2t*q2tk2t*P1*q1k1 - 8*L(vb)*q1tq2t*q2tk2t*P1*
     & q2k1 + 16*L(vb)*q1tq2t*q2tk2t*P2*m1*m2 - 8*L(vb)*q1tq2t*
     & q2tk2t*P2*q1k1 - 8*L(vb)*q1tq2t*q2tk2t*P2*q2k1 + 16*L(vb)*
     & q1tq2t*q2tk2t*P2*k1k2 + 8*L(vb)*q1tq2t*P1*m1*m2*q1k1 - 8*
     & L(vb)*q1tq2t*P1*m1*m2*q2k1 - 8*L(vb)*q1tq2t*P1*q1k1*q2k2 + 8*
     & L(vb)*q1tq2t*P1*q1k1*k1k2 + 8*L(vb)*q1tq2t*P1*q1k2*q2k1 - 8*
     & L(vb)*q1tq2t*P1*q2k1*k1k2 - 8*L(vb)*q1tq2t*P2*m1*m2*q1k2 + 8*
     & L(vb)*q1tq2t*P2*m1*m2*q2k2 - 8*L(vb)*q1tq2t*P2*q1k1*q2k2 + 8*
     & L(vb)*q1tq2t*P2*q1k2*q2k1 - 8*L(vb)*q1tq2t*P2*q1k2*k1k2 + 8*
     & L(vb)*q1tq2t*P2*q2k2*k1k2 - 8*L(vb)*q1tk1t*t2*P1*m1*m2 + 4*
     & L(vb)*q1tk1t*t2*P1*q1k2 + 4*L(vb)*q1tk1t*t2*P1*q2k2 - 8*L(vb)*
     & q1tk1t*t2*P1*k1k2 + 4*L(vb)*q1tk1t*t2*P2*q1k2 + 4*L(vb)*
     & q1tk1t*t2*P2*q2k2
      A8(3,vb) = A8(3,vb) - 4*L(vb)*q1tk2t*t2*P1*q1k1 - 4*L(vb)*
     & q1tk2t*t2*P1*q2k1 + 8*L(vb)*q1tk2t*t2*P2*m1*m2 - 4*L(vb)*
     & q1tk2t*t2*P2*q1k1 - 4*L(vb)*q1tk2t*t2*P2*q2k1 + 8*L(vb)*
     & q1tk2t*t2*P2*k1k2 + 8*L(vb)*q2tk1t*t1*P1*m1*m2 - 4*L(vb)*
     & q2tk1t*t1*P1*q1k2 - 4*L(vb)*q2tk1t*t1*P1*q2k2 + 8*L(vb)*
     & q2tk1t*t1*P1*k1k2 - 4*L(vb)*q2tk1t*t1*P2*q1k2 - 4*L(vb)*
     & q2tk1t*t1*P2*q2k2 + 4*L(vb)*q2tk2t*t1*P1*q1k1 + 4*L(vb)*
     & q2tk2t*t1*P1*q2k1 - 8*L(vb)*q2tk2t*t1*P2*m1*m2 + 4*L(vb)*
     & q2tk2t*t1*P2*q1k1 + 4*L(vb)*q2tk2t*t1*P2*q2k1 - 8*L(vb)*
     & q2tk2t*t1*P2*k1k2

      A8(4,vb) =
     &  - 8*k1(vb)*q1L*q2U*q1tq2t*P1*m2 - 8*k1(vb)*q1L*q2U*q1tq2t
     & *P2*m2 - 8*k1(vb)*q1L*q1tq2t*q1tU*P1*m2 - 8*k1(vb)*q1L*
     & q1tq2t*q1tU*P2*m2 + 8*k1(vb)*q1L*q1tq2t*q2tU*P1*m2 + 8*
     & k1(vb)*q1L*q1tq2t*q2tU*P2*m2 + 4*k1(vb)*q1L*q1tU*t2*P1*m2
     &  + 4*k1(vb)*q1L*q1tU*t2*P2*m2 - 4*k1(vb)*q1L*q2tU*t1*P1*m2
     &  - 4*k1(vb)*q1L*q2tU*t1*P2*m2 - 8*k1(vb)*q2U*q2L*q1tq2t*P1*
     & m2 - 8*k1(vb)*q2U*q2L*q1tq2t*P2*m2 - 8*k1(vb)*q2L*q1tq2t*
     & q1tU*P1*m2 - 8*k1(vb)*q2L*q1tq2t*q1tU*P2*m2 + 8*k1(vb)*q2L*
     & q1tq2t*q2tU*P1*m2 + 8*k1(vb)*q2L*q1tq2t*q2tU*P2*m2 + 4*
     & k1(vb)*q2L*q1tU*t2*P1*m2 + 4*k1(vb)*q2L*q1tU*t2*P2*m2 - 4*
     & k1(vb)*q2L*q2tU*t1*P1*m2 - 4*k1(vb)*q2L*q2tU*t1*P2*m2 + 8*
     & k2(vb)*q1L*q2U*q1tq2t*P1*m1 + 8*k2(vb)*q1L*q2U*q1tq2t*P2*
     & m1 + 8*k2(vb)*q1L*q1tq2t*q1tU*P1*m1 + 8*k2(vb)*q1L*q1tq2t*
     & q1tU*P2*m1 - 8*k2(vb)*q1L*q1tq2t*q2tU*P1*m1 - 8*k2(vb)*q1L*
     & q1tq2t*q2tU*P2*m1
      A8(4,vb) = A8(4,vb) - 4*k2(vb)*q1L*q1tU*t2*P1*m1 - 4*k2(vb)*
     & q1L*q1tU*t2*P2*m1 + 4*k2(vb)*q1L*q2tU*t1*P1*m1 + 4*k2(vb)*
     & q1L*q2tU*t1*P2*m1 + 8*k2(vb)*q2U*q2L*q1tq2t*P1*m1 + 8*
     & k2(vb)*q2U*q2L*q1tq2t*P2*m1 + 8*k2(vb)*q2L*q1tq2t*q1tU*P1*
     & m1 + 8*k2(vb)*q2L*q1tq2t*q1tU*P2*m1 - 8*k2(vb)*q2L*q1tq2t*
     & q2tU*P1*m1 - 8*k2(vb)*q2L*q1tq2t*q2tU*P2*m1 - 4*k2(vb)*q2L*
     & q1tU*t2*P1*m1 - 4*k2(vb)*q2L*q1tU*t2*P2*m1 + 4*k2(vb)*q2L*
     & q2tU*t1*P1*m1 + 4*k2(vb)*q2L*q2tU*t1*P2*m1 - 8*U(vb)*q1L*
     & q1tq2t*q1tk1t*P1*m2 - 8*U(vb)*q1L*q1tq2t*q1tk1t*P2*m2 + 8*
     & U(vb)*q1L*q1tq2t*q1tk2t*P1*m1 + 8*U(vb)*q1L*q1tq2t*q1tk2t*
     & P2*m1 + 8*U(vb)*q1L*q1tq2t*q2tk1t*P1*m2 + 8*U(vb)*q1L*
     & q1tq2t*q2tk1t*P2*m2 - 8*U(vb)*q1L*q1tq2t*q2tk2t*P1*m1 - 8*
     & U(vb)*q1L*q1tq2t*q2tk2t*P2*m1 - 8*U(vb)*q1L*q1tq2t*P1*m1*
     & q2k2 + 8*U(vb)*q1L*q1tq2t*P1*m2*q2k1 - 8*U(vb)*q1L*q1tq2t*P2
     & *m1*q2k2 + 8*U(vb)*q1L*q1tq2t*P2*m2*q2k1 + 4*U(vb)*q1L*
     & q1tk1t*t2*P1*m2
      A8(4,vb) = A8(4,vb) + 4*U(vb)*q1L*q1tk1t*t2*P2*m2 - 4*U(vb)*
     & q1L*q1tk2t*t2*P1*m1 - 4*U(vb)*q1L*q1tk2t*t2*P2*m1 - 4*U(vb)*
     & q1L*q2tk1t*t1*P1*m2 - 4*U(vb)*q1L*q2tk1t*t1*P2*m2 + 4*U(vb)*
     & q1L*q2tk2t*t1*P1*m1 + 4*U(vb)*q1L*q2tk2t*t1*P2*m1 - 8*U(vb)*
     & q2L*q1tq2t*q1tk1t*P1*m2 - 8*U(vb)*q2L*q1tq2t*q1tk1t*P2*m2
     &  + 8*U(vb)*q2L*q1tq2t*q1tk2t*P1*m1 + 8*U(vb)*q2L*q1tq2t*
     & q1tk2t*P2*m1 + 8*U(vb)*q2L*q1tq2t*q2tk1t*P1*m2 + 8*U(vb)*
     & q2L*q1tq2t*q2tk1t*P2*m2 - 8*U(vb)*q2L*q1tq2t*q2tk2t*P1*m1
     &  - 8*U(vb)*q2L*q1tq2t*q2tk2t*P2*m1 + 8*U(vb)*q2L*q1tq2t*P1*
     & m1*q1k2 - 8*U(vb)*q2L*q1tq2t*P1*m2*q1k1 + 8*U(vb)*q2L*q1tq2t
     & *P2*m1*q1k2 - 8*U(vb)*q2L*q1tq2t*P2*m2*q1k1 + 4*U(vb)*q2L*
     & q1tk1t*t2*P1*m2 + 4*U(vb)*q2L*q1tk1t*t2*P2*m2 - 4*U(vb)*q2L*
     & q1tk2t*t2*P1*m1 - 4*U(vb)*q2L*q1tk2t*t2*P2*m1 - 4*U(vb)*q2L*
     & q2tk1t*t1*P1*m2 - 4*U(vb)*q2L*q2tk1t*t1*P2*m2 + 4*U(vb)*q2L*
     & q2tk2t*t1*P1*m1 + 4*U(vb)*q2L*q2tk2t*t1*P2*m1 + 8*U(vb)*
     & q1tq2t*q1tL*P1*m1*q1k2
      A8(4,vb) = A8(4,vb) + 8*U(vb)*q1tq2t*q1tL*P1*m1*q2k2 - 8*
     & U(vb)*q1tq2t*q1tL*P1*m2*q1k1 - 8*U(vb)*q1tq2t*q1tL*P1*m2*
     & q2k1 + 8*U(vb)*q1tq2t*q1tL*P2*m1*q1k2 + 8*U(vb)*q1tq2t*q1tL*
     & P2*m1*q2k2 - 8*U(vb)*q1tq2t*q1tL*P2*m2*q1k1 - 8*U(vb)*q1tq2t*
     & q1tL*P2*m2*q2k1 - 8*U(vb)*q1tq2t*q2tL*P1*m1*q1k2 - 8*U(vb)*
     & q1tq2t*q2tL*P1*m1*q2k2 + 8*U(vb)*q1tq2t*q2tL*P1*m2*q1k1 + 8*
     & U(vb)*q1tq2t*q2tL*P1*m2*q2k1 - 8*U(vb)*q1tq2t*q2tL*P2*m1*
     & q1k2 - 8*U(vb)*q1tq2t*q2tL*P2*m1*q2k2 + 8*U(vb)*q1tq2t*q2tL*
     & P2*m2*q1k1 + 8*U(vb)*q1tq2t*q2tL*P2*m2*q2k1 - 4*U(vb)*q1tL*t2
     & *P1*m1*q1k2 - 4*U(vb)*q1tL*t2*P1*m1*q2k2 + 4*U(vb)*q1tL*t2*P1*
     & m2*q1k1 + 4*U(vb)*q1tL*t2*P1*m2*q2k1 - 4*U(vb)*q1tL*t2*P2*m1*
     & q1k2 - 4*U(vb)*q1tL*t2*P2*m1*q2k2 + 4*U(vb)*q1tL*t2*P2*m2*q1k1
     &  + 4*U(vb)*q1tL*t2*P2*m2*q2k1 + 4*U(vb)*q2tL*t1*P1*m1*q1k2 + 4
     & *U(vb)*q2tL*t1*P1*m1*q2k2 - 4*U(vb)*q2tL*t1*P1*m2*q1k1 - 4*
     & U(vb)*q2tL*t1*P1*m2*q2k1 + 4*U(vb)*q2tL*t1*P2*m1*q1k2 + 4*
     & U(vb)*q2tL*t1*P2*m1*q2k2
      A8(4,vb) = A8(4,vb) - 4*U(vb)*q2tL*t1*P2*m2*q1k1 - 4*U(vb)*
     & q2tL*t1*P2*m2*q2k1 - 8*L(vb)*q2U*q1tq2t*P1*m1*q1k2 - 8*L(vb)*
     & q2U*q1tq2t*P1*m1*q2k2 + 8*L(vb)*q2U*q1tq2t*P1*m2*q1k1 + 8*
     & L(vb)*q2U*q1tq2t*P1*m2*q2k1 - 8*L(vb)*q2U*q1tq2t*P2*m1*q1k2
     &  - 8*L(vb)*q2U*q1tq2t*P2*m1*q2k2 + 8*L(vb)*q2U*q1tq2t*P2*m2*
     & q1k1 + 8*L(vb)*q2U*q1tq2t*P2*m2*q2k1 - 8*L(vb)*q1tq2t*q1tU*
     & P1*m1*q1k2 - 8*L(vb)*q1tq2t*q1tU*P1*m1*q2k2 + 8*L(vb)*q1tq2t*
     & q1tU*P1*m2*q1k1 + 8*L(vb)*q1tq2t*q1tU*P1*m2*q2k1 - 8*L(vb)*
     & q1tq2t*q1tU*P2*m1*q1k2 - 8*L(vb)*q1tq2t*q1tU*P2*m1*q2k2 + 8*
     & L(vb)*q1tq2t*q1tU*P2*m2*q1k1 + 8*L(vb)*q1tq2t*q1tU*P2*m2*
     & q2k1 + 8*L(vb)*q1tq2t*q2tU*P1*m1*q1k2 + 8*L(vb)*q1tq2t*q2tU*
     & P1*m1*q2k2 - 8*L(vb)*q1tq2t*q2tU*P1*m2*q1k1 - 8*L(vb)*q1tq2t*
     & q2tU*P1*m2*q2k1 + 8*L(vb)*q1tq2t*q2tU*P2*m1*q1k2 + 8*L(vb)*
     & q1tq2t*q2tU*P2*m1*q2k2 - 8*L(vb)*q1tq2t*q2tU*P2*m2*q1k1 - 8*
     & L(vb)*q1tq2t*q2tU*P2*m2*q2k1 + 4*L(vb)*q1tU*t2*P1*m1*q1k2 + 4
     & *L(vb)*q1tU*t2*P1*m1*q2k2
      A8(4,vb) = A8(4,vb) - 4*L(vb)*q1tU*t2*P1*m2*q1k1 - 4*L(vb)*
     & q1tU*t2*P1*m2*q2k1 + 4*L(vb)*q1tU*t2*P2*m1*q1k2 + 4*L(vb)*
     & q1tU*t2*P2*m1*q2k2 - 4*L(vb)*q1tU*t2*P2*m2*q1k1 - 4*L(vb)*
     & q1tU*t2*P2*m2*q2k1 - 4*L(vb)*q2tU*t1*P1*m1*q1k2 - 4*L(vb)*
     & q2tU*t1*P1*m1*q2k2 + 4*L(vb)*q2tU*t1*P1*m2*q1k1 + 4*L(vb)*
     & q2tU*t1*P1*m2*q2k1 - 4*L(vb)*q2tU*t1*P2*m1*q1k2 - 4*L(vb)*
     & q2tU*t1*P2*m1*q2k2 + 4*L(vb)*q2tU*t1*P2*m2*q1k1 + 4*L(vb)*
     & q2tU*t1*P2*m2*q2k1

      AS1(1,vb) =
     & 4*k1(vb)*t1*t2*P3*m2 + 4*k1(vb)*t1*t2*P4*m2 - 4*k1(vb)*t1*t2*P7*
     & m2 - 12*k1(vb)*q1tq2t*q1tk1t*P3*m2 - 4*k1(vb)*q1tq2t*q1tk2t*
     & P3*m1 + 4*k1(vb)*q1tq2t*q2tk1t*P4*m2 - 4*k1(vb)*q1tq2t*
     & q2tk2t*P4*m1 - 16*k1(vb)*q1tq2t*q2tk2t*P4*m2 - 4*k1(vb)*
     & q1tq2t*t1*P3*m2 - 4*k1(vb)*q1tq2t*t2*P4*m2 - 4*k1(vb)*q1tq2t*
     & P3*m1*q1k2 + 4*k1(vb)*q1tq2t*P3*m2*q1q2 - 4*k1(vb)*q1tq2t*P3*
     & m2*q1k1 - 4*k1(vb)*q1tq2t*P4*m1*q2k2 + 4*k1(vb)*q1tq2t*P4*m2*
     & q1q2 + 4*k1(vb)*q1tq2t*P4*m2*q2k1 - 8*k1(vb)*q1tq2t*P4*m2*q2k2
     &  + 4*k1(vb)*q1tq2t*P7*m2*q1q2 + 4*k1(vb)*q1tq2t**2*P3*m2 + 4*
     & k1(vb)*q1tq2t**2*P4*m2 + 4*k1(vb)*q1tq2t**2*P7*m2 + 8*k1(vb)*
     & q1tk1t*q2tk1t*P3*m2 + 8*k1(vb)*q1tk1t*q2tk2t*P3*m1 - 8*
     & k1(vb)*q1tk1t*q2tk2t*P4*m2 - 16*k1(vb)*q1tk1t*q2tk2t*P7*m2
     &  + 8*k1(vb)*q1tk1t*t2*P3*m2 - 4*k1(vb)*q1tk1t*t2*P4*m2 - 8*
     & k1(vb)*q1tk1t*t2*P7*m2 + 8*k1(vb)*q1tk2t*q2tk2t*P4*m1 + 16*
     & k1(vb)*q1tk2t*q2tk2t*P4*m2
      AS1(1,vb) = AS1(1,vb) + 4*k1(vb)*q1tk2t*t2*P4*m1 + 8*k1(vb)*
     & q1tk2t*t2*P4*m2 + 4*k1(vb)*q2tk1t*t1*P3*m2 + 4*k1(vb)*q2tk2t*
     & t1*P3*m1 + 8*k1(vb)*q2tk2t*t1*P4*m2 - 8*k1(vb)*q2tk2t*t1*P7*m2
     &  - 4*k2(vb)*t1*t2*P3*m1 - 4*k2(vb)*t1*t2*P4*m1 + 4*k2(vb)*t1*t2*
     & P7*m1 + 16*k2(vb)*q1tq2t*q1tk1t*P3*m1 + 4*k2(vb)*q1tq2t*
     & q1tk1t*P3*m2 - 4*k2(vb)*q1tq2t*q1tk2t*P3*m1 + 4*k2(vb)*
     & q1tq2t*q2tk1t*P4*m2 + 12*k2(vb)*q1tq2t*q2tk2t*P4*m1 + 4*
     & k2(vb)*q1tq2t*t1*P3*m1 + 4*k2(vb)*q1tq2t*t2*P4*m1 - 4*k2(vb)*
     & q1tq2t*P3*m1*q1q2 + 8*k2(vb)*q1tq2t*P3*m1*q1k1 - 4*k2(vb)*
     & q1tq2t*P3*m1*q1k2 + 4*k2(vb)*q1tq2t*P3*m2*q1k1 - 4*k2(vb)*
     & q1tq2t*P4*m1*q1q2 + 4*k2(vb)*q1tq2t*P4*m1*q2k2 + 4*k2(vb)*
     & q1tq2t*P4*m2*q2k1 - 4*k2(vb)*q1tq2t*P7*m1*q1q2 - 4*k2(vb)*
     & q1tq2t**2*P3*m1 - 4*k2(vb)*q1tq2t**2*P4*m1 - 4*k2(vb)*
     & q1tq2t**2*P7*m1 - 16*k2(vb)*q1tk1t*q2tk1t*P3*m1 - 8*k2(vb)*
     & q1tk1t*q2tk1t*P3*m2 + 8*k2(vb)*q1tk1t*q2tk2t*P3*m1 - 8*
     & k2(vb)*q1tk1t*q2tk2t*P4*m2
      AS1(1,vb) = AS1(1,vb) + 16*k2(vb)*q1tk1t*q2tk2t*P7*m1 - 8*
     & k2(vb)*q1tk1t*t2*P3*m1 - 4*k2(vb)*q1tk1t*t2*P4*m2 + 8*k2(vb)*
     & q1tk1t*t2*P7*m1 - 8*k2(vb)*q1tk2t*q2tk2t*P4*m1 - 4*k2(vb)*
     & q1tk2t*t2*P4*m1 - 8*k2(vb)*q2tk1t*t1*P3*m1 - 4*k2(vb)*q2tk1t*
     & t1*P3*m2 + 4*k2(vb)*q2tk2t*t1*P3*m1 - 8*k2(vb)*q2tk2t*t1*P4*m1
     &  + 8*k2(vb)*q2tk2t*t1*P7*m1 + 4*kw(vb)*q1tq2t*q1tk1t*P3*m2 - 
     & 4*kw(vb)*q1tq2t*q1tk2t*P3*m1 + 4*kw(vb)*q1tq2t*q2tk1t*P4*m2
     &  - 4*kw(vb)*q1tq2t*q2tk2t*P4*m1 - 4*kw(vb)*q1tq2t*P3*m1*q1k2
     &  + 4*kw(vb)*q1tq2t*P3*m2*q1k1 - 4*kw(vb)*q1tq2t*P4*m1*q2k2 + 4
     & *kw(vb)*q1tq2t*P4*m2*q2k1 - 8*kw(vb)*q1tk1t*q2tk1t*P3*m2 + 8*
     & kw(vb)*q1tk1t*q2tk2t*P3*m1 - 8*kw(vb)*q1tk1t*q2tk2t*P4*m2 - 
     & 4*kw(vb)*q1tk1t*t2*P4*m2 + 8*kw(vb)*q1tk2t*q2tk2t*P4*m1 + 4*
     & kw(vb)*q1tk2t*t2*P4*m1 - 4*kw(vb)*q2tk1t*t1*P3*m2 + 4*kw(vb)*
     & q2tk2t*t1*P3*m1 + 4*q1(vb)*q1tq2t*q1tk1t*P3*m2 - 4*q1(vb)*
     & q1tq2t*q1tk2t*P3*m1 - 4*q1(vb)*q1tq2t*q2tk1t*P3*m2 - 4*
     & q1(vb)*q1tq2t*q2tk1t*P7*m2
      AS1(1,vb) = AS1(1,vb) + 4*q1(vb)*q1tq2t*q2tk2t*P3*m1 + 4*
     & q1(vb)*q1tq2t*q2tk2t*P7*m1 + 4*q1(vb)*q1tq2t*P3*m1*q1k2 + 4*
     & q1(vb)*q1tq2t*P3*m1*q2k2 - 4*q1(vb)*q1tq2t*P3*m2*q1k1 - 4*
     & q1(vb)*q1tq2t*P3*m2*q2k1 + 4*q1(vb)*q1tq2t*P7*m1*q2k2 - 4*
     & q1(vb)*q1tq2t*P7*m2*q2k1 - 8*q1(vb)*q1tk1t*q2tk2t*P3*m1 + 8*
     & q1(vb)*q1tk1t*q2tk2t*P7*m2 - 4*q1(vb)*q1tk1t*t2*P3*m2 + 4*
     & q1(vb)*q1tk1t*t2*P7*m2 + 8*q1(vb)*q1tk2t*q2tk1t*P3*m1 - 8*
     & q1(vb)*q1tk2t*q2tk2t*P7*m1 + 4*q1(vb)*q1tk2t*t2*P3*m1 - 4*
     & q1(vb)*q1tk2t*t2*P7*m1 + 4*q1(vb)*q2tk1t*t1*P3*m2 - 4*q1(vb)*
     & q2tk2t*t1*P3*m1 - 4*q2(vb)*q1tq2t*q1tk1t*P4*m2 - 4*q2(vb)*
     & q1tq2t*q1tk1t*P7*m2 + 4*q2(vb)*q1tq2t*q1tk2t*P4*m1 + 4*
     & q2(vb)*q1tq2t*q1tk2t*P7*m1 + 4*q2(vb)*q1tq2t*q2tk1t*P4*m2 - 
     & 4*q2(vb)*q1tq2t*q2tk2t*P4*m1 + 4*q2(vb)*q1tq2t*P4*m1*q1k2 + 4
     & *q2(vb)*q1tq2t*P4*m1*q2k2 - 4*q2(vb)*q1tq2t*P4*m2*q1k1 - 4*
     & q2(vb)*q1tq2t*P4*m2*q2k1 + 4*q2(vb)*q1tq2t*P7*m1*q1k2 - 4*
     & q2(vb)*q1tq2t*P7*m2*q1k1
      AS1(1,vb) = AS1(1,vb) + 8*q2(vb)*q1tk1t*q2tk1t*P7*m2 + 8*
     & q2(vb)*q1tk1t*q2tk2t*P4*m2 - 8*q2(vb)*q1tk1t*q2tk2t*P7*m1 + 
     & 4*q2(vb)*q1tk1t*t2*P4*m2 - 8*q2(vb)*q1tk2t*q2tk1t*P4*m2 - 4*
     & q2(vb)*q1tk2t*t2*P4*m1 - 4*q2(vb)*q2tk1t*t1*P4*m2 + 4*q2(vb)*
     & q2tk1t*t1*P7*m2 + 4*q2(vb)*q2tk2t*t1*P4*m1 - 4*q2(vb)*q2tk2t*
     & t1*P7*m1 + 4*q1t(vb)*t2*P3*m1*q1k2 - 4*q1t(vb)*t2*P3*m2*q1k1 + 4
     & *q1t(vb)*t2*P4*m1*q1k2 - 4*q1t(vb)*t2*P4*m2*q1k1 - 4*q1t(vb)*t2*
     & P7*m1*q1k2 + 4*q1t(vb)*t2*P7*m2*q1k1 - 8*q1t(vb)*q1tq2t*P3*m1*
     & q1k2 - 4*q1t(vb)*q1tq2t*P3*m1*q2k2 + 8*q1t(vb)*q1tq2t*P3*m2*
     & q1k1 + 4*q1t(vb)*q1tq2t*P3*m2*q2k1 + 4*q1t(vb)*q1tq2t*P4*m1*
     & q2k2 - 4*q1t(vb)*q1tq2t*P4*m2*q2k1 - 4*q1t(vb)*q1tq2t*P7*m1*
     & q2k2 + 4*q1t(vb)*q1tq2t*P7*m2*q2k1 - 4*q1t(vb)*q2tk1t*t1*P3*m2
     &  + 4*q1t(vb)*q2tk1t*t2*P4*m2 + 8*q1t(vb)*q2tk1t*P3*m1*q1k2 - 4
     & *q1t(vb)*q2tk1t*P3*m2*q1q2 + 4*q1t(vb)*q2tk1t*P4*m2*q1q2 - 8*
     & q1t(vb)*q2tk1t*P4*m2*q2k2 - 4*q1t(vb)*q2tk1t*P7*m2*q1q2 + 4*
     & q1t(vb)*q2tk2t*t1*P3*m1
      AS1(1,vb) = AS1(1,vb) - 4*q1t(vb)*q2tk2t*t2*P4*m1 + 4*q1t(vb)*
     & q2tk2t*P3*m1*q1q2 - 8*q1t(vb)*q2tk2t*P3*m1*q1k1 - 4*q1t(vb)*
     & q2tk2t*P4*m1*q1q2 + 8*q1t(vb)*q2tk2t*P4*m1*q1k2 + 8*q1t(vb)*
     & q2tk2t*P4*m1*q2k2 - 8*q1t(vb)*q2tk2t*P4*m2*q1k1 + 4*q1t(vb)*
     & q2tk2t*P7*m1*q1q2 - 8*q1t(vb)*q2tk2t*P7*m1*q1k2 + 8*q1t(vb)*
     & q2tk2t*P7*m2*q1k1 + 4*q2t(vb)*t1*P3*m1*q2k2 - 4*q2t(vb)*t1*P3*
     & m2*q2k1 + 4*q2t(vb)*t1*P4*m1*q2k2 - 4*q2t(vb)*t1*P4*m2*q2k1 - 4*
     & q2t(vb)*t1*P7*m1*q2k2 + 4*q2t(vb)*t1*P7*m2*q2k1 + 4*q2t(vb)*
     & q1tq2t*P3*m1*q1k2 - 4*q2t(vb)*q1tq2t*P3*m2*q1k1 - 4*q2t(vb)*
     & q1tq2t*P4*m1*q1k2 - 8*q2t(vb)*q1tq2t*P4*m1*q2k2 + 4*q2t(vb)*
     & q1tq2t*P4*m2*q1k1 + 8*q2t(vb)*q1tq2t*P4*m2*q2k1 - 4*q2t(vb)*
     & q1tq2t*P7*m1*q1k2 + 4*q2t(vb)*q1tq2t*P7*m2*q1k1 + 4*q2t(vb)*
     & q1tk1t*t1*P3*m2 - 4*q2t(vb)*q1tk1t*t2*P4*m2 + 8*q2t(vb)*
     & q1tk1t*P3*m1*q2k2 + 4*q2t(vb)*q1tk1t*P3*m2*q1q2 - 8*q2t(vb)*
     & q1tk1t*P3*m2*q1k1 - 8*q2t(vb)*q1tk1t*P3*m2*q2k1 - 4*q2t(vb)*
     & q1tk1t*P4*m2*q1q2
      AS1(1,vb) = AS1(1,vb) + 8*q2t(vb)*q1tk1t*P4*m2*q2k2 - 8*
     & q2t(vb)*q1tk1t*P7*m1*q2k2 - 4*q2t(vb)*q1tk1t*P7*m2*q1q2 + 8*
     & q2t(vb)*q1tk1t*P7*m2*q2k1 - 4*q2t(vb)*q1tk2t*t1*P3*m1 + 4*
     & q2t(vb)*q1tk2t*t2*P4*m1 - 4*q2t(vb)*q1tk2t*P3*m1*q1q2 + 8*
     & q2t(vb)*q1tk2t*P3*m1*q1k1 + 4*q2t(vb)*q1tk2t*P4*m1*q1q2 - 8*
     & q2t(vb)*q1tk2t*P4*m2*q2k1 + 4*q2t(vb)*q1tk2t*P7*m1*q1q2

      AS1(2,vb) =
     & 4*k1(vb)*q2U*q1tq2t*q1tk2t*P3 - 4*k1(vb)*q2U*q1tq2t*q1tk2t
     & *P4 - 4*k1(vb)*q2U*q1tq2t*q1tk2t*P7 - 4*k1(vb)*q2U*q1tq2t*
     & q2tk2t*P3 + 4*k1(vb)*q2U*q1tq2t*q2tk2t*P4 - 4*k1(vb)*q2U*
     & q1tq2t*q2tk2t*P7 + 4*k1(vb)*q2U*q1tq2t*P3*m1*m2 - 4*k1(vb)*
     & q2U*q1tq2t*P3*q1k2 - 4*k1(vb)*q2U*q1tq2t*P3*q2k2 + 4*k1(vb)*
     & q2U*q1tq2t*P3*k1k2 - 4*k1(vb)*q2U*q1tq2t*P4*m1*m2 - 4*k1(vb)
     & *q2U*q1tq2t*P4*q1k2 - 4*k1(vb)*q2U*q1tq2t*P4*q2k2 - 4*k1(vb)
     & *q2U*q1tq2t*P4*k1k2 - 4*k1(vb)*q2U*q1tq2t*P7*q1k2 - 4*k1(vb)
     & *q2U*q1tq2t*P7*q2k2 + 4*k1(vb)*q2U*q1tk1t*q2tk2t*P3 + 4*
     & k1(vb)*q2U*q1tk1t*q2tk2t*P4 + 8*k1(vb)*q2U*q1tk1t*q2tk2t*
     & P7 - 4*k1(vb)*q2U*q1tk2t*q2tk1t*P3 - 4*k1(vb)*q2U*q1tk2t*
     & q2tk1t*P4 + 8*k1(vb)*q2U*q1tk2t*q2tk2t*P7 - 4*k1(vb)*q2U*
     & q1tk2t*t2*P3 + 4*k1(vb)*q2U*q1tk2t*t2*P4 + 4*k1(vb)*q2U*
     & q1tk2t*t2*P7 + 4*k1(vb)*q2U*q2tk2t*t1*P3 - 4*k1(vb)*q2U*
     & q2tk2t*t1*P4
      AS1(2,vb) = AS1(2,vb) + 4*k1(vb)*q2U*q2tk2t*t1*P7 + 4*k1(vb)*
     & q1tq2t*q1tU*P3*m1*m2 - 8*k1(vb)*q1tq2t*q1tU*P3*q1k2 - 4*
     & k1(vb)*q1tq2t*q1tU*P3*q2k2 + 4*k1(vb)*q1tq2t*q1tU*P3*k1k2 - 
     & 4*k1(vb)*q1tq2t*q1tU*P4*q2k2 - 4*k1(vb)*q1tq2t*q1tU*P7*q2k2
     &  + 4*k1(vb)*q1tq2t*q2tU*P3*q1k2 + 4*k1(vb)*q1tq2t*q2tU*P4*m1
     & *m2 + 4*k1(vb)*q1tq2t*q2tU*P4*q1k2 + 8*k1(vb)*q1tq2t*q2tU*P4
     & *q2k2 + 4*k1(vb)*q1tq2t*q2tU*P4*k1k2 + 4*k1(vb)*q1tq2t*q2tU*
     & P7*q1k2 - 8*k1(vb)*q1tk1t*q2tU*P3*m1*m2 + 4*k1(vb)*q1tk1t*
     & q2tU*P3*q1k2 + 8*k1(vb)*q1tk1t*q2tU*P3*q2k2 - 8*k1(vb)*
     & q1tk1t*q2tU*P3*k1k2 + 4*k1(vb)*q1tk1t*q2tU*P4*q2k2 + 8*
     & k1(vb)*q1tk1t*q2tU*P7*q2k2 - 4*k1(vb)*q1tk2t*q2tU*t1*P3 - 4*
     & k1(vb)*q1tk2t*q2tU*t2*P4 - 4*k1(vb)*q1tk2t*q2tU*P3*q1q2 + 4*
     & k1(vb)*q1tk2t*q2tU*P3*q1k1 - 4*k1(vb)*q1tk2t*q2tU*P4*q1q2 - 
     & 4*k1(vb)*q1tk2t*q2tU*P4*q2k1 - 4*k1(vb)*q1tk2t*q2tU*P7*q1q2
     &  + 4*k1(vb)*q1tU*q2tk1t*P3*q1k2 - 4*k1(vb)*q1tU*q2tk1t*P4*
     & q2k2
      AS1(2,vb) = AS1(2,vb) + 4*k1(vb)*q1tU*q2tk2t*t1*P3 + 4*k1(vb)
     & *q1tU*q2tk2t*t2*P4 + 4*k1(vb)*q1tU*q2tk2t*P3*q1q2 - 4*k1(vb)
     & *q1tU*q2tk2t*P3*q1k1 - 8*k1(vb)*q1tU*q2tk2t*P4*m1*m2 + 4*
     & k1(vb)*q1tU*q2tk2t*P4*q1q2 - 8*k1(vb)*q1tU*q2tk2t*P4*q1k2 + 
     & 4*k1(vb)*q1tU*q2tk2t*P4*q2k1 - 8*k1(vb)*q1tU*q2tk2t*P4*q2k2
     &  - 8*k1(vb)*q1tU*q2tk2t*P4*k1k2 + 4*k1(vb)*q1tU*q2tk2t*P7*
     & q1q2 - 8*k1(vb)*q1tU*q2tk2t*P7*q1k2 + 4*k1(vb)*q1tU*t2*P3*
     & q1k2 - 4*k1(vb)*q1tU*t2*P4*m1*m2 - 4*k1(vb)*q1tU*t2*P4*q1k2 - 
     & 4*k1(vb)*q1tU*t2*P4*k1k2 - 4*k1(vb)*q1tU*t2*P7*q1k2 - 4*k1(vb)
     & *q2tU*t1*P3*m1*m2 + 4*k1(vb)*q2tU*t1*P3*q2k2 - 4*k1(vb)*q2tU*
     & t1*P3*k1k2 - 4*k1(vb)*q2tU*t1*P4*q2k2 + 4*k1(vb)*q2tU*t1*P7*
     & q2k2 - 4*k2(vb)*q2U*q1tq2t*q1tk1t*P3 + 4*k2(vb)*q2U*q1tq2t*
     & q1tk1t*P4 + 4*k2(vb)*q2U*q1tq2t*q1tk1t*P7 + 4*k2(vb)*q2U*
     & q1tq2t*q2tk1t*P3 - 4*k2(vb)*q2U*q1tq2t*q2tk1t*P4 + 4*k2(vb)
     & *q2U*q1tq2t*q2tk1t*P7 + 4*k2(vb)*q2U*q1tq2t*P3*m1*m2 + 4*
     & k2(vb)*q2U*q1tq2t*P3*q1k1
      AS1(2,vb) = AS1(2,vb) + 4*k2(vb)*q2U*q1tq2t*P3*q2k1 + 4*
     & k2(vb)*q2U*q1tq2t*P3*k1k2 - 4*k2(vb)*q2U*q1tq2t*P4*m1*m2 + 4
     & *k2(vb)*q2U*q1tq2t*P4*q1k1 + 4*k2(vb)*q2U*q1tq2t*P4*q2k1 - 4
     & *k2(vb)*q2U*q1tq2t*P4*k1k2 + 4*k2(vb)*q2U*q1tq2t*P7*q1k1 + 4
     & *k2(vb)*q2U*q1tq2t*P7*q2k1 - 8*k2(vb)*q2U*q1tk1t*q2tk1t*P7
     &  - 4*k2(vb)*q2U*q1tk1t*q2tk2t*P3 - 4*k2(vb)*q2U*q1tk1t*
     & q2tk2t*P4 - 8*k2(vb)*q2U*q1tk1t*q2tk2t*P7 + 4*k2(vb)*q2U*
     & q1tk1t*t2*P3 - 4*k2(vb)*q2U*q1tk1t*t2*P4 - 4*k2(vb)*q2U*
     & q1tk1t*t2*P7 + 4*k2(vb)*q2U*q1tk2t*q2tk1t*P3 + 4*k2(vb)*q2U
     & *q1tk2t*q2tk1t*P4 - 4*k2(vb)*q2U*q2tk1t*t1*P3 + 4*k2(vb)*
     & q2U*q2tk1t*t1*P4 - 4*k2(vb)*q2U*q2tk1t*t1*P7 + 4*k2(vb)*
     & q1tq2t*q1tU*P3*m1*m2 + 8*k2(vb)*q1tq2t*q1tU*P3*q1k1 + 4*
     & k2(vb)*q1tq2t*q1tU*P3*q2k1 + 4*k2(vb)*q1tq2t*q1tU*P3*k1k2 + 
     & 4*k2(vb)*q1tq2t*q1tU*P4*q2k1 + 4*k2(vb)*q1tq2t*q1tU*P7*q2k1
     &  - 4*k2(vb)*q1tq2t*q2tU*P3*q1k1 + 4*k2(vb)*q1tq2t*q2tU*P4*m1
     & *m2
      AS1(2,vb) = AS1(2,vb) - 4*k2(vb)*q1tq2t*q2tU*P4*q1k1 - 8*
     & k2(vb)*q1tq2t*q2tU*P4*q2k1 + 4*k2(vb)*q1tq2t*q2tU*P4*k1k2 - 
     & 4*k2(vb)*q1tq2t*q2tU*P7*q1k1 + 4*k2(vb)*q1tk1t*q2tU*t1*P3 + 
     & 4*k2(vb)*q1tk1t*q2tU*t2*P4 - 8*k2(vb)*q1tk1t*q2tU*P3*m1*m2
     &  + 4*k2(vb)*q1tk1t*q2tU*P3*q1q2 - 8*k2(vb)*q1tk1t*q2tU*P3*
     & q1k1 + 4*k2(vb)*q1tk1t*q2tU*P3*q1k2 - 8*k2(vb)*q1tk1t*q2tU*
     & P3*q2k1 - 8*k2(vb)*q1tk1t*q2tU*P3*k1k2 + 4*k2(vb)*q1tk1t*
     & q2tU*P4*q1q2 - 4*k2(vb)*q1tk1t*q2tU*P4*q2k2 + 4*k2(vb)*
     & q1tk1t*q2tU*P7*q1q2 - 8*k2(vb)*q1tk1t*q2tU*P7*q2k1 - 4*
     & k2(vb)*q1tk2t*q2tU*P3*q1k1 + 4*k2(vb)*q1tk2t*q2tU*P4*q2k1 - 
     & 4*k2(vb)*q1tU*q2tk1t*t1*P3 - 4*k2(vb)*q1tU*q2tk1t*t2*P4 - 4*
     & k2(vb)*q1tU*q2tk1t*P3*q1q2 - 4*k2(vb)*q1tU*q2tk1t*P3*q1k2 -
     & 4*k2(vb)*q1tU*q2tk1t*P4*q1q2 + 4*k2(vb)*q1tU*q2tk1t*P4*q2k2
     &  - 4*k2(vb)*q1tU*q2tk1t*P7*q1q2 + 4*k2(vb)*q1tU*q2tk2t*P3*
     & q1k1 - 8*k2(vb)*q1tU*q2tk2t*P4*m1*m2 + 8*k2(vb)*q1tU*q2tk2t*
     & P4*q1k1
      AS1(2,vb) = AS1(2,vb) + 4*k2(vb)*q1tU*q2tk2t*P4*q2k1 - 8*
     & k2(vb)*q1tU*q2tk2t*P4*k1k2 + 8*k2(vb)*q1tU*q2tk2t*P7*q1k1 - 
     & 4*k2(vb)*q1tU*t2*P3*q1k1 - 4*k2(vb)*q1tU*t2*P4*m1*m2 + 4*
     & k2(vb)*q1tU*t2*P4*q1k1 - 4*k2(vb)*q1tU*t2*P4*k1k2 + 4*k2(vb)*
     & q1tU*t2*P7*q1k1 - 4*k2(vb)*q2tU*t1*P3*m1*m2 - 4*k2(vb)*q2tU*
     & t1*P3*q2k1 - 4*k2(vb)*q2tU*t1*P3*k1k2 + 4*k2(vb)*q2tU*t1*P4*
     & q2k1 - 4*k2(vb)*q2tU*t1*P7*q2k1 + 4*kw(vb)*q2U*q1tq2t*P3*m1*
     & m2 + 4*kw(vb)*q2U*q1tq2t*P3*k1k2 - 4*kw(vb)*q2U*q1tq2t*P4*m1
     & *m2 - 4*kw(vb)*q2U*q1tq2t*P4*k1k2 - 4*kw(vb)*q2U*q1tk1t*
     & q2tk2t*P3 + 4*kw(vb)*q2U*q1tk1t*q2tk2t*P4 + 4*kw(vb)*q2U*
     & q1tk2t*q2tk1t*P3 - 4*kw(vb)*q2U*q1tk2t*q2tk1t*P4 + 4*kw(vb)
     & *q1tq2t*q1tU*P3*m1*m2 + 4*kw(vb)*q1tq2t*q1tU*P3*k1k2 + 4*
     & kw(vb)*q1tq2t*q2tU*P4*m1*m2 + 4*kw(vb)*q1tq2t*q2tU*P4*k1k2
     &  - 8*kw(vb)*q1tk1t*q2tU*P3*m1*m2 + 4*kw(vb)*q1tk1t*q2tU*P3*
     & q1k2 - 8*kw(vb)*q1tk1t*q2tU*P3*k1k2 + 4*kw(vb)*q1tk1t*q2tU*
     & P4*q2k2
      AS1(2,vb) = AS1(2,vb) - 4*kw(vb)*q1tk2t*q2tU*P3*q1k1 - 4*
     & kw(vb)*q1tk2t*q2tU*P4*q2k1 - 4*kw(vb)*q1tU*q2tk1t*P3*q1k2 - 
     & 4*kw(vb)*q1tU*q2tk1t*P4*q2k2 + 4*kw(vb)*q1tU*q2tk2t*P3*q1k1
     &  - 8*kw(vb)*q1tU*q2tk2t*P4*m1*m2 + 4*kw(vb)*q1tU*q2tk2t*P4*
     & q2k1 - 8*kw(vb)*q1tU*q2tk2t*P4*k1k2 - 4*kw(vb)*q1tU*t2*P4*m1*
     & m2 - 4*kw(vb)*q1tU*t2*P4*k1k2 - 4*kw(vb)*q2tU*t1*P3*m1*m2 - 4*
     & kw(vb)*q2tU*t1*P3*k1k2 + 4*q1(vb)*q2U*q1tq2t*P7*m1*m2 + 4*
     & q1(vb)*q2U*q1tq2t*P7*k1k2 - 4*q1(vb)*q2U*q1tk1t*q2tk2t*P7
     &  + 4*q1(vb)*q2U*q1tk2t*q2tk1t*P7 + 4*q1(vb)*q1tq2t*q1tU*P3*
     & m1*m2 + 4*q1(vb)*q1tq2t*q1tU*P3*k1k2 - 4*q1(vb)*q1tq2t*q2tU*
     & P3*m1*m2 - 4*q1(vb)*q1tq2t*q2tU*P3*k1k2 - 4*q1(vb)*q1tq2t*
     & q2tU*P7*m1*m2 - 4*q1(vb)*q1tq2t*q2tU*P7*k1k2 + 8*q1(vb)*
     & q1tk1t*q2tU*P3*m1*m2 - 4*q1(vb)*q1tk1t*q2tU*P3*q1k2 - 4*
     & q1(vb)*q1tk1t*q2tU*P3*q2k2 + 8*q1(vb)*q1tk1t*q2tU*P3*k1k2 - 
     & 4*q1(vb)*q1tk1t*q2tU*P7*q2k2 + 4*q1(vb)*q1tk2t*q2tU*P3*q1k1
     &  + 4*q1(vb)*q1tk2t*q2tU*P3*q2k1
      AS1(2,vb) = AS1(2,vb) + 4*q1(vb)*q1tk2t*q2tU*P7*q2k1 - 8*
     & q1(vb)*q1tU*q2tk1t*P3*m1*m2 + 4*q1(vb)*q1tU*q2tk1t*P3*q1k2
     &  + 4*q1(vb)*q1tU*q2tk1t*P3*q2k2 - 8*q1(vb)*q1tU*q2tk1t*P3*
     & k1k2 + 4*q1(vb)*q1tU*q2tk1t*P7*q2k2 - 4*q1(vb)*q1tU*q2tk2t*
     & P3*q1k1 - 4*q1(vb)*q1tU*q2tk2t*P3*q2k1 + 8*q1(vb)*q1tU*
     & q2tk2t*P7*m1*m2 - 4*q1(vb)*q1tU*q2tk2t*P7*q2k1 + 8*q1(vb)*
     & q1tU*q2tk2t*P7*k1k2 - 4*q1(vb)*q1tU*t2*P3*m1*m2 - 4*q1(vb)*
     & q1tU*t2*P3*k1k2 + 4*q1(vb)*q1tU*t2*P7*m1*m2 + 4*q1(vb)*q1tU*
     & t2*P7*k1k2 + 4*q1(vb)*q2tU*t1*P3*m1*m2 + 4*q1(vb)*q2tU*t1*P3*
     & k1k2 - 4*q2(vb)*q2U*q1tq2t*P7*m1*m2 - 4*q2(vb)*q2U*q1tq2t*P7
     & *k1k2 + 4*q2(vb)*q2U*q1tk1t*q2tk2t*P7 - 4*q2(vb)*q2U*q1tk2t
     & *q2tk1t*P7 - 4*q2(vb)*q1tq2t*q1tU*P4*m1*m2 - 4*q2(vb)*q1tq2t
     & *q1tU*P4*k1k2 - 4*q2(vb)*q1tq2t*q1tU*P7*m1*m2 - 4*q2(vb)*
     & q1tq2t*q1tU*P7*k1k2 + 4*q2(vb)*q1tq2t*q2tU*P4*m1*m2 + 4*
     & q2(vb)*q1tq2t*q2tU*P4*k1k2 - 4*q2(vb)*q1tk1t*q2tU*P4*q1k2 - 
     & 4*q2(vb)*q1tk1t*q2tU*P4*q2k2
      AS1(2,vb) = AS1(2,vb) + 8*q2(vb)*q1tk1t*q2tU*P7*m1*m2 - 4*
     & q2(vb)*q1tk1t*q2tU*P7*q1k2 + 8*q2(vb)*q1tk1t*q2tU*P7*k1k2 - 
     & 8*q2(vb)*q1tk2t*q2tU*P4*m1*m2 + 4*q2(vb)*q1tk2t*q2tU*P4*q1k1
     &  + 4*q2(vb)*q1tk2t*q2tU*P4*q2k1 - 8*q2(vb)*q1tk2t*q2tU*P4*
     & k1k2 + 4*q2(vb)*q1tk2t*q2tU*P7*q1k1 + 4*q2(vb)*q1tU*q2tk1t*
     & P4*q1k2 + 4*q2(vb)*q1tU*q2tk1t*P4*q2k2 + 4*q2(vb)*q1tU*
     & q2tk1t*P7*q1k2 + 8*q2(vb)*q1tU*q2tk2t*P4*m1*m2 - 4*q2(vb)*
     & q1tU*q2tk2t*P4*q1k1 - 4*q2(vb)*q1tU*q2tk2t*P4*q2k1 + 8*
     & q2(vb)*q1tU*q2tk2t*P4*k1k2 - 4*q2(vb)*q1tU*q2tk2t*P7*q1k1 + 
     & 4*q2(vb)*q1tU*t2*P4*m1*m2 + 4*q2(vb)*q1tU*t2*P4*k1k2 - 4*
     & q2(vb)*q2tU*t1*P4*m1*m2 - 4*q2(vb)*q2tU*t1*P4*k1k2 + 4*q2(vb)*
     & q2tU*t1*P7*m1*m2 + 4*q2(vb)*q2tU*t1*P7*k1k2 + 4*q1t(vb)*q2U*
     & q1tq2t*P3*m1*m2 + 4*q1t(vb)*q2U*q1tq2t*P3*k1k2 + 4*q1t(vb)*
     & q2U*q1tq2t*P4*m1*m2 + 4*q1t(vb)*q2U*q1tq2t*P4*k1k2 - 4*
     & q1t(vb)*q2U*q1tq2t*P7*m1*m2 - 4*q1t(vb)*q2U*q1tq2t*P7*k1k2
     &  - 8*q1t(vb)*q2U*q2tk1t*P3*m1*m2
      AS1(2,vb) = AS1(2,vb) + 4*q1t(vb)*q2U*q2tk1t*P3*q1k2 + 4*
     & q1t(vb)*q2U*q2tk1t*P3*q2k2 - 8*q1t(vb)*q2U*q2tk1t*P3*k1k2 - 
     & 4*q1t(vb)*q2U*q2tk1t*P4*q1k2 - 4*q1t(vb)*q2U*q2tk1t*P4*q2k2
     &  + 4*q1t(vb)*q2U*q2tk1t*P7*q1k2 + 4*q1t(vb)*q2U*q2tk1t*P7*
     & q2k2 - 4*q1t(vb)*q2U*q2tk2t*P3*q1k1 - 4*q1t(vb)*q2U*q2tk2t*
     & P3*q2k1 - 8*q1t(vb)*q2U*q2tk2t*P4*m1*m2 + 4*q1t(vb)*q2U*
     & q2tk2t*P4*q1k1 + 4*q1t(vb)*q2U*q2tk2t*P4*q2k1 - 8*q1t(vb)*
     & q2U*q2tk2t*P4*k1k2 + 8*q1t(vb)*q2U*q2tk2t*P7*m1*m2 - 4*
     & q1t(vb)*q2U*q2tk2t*P7*q1k1 - 4*q1t(vb)*q2U*q2tk2t*P7*q2k1 + 
     & 8*q1t(vb)*q2U*q2tk2t*P7*k1k2 - 4*q1t(vb)*q2U*t2*P3*m1*m2 - 4*
     & q1t(vb)*q2U*t2*P3*k1k2 - 4*q1t(vb)*q2U*t2*P4*m1*m2 - 4*q1t(vb)
     & *q2U*t2*P4*k1k2 + 4*q1t(vb)*q2U*t2*P7*m1*m2 + 4*q1t(vb)*q2U*
     & t2*P7*k1k2 - 4*q1t(vb)*q2tU*t1*P3*m1*m2 - 4*q1t(vb)*q2tU*t1*P3
     & *k1k2 + 4*q1t(vb)*q2tU*t2*P4*m1*m2 + 4*q1t(vb)*q2tU*t2*P4*k1k2
     &  - 4*q1t(vb)*q2tU*P3*m1*m2*q1q2 + 8*q1t(vb)*q2tU*P3*m1*m2*q1k1
     &  - 4*q1t(vb)*q2tU*P3*q1q2*k1k2
      AS1(2,vb) = AS1(2,vb) - 4*q1t(vb)*q2tU*P3*q1k1*q2k2 + 8*
     & q1t(vb)*q2tU*P3*q1k1*k1k2 + 4*q1t(vb)*q2tU*P3*q1k2*q2k1 + 4*
     & q1t(vb)*q2tU*P4*m1*m2*q1q2 - 8*q1t(vb)*q2tU*P4*m1*m2*q2k2 + 4*
     & q1t(vb)*q2tU*P4*q1q2*k1k2 + 4*q1t(vb)*q2tU*P4*q1k1*q2k2 - 4*
     & q1t(vb)*q2tU*P4*q1k2*q2k1 - 8*q1t(vb)*q2tU*P4*q2k2*k1k2 - 4*
     & q1t(vb)*q2tU*P7*m1*m2*q1q2 - 4*q1t(vb)*q2tU*P7*q1q2*k1k2 - 4*
     & q1t(vb)*q2tU*P7*q1k1*q2k2 + 4*q1t(vb)*q2tU*P7*q1k2*q2k1 - 4*
     & q2t(vb)*q2U*q1tq2t*P3*m1*m2 - 4*q2t(vb)*q2U*q1tq2t*P3*k1k2
     &  - 4*q2t(vb)*q2U*q1tq2t*P4*m1*m2 - 4*q2t(vb)*q2U*q1tq2t*P4*
     & k1k2 + 4*q2t(vb)*q2U*q1tq2t*P7*m1*m2 + 4*q2t(vb)*q2U*q1tq2t*
     & P7*k1k2 + 8*q2t(vb)*q2U*q1tk1t*P3*m1*m2 - 4*q2t(vb)*q2U*
     & q1tk1t*P3*q1k2 - 4*q2t(vb)*q2U*q1tk1t*P3*q2k2 + 8*q2t(vb)*
     & q2U*q1tk1t*P3*k1k2 + 4*q2t(vb)*q2U*q1tk1t*P4*q1k2 + 4*
     & q2t(vb)*q2U*q1tk1t*P4*q2k2 - 8*q2t(vb)*q2U*q1tk1t*P7*m1*m2
     &  + 4*q2t(vb)*q2U*q1tk1t*P7*q1k2 + 4*q2t(vb)*q2U*q1tk1t*P7*
     & q2k2
      AS1(2,vb) = AS1(2,vb) - 8*q2t(vb)*q2U*q1tk1t*P7*k1k2 + 4*
     & q2t(vb)*q2U*q1tk2t*P3*q1k1 + 4*q2t(vb)*q2U*q1tk2t*P3*q2k1 + 
     & 8*q2t(vb)*q2U*q1tk2t*P4*m1*m2 - 4*q2t(vb)*q2U*q1tk2t*P4*q1k1
     &  - 4*q2t(vb)*q2U*q1tk2t*P4*q2k1 + 8*q2t(vb)*q2U*q1tk2t*P4*
     & k1k2 - 4*q2t(vb)*q2U*q1tk2t*P7*q1k1 - 4*q2t(vb)*q2U*q1tk2t*
     & P7*q2k1 + 4*q2t(vb)*q2U*t1*P3*m1*m2 + 4*q2t(vb)*q2U*t1*P3*k1k2
     &  + 4*q2t(vb)*q2U*t1*P4*m1*m2 + 4*q2t(vb)*q2U*t1*P4*k1k2 - 4*
     & q2t(vb)*q2U*t1*P7*m1*m2 - 4*q2t(vb)*q2U*t1*P7*k1k2 + 4*q2t(vb)
     & *q1tU*t1*P3*m1*m2 + 4*q2t(vb)*q1tU*t1*P3*k1k2 - 4*q2t(vb)*
     & q1tU*t2*P4*m1*m2 - 4*q2t(vb)*q1tU*t2*P4*k1k2 + 4*q2t(vb)*q1tU
     & *P3*m1*m2*q1q2 - 8*q2t(vb)*q1tU*P3*m1*m2*q1k1 + 4*q2t(vb)*q1tU
     & *P3*q1q2*k1k2 + 4*q2t(vb)*q1tU*P3*q1k1*q2k2 - 8*q2t(vb)*q1tU*
     & P3*q1k1*k1k2 - 4*q2t(vb)*q1tU*P3*q1k2*q2k1 - 4*q2t(vb)*q1tU*P4
     & *m1*m2*q1q2 + 8*q2t(vb)*q1tU*P4*m1*m2*q2k2 - 4*q2t(vb)*q1tU*P4
     & *q1q2*k1k2 - 4*q2t(vb)*q1tU*P4*q1k1*q2k2 + 4*q2t(vb)*q1tU*P4*
     & q1k2*q2k1
      AS1(2,vb) = AS1(2,vb) + 8*q2t(vb)*q1tU*P4*q2k2*k1k2 - 4*
     & q2t(vb)*q1tU*P7*m1*m2*q1q2 - 4*q2t(vb)*q1tU*P7*q1q2*k1k2 - 4*
     & q2t(vb)*q1tU*P7*q1k1*q2k2 + 4*q2t(vb)*q1tU*P7*q1k2*q2k1 - 4*
     & U(vb)*t1*t2*P3*m1*m2 - 4*U(vb)*t1*t2*P3*k1k2 - 4*U(vb)*t1*t2*P4*
     & m1*m2 - 4*U(vb)*t1*t2*P4*k1k2 + 4*U(vb)*t1*t2*P7*m1*m2 + 4*U(vb)
     & *t1*t2*P7*k1k2 + 16*U(vb)*q1tq2t*q1tk1t*P3*m1*m2 - 8*U(vb)*
     & q1tq2t*q1tk1t*P3*q1k2 - 4*U(vb)*q1tq2t*q1tk1t*P3*q2k2 + 16*
     & U(vb)*q1tq2t*q1tk1t*P3*k1k2 - 4*U(vb)*q1tq2t*q1tk1t*P4*q2k2
     &  - 4*U(vb)*q1tq2t*q1tk1t*P7*q2k2 + 8*U(vb)*q1tq2t*q1tk2t*P3*
     & q1k1 + 4*U(vb)*q1tq2t*q1tk2t*P3*q2k1 + 4*U(vb)*q1tq2t*q1tk2t
     & *P4*q2k1 + 4*U(vb)*q1tq2t*q1tk2t*P7*q2k1 + 4*U(vb)*q1tq2t*
     & q2tk1t*P3*q1k2 + 4*U(vb)*q1tq2t*q2tk1t*P4*q1k2 + 8*U(vb)*
     & q1tq2t*q2tk1t*P4*q2k2 + 4*U(vb)*q1tq2t*q2tk1t*P7*q1k2 - 4*
     & U(vb)*q1tq2t*q2tk2t*P3*q1k1 + 16*U(vb)*q1tq2t*q2tk2t*P4*m1*
     & m2 - 4*U(vb)*q1tq2t*q2tk2t*P4*q1k1 - 8*U(vb)*q1tq2t*q2tk2t*
     & P4*q2k1
      AS1(2,vb) = AS1(2,vb) + 16*U(vb)*q1tq2t*q2tk2t*P4*k1k2 - 4*
     & U(vb)*q1tq2t*q2tk2t*P7*q1k1 + 4*U(vb)*q1tq2t*t1*P3*m1*m2 + 4*
     & U(vb)*q1tq2t*t1*P3*k1k2 + 4*U(vb)*q1tq2t*t2*P4*m1*m2 + 4*U(vb)
     & *q1tq2t*t2*P4*k1k2 - 4*U(vb)*q1tq2t*P3*m1*m2*q1q2 + 8*U(vb)*
     & q1tq2t*P3*m1*m2*q1k1 - 4*U(vb)*q1tq2t*P3*q1q2*k1k2 - 4*U(vb)*
     & q1tq2t*P3*q1k1*q2k2 + 8*U(vb)*q1tq2t*P3*q1k1*k1k2 + 4*U(vb)*
     & q1tq2t*P3*q1k2*q2k1 - 4*U(vb)*q1tq2t*P4*m1*m2*q1q2 + 8*U(vb)*
     & q1tq2t*P4*m1*m2*q2k2 - 4*U(vb)*q1tq2t*P4*q1q2*k1k2 - 4*U(vb)*
     & q1tq2t*P4*q1k1*q2k2 + 4*U(vb)*q1tq2t*P4*q1k2*q2k1 + 8*U(vb)*
     & q1tq2t*P4*q2k2*k1k2 - 4*U(vb)*q1tq2t*P7*m1*m2*q1q2 - 4*U(vb)*
     & q1tq2t*P7*q1q2*k1k2 - 4*U(vb)*q1tq2t*P7*q1k1*q2k2 + 4*U(vb)*
     & q1tq2t*P7*q1k2*q2k1 - 4*U(vb)*q1tq2t**2*P3*m1*m2 - 4*U(vb)*
     & q1tq2t**2*P3*k1k2 - 4*U(vb)*q1tq2t**2*P4*m1*m2 - 4*U(vb)*
     & q1tq2t**2*P4*k1k2 - 4*U(vb)*q1tq2t**2*P7*m1*m2 - 4*U(vb)*
     & q1tq2t**2*P7*k1k2 - 16*U(vb)*q1tk1t*q2tk1t*P3*m1*m2 + 8*U(vb)
     & *q1tk1t*q2tk1t*P3*q1k2
      AS1(2,vb) = AS1(2,vb) + 8*U(vb)*q1tk1t*q2tk1t*P3*q2k2 - 16*
     & U(vb)*q1tk1t*q2tk1t*P3*k1k2 + 8*U(vb)*q1tk1t*q2tk1t*P7*q2k2
     &  + 4*U(vb)*q1tk1t*q2tk2t*t1*P3 + 4*U(vb)*q1tk1t*q2tk2t*t2*P4
     &  + 4*U(vb)*q1tk1t*q2tk2t*P3*q1q2 - 8*U(vb)*q1tk1t*q2tk2t*P3*
     & q1k1 - 8*U(vb)*q1tk1t*q2tk2t*P3*q2k1 + 4*U(vb)*q1tk1t*q2tk2t
     & *P4*q1q2 - 8*U(vb)*q1tk1t*q2tk2t*P4*q1k2 - 8*U(vb)*q1tk1t*
     & q2tk2t*P4*q2k2 + 16*U(vb)*q1tk1t*q2tk2t*P7*m1*m2 + 4*U(vb)*
     & q1tk1t*q2tk2t*P7*q1q2 - 8*U(vb)*q1tk1t*q2tk2t*P7*q1k2 - 8*
     & U(vb)*q1tk1t*q2tk2t*P7*q2k1 + 16*U(vb)*q1tk1t*q2tk2t*P7*k1k2
     &  - 8*U(vb)*q1tk1t*t2*P3*m1*m2 + 4*U(vb)*q1tk1t*t2*P3*q1k2 - 8*
     & U(vb)*q1tk1t*t2*P3*k1k2 - 4*U(vb)*q1tk1t*t2*P4*q1k2 + 8*U(vb)*
     & q1tk1t*t2*P7*m1*m2 - 4*U(vb)*q1tk1t*t2*P7*q1k2 + 8*U(vb)*
     & q1tk1t*t2*P7*k1k2 - 4*U(vb)*q1tk2t*q2tk1t*t1*P3 - 4*U(vb)*
     & q1tk2t*q2tk1t*t2*P4 - 4*U(vb)*q1tk2t*q2tk1t*P3*q1q2 - 4*
     & U(vb)*q1tk2t*q2tk1t*P4*q1q2 - 4*U(vb)*q1tk2t*q2tk1t*P7*q1q2
     &  - 16*U(vb)*q1tk2t*q2tk2t*P4*m1*m2
      AS1(2,vb) = AS1(2,vb) + 8*U(vb)*q1tk2t*q2tk2t*P4*q1k1 + 8*
     & U(vb)*q1tk2t*q2tk2t*P4*q2k1 - 16*U(vb)*q1tk2t*q2tk2t*P4*k1k2
     &  + 8*U(vb)*q1tk2t*q2tk2t*P7*q1k1 - 4*U(vb)*q1tk2t*t2*P3*q1k1
     &  - 8*U(vb)*q1tk2t*t2*P4*m1*m2 + 4*U(vb)*q1tk2t*t2*P4*q1k1 - 8*
     & U(vb)*q1tk2t*t2*P4*k1k2 + 4*U(vb)*q1tk2t*t2*P7*q1k1 - 8*U(vb)*
     & q2tk1t*t1*P3*m1*m2 + 4*U(vb)*q2tk1t*t1*P3*q2k2 - 8*U(vb)*
     & q2tk1t*t1*P3*k1k2 - 4*U(vb)*q2tk1t*t1*P4*q2k2 + 4*U(vb)*
     & q2tk1t*t1*P7*q2k2 - 4*U(vb)*q2tk2t*t1*P3*q2k1 - 8*U(vb)*
     & q2tk2t*t1*P4*m1*m2 + 4*U(vb)*q2tk2t*t1*P4*q2k1 - 8*U(vb)*
     & q2tk2t*t1*P4*k1k2 + 8*U(vb)*q2tk2t*t1*P7*m1*m2 - 4*U(vb)*
     & q2tk2t*t1*P7*q2k1 + 8*U(vb)*q2tk2t*t1*P7*k1k2

      AS1(3,vb) =
     &  - 8*k1(vb)*q1L*q1tq2t*q1tk2t*P3 + 4*k1(vb)*q1L*q1tq2t*
     & q2tk2t*P3 + 4*k1(vb)*q1L*q1tq2t*q2tk2t*P4 + 4*k1(vb)*q1L*
     & q1tq2t*q2tk2t*P7 - 4*k1(vb)*q1L*q1tq2t*P3*m1*m2 + 4*k1(vb)*
     & q1L*q1tq2t*P3*q2k2 - 4*k1(vb)*q1L*q1tq2t*P3*k1k2 + 4*k1(vb)*
     & q1L*q1tq2t*P4*q2k2 + 4*k1(vb)*q1L*q1tq2t*P7*q2k2 + 4*k1(vb)*
     & q1L*q1tk1t*q2tk2t*P3 + 4*k1(vb)*q1L*q1tk2t*q2tk1t*P3 - 8*
     & k1(vb)*q1L*q1tk2t*q2tk2t*P4 - 8*k1(vb)*q1L*q1tk2t*q2tk2t*
     & P7 + 4*k1(vb)*q1L*q1tk2t*t2*P3 - 4*k1(vb)*q1L*q1tk2t*t2*P4
     &  - 4*k1(vb)*q1L*q1tk2t*t2*P7 - 4*k1(vb)*q2L*q1tq2t*q1tk2t*
     & P3 - 4*k1(vb)*q2L*q1tq2t*q1tk2t*P4 - 4*k1(vb)*q2L*q1tq2t*
     & q1tk2t*P7 + 8*k1(vb)*q2L*q1tq2t*q2tk2t*P4 - 4*k1(vb)*q2L*
     & q1tq2t*P3*q1k2 - 4*k1(vb)*q2L*q1tq2t*P4*m1*m2 - 4*k1(vb)*q2L
     & *q1tq2t*P4*q1k2 - 4*k1(vb)*q2L*q1tq2t*P4*k1k2 - 4*k1(vb)*q2L
     & *q1tq2t*P7*q1k2 + 8*k1(vb)*q2L*q1tk1t*q2tk2t*P3 + 4*k1(vb)*
     & q2L*q1tk1t*q2tk2t*P4
      AS1(3,vb) = AS1(3,vb) + 8*k1(vb)*q2L*q1tk1t*q2tk2t*P7 - 4*
     & k1(vb)*q2L*q1tk2t*q2tk1t*P4 - 8*k1(vb)*q2L*q1tk2t*q2tk2t*
     & P4 + 4*k1(vb)*q2L*q2tk2t*t1*P3 - 4*k1(vb)*q2L*q2tk2t*t1*P4
     &  + 4*k1(vb)*q2L*q2tk2t*t1*P7 + 4*k1(vb)*q1tq2t*q1tL*P3*m1*m2
     &  - 8*k1(vb)*q1tq2t*q1tL*P3*q1k2 - 4*k1(vb)*q1tq2t*q1tL*P3*
     & q2k2 + 4*k1(vb)*q1tq2t*q1tL*P3*k1k2 - 4*k1(vb)*q1tq2t*q1tL*
     & P4*q2k2 - 4*k1(vb)*q1tq2t*q1tL*P7*q2k2 + 4*k1(vb)*q1tq2t*
     & q2tL*P3*q1k2 + 4*k1(vb)*q1tq2t*q2tL*P4*m1*m2 + 4*k1(vb)*
     & q1tq2t*q2tL*P4*q1k2 + 8*k1(vb)*q1tq2t*q2tL*P4*q2k2 + 4*
     & k1(vb)*q1tq2t*q2tL*P4*k1k2 + 4*k1(vb)*q1tq2t*q2tL*P7*q1k2 - 
     & 8*k1(vb)*q1tk1t*q2tL*P3*m1*m2 + 4*k1(vb)*q1tk1t*q2tL*P3*q1k2
     &  + 8*k1(vb)*q1tk1t*q2tL*P3*q2k2 - 8*k1(vb)*q1tk1t*q2tL*P3*
     & k1k2 + 4*k1(vb)*q1tk1t*q2tL*P4*q2k2 + 8*k1(vb)*q1tk1t*q2tL*
     & P7*q2k2 - 4*k1(vb)*q1tk2t*q2tL*t1*P3 - 4*k1(vb)*q1tk2t*q2tL*
     & t2*P4 - 4*k1(vb)*q1tk2t*q2tL*P3*q1q2 + 4*k1(vb)*q1tk2t*q2tL*
     & P3*q1k1
      AS1(3,vb) = AS1(3,vb) - 4*k1(vb)*q1tk2t*q2tL*P4*q1q2 - 4*
     & k1(vb)*q1tk2t*q2tL*P4*q2k1 - 4*k1(vb)*q1tk2t*q2tL*P7*q1q2 + 
     & 4*k1(vb)*q1tL*q2tk1t*P3*q1k2 - 4*k1(vb)*q1tL*q2tk1t*P4*q2k2
     &  + 4*k1(vb)*q1tL*q2tk2t*t1*P3 + 4*k1(vb)*q1tL*q2tk2t*t2*P4
     &  + 4*k1(vb)*q1tL*q2tk2t*P3*q1q2 - 4*k1(vb)*q1tL*q2tk2t*P3*
     & q1k1 - 8*k1(vb)*q1tL*q2tk2t*P4*m1*m2 + 4*k1(vb)*q1tL*q2tk2t*
     & P4*q1q2 - 8*k1(vb)*q1tL*q2tk2t*P4*q1k2 + 4*k1(vb)*q1tL*
     & q2tk2t*P4*q2k1 - 8*k1(vb)*q1tL*q2tk2t*P4*q2k2 - 8*k1(vb)*
     & q1tL*q2tk2t*P4*k1k2 + 4*k1(vb)*q1tL*q2tk2t*P7*q1q2 - 8*
     & k1(vb)*q1tL*q2tk2t*P7*q1k2 + 4*k1(vb)*q1tL*t2*P3*q1k2 - 4*
     & k1(vb)*q1tL*t2*P4*m1*m2 - 4*k1(vb)*q1tL*t2*P4*q1k2 - 4*k1(vb)*
     & q1tL*t2*P4*k1k2 - 4*k1(vb)*q1tL*t2*P7*q1k2 - 4*k1(vb)*q2tL*t1
     & *P3*m1*m2 + 4*k1(vb)*q2tL*t1*P3*q2k2 - 4*k1(vb)*q2tL*t1*P3*
     & k1k2 - 4*k1(vb)*q2tL*t1*P4*q2k2 + 4*k1(vb)*q2tL*t1*P7*q2k2 + 8
     & *k2(vb)*q1L*q1tq2t*q1tk1t*P3 - 4*k2(vb)*q1L*q1tq2t*q2tk1t*
     & P3
      AS1(3,vb) = AS1(3,vb) - 4*k2(vb)*q1L*q1tq2t*q2tk1t*P4 - 4*
     & k2(vb)*q1L*q1tq2t*q2tk1t*P7 - 4*k2(vb)*q1L*q1tq2t*P3*m1*m2
     &  - 4*k2(vb)*q1L*q1tq2t*P3*q2k1 - 4*k2(vb)*q1L*q1tq2t*P3*k1k2
     &  - 4*k2(vb)*q1L*q1tq2t*P4*q2k1 - 4*k2(vb)*q1L*q1tq2t*P7*q2k1
     &  - 8*k2(vb)*q1L*q1tk1t*q2tk1t*P3 + 4*k2(vb)*q1L*q1tk1t*
     & q2tk2t*P3 + 8*k2(vb)*q1L*q1tk1t*q2tk2t*P4 + 8*k2(vb)*q1L*
     & q1tk1t*q2tk2t*P7 - 4*k2(vb)*q1L*q1tk1t*t2*P3 + 4*k2(vb)*q1L
     & *q1tk1t*t2*P4 + 4*k2(vb)*q1L*q1tk1t*t2*P7 - 4*k2(vb)*q1L*
     & q1tk2t*q2tk1t*P3 + 4*k2(vb)*q2L*q1tq2t*q1tk1t*P3 + 4*k2(vb)
     & *q2L*q1tq2t*q1tk1t*P4 + 4*k2(vb)*q2L*q1tq2t*q1tk1t*P7 - 8*
     & k2(vb)*q2L*q1tq2t*q2tk1t*P4 + 4*k2(vb)*q2L*q1tq2t*P3*q1k1
     &  - 4*k2(vb)*q2L*q1tq2t*P4*m1*m2 + 4*k2(vb)*q2L*q1tq2t*P4*
     & q1k1 - 4*k2(vb)*q2L*q1tq2t*P4*k1k2 + 4*k2(vb)*q2L*q1tq2t*P7*
     & q1k1 - 8*k2(vb)*q2L*q1tk1t*q2tk1t*P3 - 8*k2(vb)*q2L*q1tk1t*
     & q2tk1t*P7 + 4*k2(vb)*q2L*q1tk1t*q2tk2t*P4 + 4*k2(vb)*q2L*
     & q1tk2t*q2tk1t*P4
      AS1(3,vb) = AS1(3,vb) - 4*k2(vb)*q2L*q2tk1t*t1*P3 + 4*k2(vb)*
     & q2L*q2tk1t*t1*P4 - 4*k2(vb)*q2L*q2tk1t*t1*P7 + 4*k2(vb)*
     & q1tq2t*q1tL*P3*m1*m2 + 8*k2(vb)*q1tq2t*q1tL*P3*q1k1 + 4*
     & k2(vb)*q1tq2t*q1tL*P3*q2k1 + 4*k2(vb)*q1tq2t*q1tL*P3*k1k2 + 
     & 4*k2(vb)*q1tq2t*q1tL*P4*q2k1 + 4*k2(vb)*q1tq2t*q1tL*P7*q2k1
     &  - 4*k2(vb)*q1tq2t*q2tL*P3*q1k1 + 4*k2(vb)*q1tq2t*q2tL*P4*m1
     & *m2 - 4*k2(vb)*q1tq2t*q2tL*P4*q1k1 - 8*k2(vb)*q1tq2t*q2tL*P4
     & *q2k1 + 4*k2(vb)*q1tq2t*q2tL*P4*k1k2 - 4*k2(vb)*q1tq2t*q2tL*
     & P7*q1k1 + 4*k2(vb)*q1tk1t*q2tL*t1*P3 + 4*k2(vb)*q1tk1t*q2tL*
     & t2*P4 - 8*k2(vb)*q1tk1t*q2tL*P3*m1*m2 + 4*k2(vb)*q1tk1t*q2tL
     & *P3*q1q2 - 8*k2(vb)*q1tk1t*q2tL*P3*q1k1 + 4*k2(vb)*q1tk1t*
     & q2tL*P3*q1k2 - 8*k2(vb)*q1tk1t*q2tL*P3*q2k1 - 8*k2(vb)*
     & q1tk1t*q2tL*P3*k1k2 + 4*k2(vb)*q1tk1t*q2tL*P4*q1q2 - 4*
     & k2(vb)*q1tk1t*q2tL*P4*q2k2 + 4*k2(vb)*q1tk1t*q2tL*P7*q1q2 - 
     & 8*k2(vb)*q1tk1t*q2tL*P7*q2k1 - 4*k2(vb)*q1tk2t*q2tL*P3*q1k1
     &  + 4*k2(vb)*q1tk2t*q2tL*P4*q2k1
      AS1(3,vb) = AS1(3,vb) - 4*k2(vb)*q1tL*q2tk1t*t1*P3 - 4*k2(vb)
     & *q1tL*q2tk1t*t2*P4 - 4*k2(vb)*q1tL*q2tk1t*P3*q1q2 - 4*k2(vb)
     & *q1tL*q2tk1t*P3*q1k2 - 4*k2(vb)*q1tL*q2tk1t*P4*q1q2 + 4*
     & k2(vb)*q1tL*q2tk1t*P4*q2k2 - 4*k2(vb)*q1tL*q2tk1t*P7*q1q2 + 
     & 4*k2(vb)*q1tL*q2tk2t*P3*q1k1 - 8*k2(vb)*q1tL*q2tk2t*P4*m1*m2
     &  + 8*k2(vb)*q1tL*q2tk2t*P4*q1k1 + 4*k2(vb)*q1tL*q2tk2t*P4*
     & q2k1 - 8*k2(vb)*q1tL*q2tk2t*P4*k1k2 + 8*k2(vb)*q1tL*q2tk2t*
     & P7*q1k1 - 4*k2(vb)*q1tL*t2*P3*q1k1 - 4*k2(vb)*q1tL*t2*P4*m1*m2
     &  + 4*k2(vb)*q1tL*t2*P4*q1k1 - 4*k2(vb)*q1tL*t2*P4*k1k2 + 4*
     & k2(vb)*q1tL*t2*P7*q1k1 - 4*k2(vb)*q2tL*t1*P3*m1*m2 - 4*k2(vb)*
     & q2tL*t1*P3*q2k1 - 4*k2(vb)*q2tL*t1*P3*k1k2 + 4*k2(vb)*q2tL*t1
     & *P4*q2k1 - 4*k2(vb)*q2tL*t1*P7*q2k1 - 4*kw(vb)*q1L*q1tq2t*P3*
     & m1*m2 - 4*kw(vb)*q1L*q1tq2t*P3*k1k2 + 4*kw(vb)*q1L*q1tk1t*
     & q2tk2t*P3 - 4*kw(vb)*q1L*q1tk2t*q2tk1t*P3 - 4*kw(vb)*q2L*
     & q1tq2t*P4*m1*m2 - 4*kw(vb)*q2L*q1tq2t*P4*k1k2 + 4*kw(vb)*q2L
     & *q1tk1t*q2tk2t*P4
      AS1(3,vb) = AS1(3,vb) - 4*kw(vb)*q2L*q1tk2t*q2tk1t*P4 + 4*
     & kw(vb)*q1tq2t*q1tL*P3*m1*m2 + 4*kw(vb)*q1tq2t*q1tL*P3*k1k2
     &  + 4*kw(vb)*q1tq2t*q2tL*P4*m1*m2 + 4*kw(vb)*q1tq2t*q2tL*P4*
     & k1k2 - 8*kw(vb)*q1tk1t*q2tL*P3*m1*m2 + 4*kw(vb)*q1tk1t*q2tL*
     & P3*q1k2 - 8*kw(vb)*q1tk1t*q2tL*P3*k1k2 + 4*kw(vb)*q1tk1t*
     & q2tL*P4*q2k2 - 4*kw(vb)*q1tk2t*q2tL*P3*q1k1 - 4*kw(vb)*
     & q1tk2t*q2tL*P4*q2k1 - 4*kw(vb)*q1tL*q2tk1t*P3*q1k2 - 4*
     & kw(vb)*q1tL*q2tk1t*P4*q2k2 + 4*kw(vb)*q1tL*q2tk2t*P3*q1k1 - 
     & 8*kw(vb)*q1tL*q2tk2t*P4*m1*m2 + 4*kw(vb)*q1tL*q2tk2t*P4*q2k1
     &  - 8*kw(vb)*q1tL*q2tk2t*P4*k1k2 - 4*kw(vb)*q1tL*t2*P4*m1*m2
     &  - 4*kw(vb)*q1tL*t2*P4*k1k2 - 4*kw(vb)*q2tL*t1*P3*m1*m2 - 4*
     & kw(vb)*q2tL*t1*P3*k1k2 + 4*q1(vb)*q1L*q1tq2t*P3*m1*m2 + 4*
     & q1(vb)*q1L*q1tq2t*P3*k1k2 - 4*q1(vb)*q1L*q1tk1t*q2tk2t*P3
     &  + 4*q1(vb)*q1L*q1tk2t*q2tk1t*P3 + 4*q1(vb)*q2L*q1tq2t*P3*
     & m1*m2 + 4*q1(vb)*q2L*q1tq2t*P3*k1k2 + 4*q1(vb)*q2L*q1tq2t*P7
     & *m1*m2
      AS1(3,vb) = AS1(3,vb) + 4*q1(vb)*q2L*q1tq2t*P7*k1k2 - 4*
     & q1(vb)*q2L*q1tk1t*q2tk2t*P3 - 4*q1(vb)*q2L*q1tk1t*q2tk2t*
     & P7 + 4*q1(vb)*q2L*q1tk2t*q2tk1t*P3 + 4*q1(vb)*q2L*q1tk2t*
     & q2tk1t*P7 + 4*q1(vb)*q1tq2t*q1tL*P3*m1*m2 + 4*q1(vb)*q1tq2t*
     & q1tL*P3*k1k2 - 4*q1(vb)*q1tq2t*q2tL*P3*m1*m2 - 4*q1(vb)*
     & q1tq2t*q2tL*P3*k1k2 - 4*q1(vb)*q1tq2t*q2tL*P7*m1*m2 - 4*
     & q1(vb)*q1tq2t*q2tL*P7*k1k2 + 8*q1(vb)*q1tk1t*q2tL*P3*m1*m2
     &  - 4*q1(vb)*q1tk1t*q2tL*P3*q1k2 - 4*q1(vb)*q1tk1t*q2tL*P3*
     & q2k2 + 8*q1(vb)*q1tk1t*q2tL*P3*k1k2 - 4*q1(vb)*q1tk1t*q2tL*
     & P7*q2k2 + 4*q1(vb)*q1tk2t*q2tL*P3*q1k1 + 4*q1(vb)*q1tk2t*
     & q2tL*P3*q2k1 + 4*q1(vb)*q1tk2t*q2tL*P7*q2k1 - 8*q1(vb)*q1tL*
     & q2tk1t*P3*m1*m2 + 4*q1(vb)*q1tL*q2tk1t*P3*q1k2 + 4*q1(vb)*
     & q1tL*q2tk1t*P3*q2k2 - 8*q1(vb)*q1tL*q2tk1t*P3*k1k2 + 4*
     & q1(vb)*q1tL*q2tk1t*P7*q2k2 - 4*q1(vb)*q1tL*q2tk2t*P3*q1k1 - 
     & 4*q1(vb)*q1tL*q2tk2t*P3*q2k1 + 8*q1(vb)*q1tL*q2tk2t*P7*m1*m2
     &  - 4*q1(vb)*q1tL*q2tk2t*P7*q2k1
      AS1(3,vb) = AS1(3,vb) + 8*q1(vb)*q1tL*q2tk2t*P7*k1k2 - 4*
     & q1(vb)*q1tL*t2*P3*m1*m2 - 4*q1(vb)*q1tL*t2*P3*k1k2 + 4*q1(vb)*
     & q1tL*t2*P7*m1*m2 + 4*q1(vb)*q1tL*t2*P7*k1k2 + 4*q1(vb)*q2tL*
     & t1*P3*m1*m2 + 4*q1(vb)*q2tL*t1*P3*k1k2 + 4*q2(vb)*q1L*q1tq2t*
     & P4*m1*m2 + 4*q2(vb)*q1L*q1tq2t*P4*k1k2 + 4*q2(vb)*q1L*q1tq2t
     & *P7*m1*m2 + 4*q2(vb)*q1L*q1tq2t*P7*k1k2 - 4*q2(vb)*q1L*
     & q1tk1t*q2tk2t*P4 - 4*q2(vb)*q1L*q1tk1t*q2tk2t*P7 + 4*q2(vb)
     & *q1L*q1tk2t*q2tk1t*P4 + 4*q2(vb)*q1L*q1tk2t*q2tk1t*P7 + 4*
     & q2(vb)*q2L*q1tq2t*P4*m1*m2 + 4*q2(vb)*q2L*q1tq2t*P4*k1k2 - 4
     & *q2(vb)*q2L*q1tk1t*q2tk2t*P4 + 4*q2(vb)*q2L*q1tk2t*q2tk1t*
     & P4 - 4*q2(vb)*q1tq2t*q1tL*P4*m1*m2 - 4*q2(vb)*q1tq2t*q1tL*P4
     & *k1k2 - 4*q2(vb)*q1tq2t*q1tL*P7*m1*m2 - 4*q2(vb)*q1tq2t*q1tL
     & *P7*k1k2 + 4*q2(vb)*q1tq2t*q2tL*P4*m1*m2 + 4*q2(vb)*q1tq2t*
     & q2tL*P4*k1k2 - 4*q2(vb)*q1tk1t*q2tL*P4*q1k2 - 4*q2(vb)*
     & q1tk1t*q2tL*P4*q2k2 + 8*q2(vb)*q1tk1t*q2tL*P7*m1*m2 - 4*
     & q2(vb)*q1tk1t*q2tL*P7*q1k2
      AS1(3,vb) = AS1(3,vb) + 8*q2(vb)*q1tk1t*q2tL*P7*k1k2 - 8*
     & q2(vb)*q1tk2t*q2tL*P4*m1*m2 + 4*q2(vb)*q1tk2t*q2tL*P4*q1k1
     &  + 4*q2(vb)*q1tk2t*q2tL*P4*q2k1 - 8*q2(vb)*q1tk2t*q2tL*P4*
     & k1k2 + 4*q2(vb)*q1tk2t*q2tL*P7*q1k1 + 4*q2(vb)*q1tL*q2tk1t*
     & P4*q1k2 + 4*q2(vb)*q1tL*q2tk1t*P4*q2k2 + 4*q2(vb)*q1tL*
     & q2tk1t*P7*q1k2 + 8*q2(vb)*q1tL*q2tk2t*P4*m1*m2 - 4*q2(vb)*
     & q1tL*q2tk2t*P4*q1k1 - 4*q2(vb)*q1tL*q2tk2t*P4*q2k1 + 8*
     & q2(vb)*q1tL*q2tk2t*P4*k1k2 - 4*q2(vb)*q1tL*q2tk2t*P7*q1k1 + 
     & 4*q2(vb)*q1tL*t2*P4*m1*m2 + 4*q2(vb)*q1tL*t2*P4*k1k2 - 4*
     & q2(vb)*q2tL*t1*P4*m1*m2 - 4*q2(vb)*q2tL*t1*P4*k1k2 + 4*q2(vb)*
     & q2tL*t1*P7*m1*m2 + 4*q2(vb)*q2tL*t1*P7*k1k2 - 8*q1t(vb)*q1L*
     & q1tq2t*P3*m1*m2 - 8*q1t(vb)*q1L*q1tq2t*P3*k1k2 + 8*q1t(vb)*
     & q1L*q2tk1t*P3*m1*m2 - 4*q1t(vb)*q1L*q2tk1t*P3*q2k2 + 8*
     & q1t(vb)*q1L*q2tk1t*P3*k1k2 + 4*q1t(vb)*q1L*q2tk1t*P4*q2k2 - 
     & 4*q1t(vb)*q1L*q2tk1t*P7*q2k2 + 4*q1t(vb)*q1L*q2tk2t*P3*q2k1
     &  + 8*q1t(vb)*q1L*q2tk2t*P4*m1*m2
      AS1(3,vb) = AS1(3,vb) - 4*q1t(vb)*q1L*q2tk2t*P4*q2k1 + 8*
     & q1t(vb)*q1L*q2tk2t*P4*k1k2 - 8*q1t(vb)*q1L*q2tk2t*P7*m1*m2
     &  + 4*q1t(vb)*q1L*q2tk2t*P7*q2k1 - 8*q1t(vb)*q1L*q2tk2t*P7*
     & k1k2 + 4*q1t(vb)*q1L*t2*P3*m1*m2 + 4*q1t(vb)*q1L*t2*P3*k1k2 + 
     & 4*q1t(vb)*q1L*t2*P4*m1*m2 + 4*q1t(vb)*q1L*t2*P4*k1k2 - 4*
     & q1t(vb)*q1L*t2*P7*m1*m2 - 4*q1t(vb)*q1L*t2*P7*k1k2 - 4*q1t(vb)
     & *q2L*q1tq2t*P3*m1*m2 - 4*q1t(vb)*q2L*q1tq2t*P3*k1k2 + 4*
     & q1t(vb)*q2L*q1tq2t*P4*m1*m2 + 4*q1t(vb)*q2L*q1tq2t*P4*k1k2
     &  - 4*q1t(vb)*q2L*q1tq2t*P7*m1*m2 - 4*q1t(vb)*q2L*q1tq2t*P7*
     & k1k2 + 4*q1t(vb)*q2L*q2tk1t*P3*q1k2 - 4*q1t(vb)*q2L*q2tk1t*
     & P4*q1k2 + 4*q1t(vb)*q2L*q2tk1t*P7*q1k2 - 4*q1t(vb)*q2L*
     & q2tk2t*P3*q1k1 + 4*q1t(vb)*q2L*q2tk2t*P4*q1k1 - 4*q1t(vb)*
     & q2L*q2tk2t*P7*q1k1 - 4*q1t(vb)*q2tL*t1*P3*m1*m2 - 4*q1t(vb)*
     & q2tL*t1*P3*k1k2 + 4*q1t(vb)*q2tL*t2*P4*m1*m2 + 4*q1t(vb)*q2tL
     & *t2*P4*k1k2 - 4*q1t(vb)*q2tL*P3*m1*m2*q1q2 + 8*q1t(vb)*q2tL*P3
     & *m1*m2*q1k1
      AS1(3,vb) = AS1(3,vb) - 4*q1t(vb)*q2tL*P3*q1q2*k1k2 - 4*
     & q1t(vb)*q2tL*P3*q1k1*q2k2 + 8*q1t(vb)*q2tL*P3*q1k1*k1k2 + 4*
     & q1t(vb)*q2tL*P3*q1k2*q2k1 + 4*q1t(vb)*q2tL*P4*m1*m2*q1q2 - 8*
     & q1t(vb)*q2tL*P4*m1*m2*q2k2 + 4*q1t(vb)*q2tL*P4*q1q2*k1k2 + 4*
     & q1t(vb)*q2tL*P4*q1k1*q2k2 - 4*q1t(vb)*q2tL*P4*q1k2*q2k1 - 8*
     & q1t(vb)*q2tL*P4*q2k2*k1k2 - 4*q1t(vb)*q2tL*P7*m1*m2*q1q2 - 4*
     & q1t(vb)*q2tL*P7*q1q2*k1k2 - 4*q1t(vb)*q2tL*P7*q1k1*q2k2 + 4*
     & q1t(vb)*q2tL*P7*q1k2*q2k1 + 4*q2t(vb)*q1L*q1tq2t*P3*m1*m2 + 4
     & *q2t(vb)*q1L*q1tq2t*P3*k1k2 - 4*q2t(vb)*q1L*q1tq2t*P4*m1*m2
     &  - 4*q2t(vb)*q1L*q1tq2t*P4*k1k2 - 4*q2t(vb)*q1L*q1tq2t*P7*m1
     & *m2 - 4*q2t(vb)*q1L*q1tq2t*P7*k1k2 + 4*q2t(vb)*q1L*q1tk1t*P3
     & *q2k2 - 4*q2t(vb)*q1L*q1tk1t*P4*q2k2 - 4*q2t(vb)*q1L*q1tk1t*
     & P7*q2k2 - 4*q2t(vb)*q1L*q1tk2t*P3*q2k1 + 4*q2t(vb)*q1L*
     & q1tk2t*P4*q2k1 + 4*q2t(vb)*q1L*q1tk2t*P7*q2k1 - 8*q2t(vb)*
     & q2L*q1tq2t*P4*m1*m2 - 8*q2t(vb)*q2L*q1tq2t*P4*k1k2 + 8*
     & q2t(vb)*q2L*q1tk1t*P3*m1*m2
      AS1(3,vb) = AS1(3,vb) - 4*q2t(vb)*q2L*q1tk1t*P3*q1k2 + 8*
     & q2t(vb)*q2L*q1tk1t*P3*k1k2 + 4*q2t(vb)*q2L*q1tk1t*P4*q1k2 - 
     & 8*q2t(vb)*q2L*q1tk1t*P7*m1*m2 + 4*q2t(vb)*q2L*q1tk1t*P7*q1k2
     &  - 8*q2t(vb)*q2L*q1tk1t*P7*k1k2 + 4*q2t(vb)*q2L*q1tk2t*P3*
     & q1k1 + 8*q2t(vb)*q2L*q1tk2t*P4*m1*m2 - 4*q2t(vb)*q2L*q1tk2t*
     & P4*q1k1 + 8*q2t(vb)*q2L*q1tk2t*P4*k1k2 - 4*q2t(vb)*q2L*
     & q1tk2t*P7*q1k1 + 4*q2t(vb)*q2L*t1*P3*m1*m2 + 4*q2t(vb)*q2L*t1
     & *P3*k1k2 + 4*q2t(vb)*q2L*t1*P4*m1*m2 + 4*q2t(vb)*q2L*t1*P4*
     & k1k2 - 4*q2t(vb)*q2L*t1*P7*m1*m2 - 4*q2t(vb)*q2L*t1*P7*k1k2 + 
     & 4*q2t(vb)*q1tL*t1*P3*m1*m2 + 4*q2t(vb)*q1tL*t1*P3*k1k2 - 4*
     & q2t(vb)*q1tL*t2*P4*m1*m2 - 4*q2t(vb)*q1tL*t2*P4*k1k2 + 4*
     & q2t(vb)*q1tL*P3*m1*m2*q1q2 - 8*q2t(vb)*q1tL*P3*m1*m2*q1k1 + 4*
     & q2t(vb)*q1tL*P3*q1q2*k1k2 + 4*q2t(vb)*q1tL*P3*q1k1*q2k2 - 8*
     & q2t(vb)*q1tL*P3*q1k1*k1k2 - 4*q2t(vb)*q1tL*P3*q1k2*q2k1 - 4*
     & q2t(vb)*q1tL*P4*m1*m2*q1q2 + 8*q2t(vb)*q1tL*P4*m1*m2*q2k2 - 4*
     & q2t(vb)*q1tL*P4*q1q2*k1k2
      AS1(3,vb) = AS1(3,vb) - 4*q2t(vb)*q1tL*P4*q1k1*q2k2 + 4*
     & q2t(vb)*q1tL*P4*q1k2*q2k1 + 8*q2t(vb)*q1tL*P4*q2k2*k1k2 - 4*
     & q2t(vb)*q1tL*P7*m1*m2*q1q2 - 4*q2t(vb)*q1tL*P7*q1q2*k1k2 - 4*
     & q2t(vb)*q1tL*P7*q1k1*q2k2 + 4*q2t(vb)*q1tL*P7*q1k2*q2k1 - 4*
     & L(vb)*t1*t2*P3*m1*m2 - 4*L(vb)*t1*t2*P3*k1k2 - 4*L(vb)*t1*t2*P4*
     & m1*m2 - 4*L(vb)*t1*t2*P4*k1k2 + 4*L(vb)*t1*t2*P7*m1*m2 + 4*L(vb)
     & *t1*t2*P7*k1k2 + 16*L(vb)*q1tq2t*q1tk1t*P3*m1*m2 - 8*L(vb)*
     & q1tq2t*q1tk1t*P3*q1k2 - 4*L(vb)*q1tq2t*q1tk1t*P3*q2k2 + 16*
     & L(vb)*q1tq2t*q1tk1t*P3*k1k2 - 4*L(vb)*q1tq2t*q1tk1t*P4*q2k2
     &  - 4*L(vb)*q1tq2t*q1tk1t*P7*q2k2 + 8*L(vb)*q1tq2t*q1tk2t*P3*
     & q1k1 + 4*L(vb)*q1tq2t*q1tk2t*P3*q2k1 + 4*L(vb)*q1tq2t*q1tk2t
     & *P4*q2k1 + 4*L(vb)*q1tq2t*q1tk2t*P7*q2k1 + 4*L(vb)*q1tq2t*
     & q2tk1t*P3*q1k2 + 4*L(vb)*q1tq2t*q2tk1t*P4*q1k2 + 8*L(vb)*
     & q1tq2t*q2tk1t*P4*q2k2 + 4*L(vb)*q1tq2t*q2tk1t*P7*q1k2 - 4*
     & L(vb)*q1tq2t*q2tk2t*P3*q1k1 + 16*L(vb)*q1tq2t*q2tk2t*P4*m1*
     & m2
      AS1(3,vb) = AS1(3,vb) - 4*L(vb)*q1tq2t*q2tk2t*P4*q1k1 - 8*
     & L(vb)*q1tq2t*q2tk2t*P4*q2k1 + 16*L(vb)*q1tq2t*q2tk2t*P4*k1k2
     &  - 4*L(vb)*q1tq2t*q2tk2t*P7*q1k1 + 4*L(vb)*q1tq2t*t1*P3*m1*m2
     &  + 4*L(vb)*q1tq2t*t1*P3*k1k2 + 4*L(vb)*q1tq2t*t2*P4*m1*m2 + 4*
     & L(vb)*q1tq2t*t2*P4*k1k2 - 4*L(vb)*q1tq2t*P3*m1*m2*q1q2 + 8*
     & L(vb)*q1tq2t*P3*m1*m2*q1k1 - 4*L(vb)*q1tq2t*P3*q1q2*k1k2 - 4*
     & L(vb)*q1tq2t*P3*q1k1*q2k2 + 8*L(vb)*q1tq2t*P3*q1k1*k1k2 + 4*
     & L(vb)*q1tq2t*P3*q1k2*q2k1 - 4*L(vb)*q1tq2t*P4*m1*m2*q1q2 + 8*
     & L(vb)*q1tq2t*P4*m1*m2*q2k2 - 4*L(vb)*q1tq2t*P4*q1q2*k1k2 - 4*
     & L(vb)*q1tq2t*P4*q1k1*q2k2 + 4*L(vb)*q1tq2t*P4*q1k2*q2k1 + 8*
     & L(vb)*q1tq2t*P4*q2k2*k1k2 - 4*L(vb)*q1tq2t*P7*m1*m2*q1q2 - 4*
     & L(vb)*q1tq2t*P7*q1q2*k1k2 - 4*L(vb)*q1tq2t*P7*q1k1*q2k2 + 4*
     & L(vb)*q1tq2t*P7*q1k2*q2k1 - 4*L(vb)*q1tq2t**2*P3*m1*m2 - 4*
     & L(vb)*q1tq2t**2*P3*k1k2 - 4*L(vb)*q1tq2t**2*P4*m1*m2 - 4*L(vb)
     & *q1tq2t**2*P4*k1k2 - 4*L(vb)*q1tq2t**2*P7*m1*m2 - 4*L(vb)*
     & q1tq2t**2*P7*k1k2
      AS1(3,vb) = AS1(3,vb) - 16*L(vb)*q1tk1t*q2tk1t*P3*m1*m2 + 8*
     & L(vb)*q1tk1t*q2tk1t*P3*q1k2 + 8*L(vb)*q1tk1t*q2tk1t*P3*q2k2
     &  - 16*L(vb)*q1tk1t*q2tk1t*P3*k1k2 + 8*L(vb)*q1tk1t*q2tk1t*P7
     & *q2k2 + 4*L(vb)*q1tk1t*q2tk2t*t1*P3 + 4*L(vb)*q1tk1t*q2tk2t*
     & t2*P4 + 4*L(vb)*q1tk1t*q2tk2t*P3*q1q2 - 8*L(vb)*q1tk1t*
     & q2tk2t*P3*q1k1 - 8*L(vb)*q1tk1t*q2tk2t*P3*q2k1 + 4*L(vb)*
     & q1tk1t*q2tk2t*P4*q1q2 - 8*L(vb)*q1tk1t*q2tk2t*P4*q1k2 - 8*
     & L(vb)*q1tk1t*q2tk2t*P4*q2k2 + 16*L(vb)*q1tk1t*q2tk2t*P7*m1*
     & m2 + 4*L(vb)*q1tk1t*q2tk2t*P7*q1q2 - 8*L(vb)*q1tk1t*q2tk2t*
     & P7*q1k2 - 8*L(vb)*q1tk1t*q2tk2t*P7*q2k1 + 16*L(vb)*q1tk1t*
     & q2tk2t*P7*k1k2 - 8*L(vb)*q1tk1t*t2*P3*m1*m2 + 4*L(vb)*q1tk1t*
     & t2*P3*q1k2 - 8*L(vb)*q1tk1t*t2*P3*k1k2 - 4*L(vb)*q1tk1t*t2*P4*
     & q1k2 + 8*L(vb)*q1tk1t*t2*P7*m1*m2 - 4*L(vb)*q1tk1t*t2*P7*q1k2
     &  + 8*L(vb)*q1tk1t*t2*P7*k1k2 - 4*L(vb)*q1tk2t*q2tk1t*t1*P3 - 
     & 4*L(vb)*q1tk2t*q2tk1t*t2*P4 - 4*L(vb)*q1tk2t*q2tk1t*P3*q1q2
     &  - 4*L(vb)*q1tk2t*q2tk1t*P4*q1q2
      AS1(3,vb) = AS1(3,vb) - 4*L(vb)*q1tk2t*q2tk1t*P7*q1q2 - 16*
     & L(vb)*q1tk2t*q2tk2t*P4*m1*m2 + 8*L(vb)*q1tk2t*q2tk2t*P4*q1k1
     &  + 8*L(vb)*q1tk2t*q2tk2t*P4*q2k1 - 16*L(vb)*q1tk2t*q2tk2t*P4
     & *k1k2 + 8*L(vb)*q1tk2t*q2tk2t*P7*q1k1 - 4*L(vb)*q1tk2t*t2*P3*
     & q1k1 - 8*L(vb)*q1tk2t*t2*P4*m1*m2 + 4*L(vb)*q1tk2t*t2*P4*q1k1
     &  - 8*L(vb)*q1tk2t*t2*P4*k1k2 + 4*L(vb)*q1tk2t*t2*P7*q1k1 - 8*
     & L(vb)*q2tk1t*t1*P3*m1*m2 + 4*L(vb)*q2tk1t*t1*P3*q2k2 - 8*L(vb)
     & *q2tk1t*t1*P3*k1k2 - 4*L(vb)*q2tk1t*t1*P4*q2k2 + 4*L(vb)*
     & q2tk1t*t1*P7*q2k2 - 4*L(vb)*q2tk2t*t1*P3*q2k1 - 8*L(vb)*
     & q2tk2t*t1*P4*m1*m2 + 4*L(vb)*q2tk2t*t1*P4*q2k1 - 8*L(vb)*
     & q2tk2t*t1*P4*k1k2 + 8*L(vb)*q2tk2t*t1*P7*m1*m2 - 4*L(vb)*
     & q2tk2t*t1*P7*q2k1 + 8*L(vb)*q2tk2t*t1*P7*k1k2

      AS1(4,vb) =
     &  - 4*k1(vb)*q1L*q2U*q1tq2t*P3*m2 - 4*k1(vb)*q1L*q2U*q1tq2t
     & *P4*m2 - 4*k1(vb)*q1L*q2U*q1tq2t*P7*m2 - 8*k1(vb)*q1L*
     & q1tq2t*q1tU*P3*m2 + 4*k1(vb)*q1L*q1tq2t*q2tU*P3*m2 + 4*
     & k1(vb)*q1L*q1tq2t*q2tU*P4*m2 + 4*k1(vb)*q1L*q1tq2t*q2tU*P7
     & *m2 + 4*k1(vb)*q1L*q1tk1t*q2tU*P3*m2 - 4*k1(vb)*q1L*q1tk2t*
     & q2tU*P3*m1 + 4*k1(vb)*q1L*q1tU*q2tk1t*P3*m2 + 4*k1(vb)*q1L*
     & q1tU*q2tk2t*P3*m1 - 8*k1(vb)*q1L*q1tU*q2tk2t*P4*m2 - 8*
     & k1(vb)*q1L*q1tU*q2tk2t*P7*m2 + 4*k1(vb)*q1L*q1tU*t2*P3*m2
     &  - 4*k1(vb)*q1L*q1tU*t2*P4*m2 - 4*k1(vb)*q1L*q1tU*t2*P7*m2
     &  - 4*k1(vb)*q2U*q2L*q1tq2t*P3*m2 - 4*k1(vb)*q2U*q2L*q1tq2t
     & *P4*m2 - 4*k1(vb)*q2U*q2L*q1tq2t*P7*m2 - 4*k1(vb)*q2U*
     & q1tq2t*q1tL*P3*m2 + 4*k1(vb)*q2U*q1tq2t*q1tL*P4*m2 + 4*
     & k1(vb)*q2U*q1tq2t*q1tL*P7*m2 + 4*k1(vb)*q2U*q1tq2t*q2tL*P3
     & *m2 - 4*k1(vb)*q2U*q1tq2t*q2tL*P4*m2 + 4*k1(vb)*q2U*q1tq2t*
     & q2tL*P7*m2
      AS1(4,vb) = AS1(4,vb) - 4*k1(vb)*q2U*q1tk1t*q2tL*P3*m2 - 4*
     & k1(vb)*q2U*q1tk1t*q2tL*P4*m2 - 8*k1(vb)*q2U*q1tk1t*q2tL*P7
     & *m2 - 4*k1(vb)*q2U*q1tk2t*q2tL*P3*m1 + 4*k1(vb)*q2U*q1tk2t*
     & q2tL*P4*m1 + 8*k1(vb)*q2U*q1tk2t*q2tL*P4*m2 + 4*k1(vb)*q2U*
     & q1tL*q2tk1t*P3*m2 + 4*k1(vb)*q2U*q1tL*q2tk1t*P4*m2 + 4*
     & k1(vb)*q2U*q1tL*q2tk2t*P3*m1 - 4*k1(vb)*q2U*q1tL*q2tk2t*P4
     & *m1 - 8*k1(vb)*q2U*q1tL*q2tk2t*P4*m2 - 8*k1(vb)*q2U*q1tL*
     & q2tk2t*P7*m2 + 4*k1(vb)*q2U*q1tL*t2*P3*m2 - 4*k1(vb)*q2U*
     & q1tL*t2*P4*m2 - 4*k1(vb)*q2U*q1tL*t2*P7*m2 - 4*k1(vb)*q2U*
     & q2tL*t1*P3*m2 + 4*k1(vb)*q2U*q2tL*t1*P4*m2 - 4*k1(vb)*q2U*
     & q2tL*t1*P7*m2 - 4*k1(vb)*q2L*q1tq2t*q1tU*P3*m2 - 4*k1(vb)*
     & q2L*q1tq2t*q1tU*P4*m2 - 4*k1(vb)*q2L*q1tq2t*q1tU*P7*m2 + 8
     & *k1(vb)*q2L*q1tq2t*q2tU*P4*m2 + 8*k1(vb)*q2L*q1tk1t*q2tU*
     & P3*m2 + 4*k1(vb)*q2L*q1tk1t*q2tU*P4*m2 + 8*k1(vb)*q2L*
     & q1tk1t*q2tU*P7*m2 - 4*k1(vb)*q2L*q1tk2t*q2tU*P4*m1 - 8*
     & k1(vb)*q2L*q1tk2t*q2tU*P4*m2
      AS1(4,vb) = AS1(4,vb) - 4*k1(vb)*q2L*q1tU*q2tk1t*P4*m2 + 4*
     & k1(vb)*q2L*q1tU*q2tk2t*P4*m1 + 4*k1(vb)*q2L*q2tU*t1*P3*m2
     &  - 4*k1(vb)*q2L*q2tU*t1*P4*m2 + 4*k1(vb)*q2L*q2tU*t1*P7*m2
     &  - 4*k1(vb)*q1tU*q2tL*t1*P3*m2 - 4*k1(vb)*q1tU*q2tL*t2*P4*m2
     &  + 4*k1(vb)*q1tU*q2tL*P3*m1*q1k2 - 4*k1(vb)*q1tU*q2tL*P3*m2*
     & q1q2 + 4*k1(vb)*q1tU*q2tL*P3*m2*q1k1 + 4*k1(vb)*q1tU*q2tL*P4
     & *m1*q2k2 - 4*k1(vb)*q1tU*q2tL*P4*m2*q1q2 - 4*k1(vb)*q1tU*
     & q2tL*P4*m2*q2k1 + 8*k1(vb)*q1tU*q2tL*P4*m2*q2k2 - 4*k1(vb)*
     & q1tU*q2tL*P7*m2*q1q2 + 4*k1(vb)*q1tL*q2tU*t1*P3*m2 + 4*
     & k1(vb)*q1tL*q2tU*t2*P4*m2 - 4*k1(vb)*q1tL*q2tU*P3*m1*q1k2 + 
     & 4*k1(vb)*q1tL*q2tU*P3*m2*q1q2 - 4*k1(vb)*q1tL*q2tU*P3*m2*
     & q1k1 - 4*k1(vb)*q1tL*q2tU*P4*m1*q2k2 + 4*k1(vb)*q1tL*q2tU*P4
     & *m2*q1q2 + 4*k1(vb)*q1tL*q2tU*P4*m2*q2k1 - 8*k1(vb)*q1tL*
     & q2tU*P4*m2*q2k2 + 4*k1(vb)*q1tL*q2tU*P7*m2*q1q2 + 4*k2(vb)*
     & q1L*q2U*q1tq2t*P3*m1 + 4*k2(vb)*q1L*q2U*q1tq2t*P4*m1 + 4*
     & k2(vb)*q1L*q2U*q1tq2t*P7*m1
      AS1(4,vb) = AS1(4,vb) + 8*k2(vb)*q1L*q1tq2t*q1tU*P3*m1 - 4*
     & k2(vb)*q1L*q1tq2t*q2tU*P3*m1 - 4*k2(vb)*q1L*q1tq2t*q2tU*P4
     & *m1 - 4*k2(vb)*q1L*q1tq2t*q2tU*P7*m1 + 4*k2(vb)*q1L*q1tk1t*
     & q2tU*P3*m2 - 4*k2(vb)*q1L*q1tk2t*q2tU*P3*m1 - 8*k2(vb)*q1L*
     & q1tU*q2tk1t*P3*m1 - 4*k2(vb)*q1L*q1tU*q2tk1t*P3*m2 + 4*
     & k2(vb)*q1L*q1tU*q2tk2t*P3*m1 + 8*k2(vb)*q1L*q1tU*q2tk2t*P4
     & *m1 + 8*k2(vb)*q1L*q1tU*q2tk2t*P7*m1 - 4*k2(vb)*q1L*q1tU*t2
     & *P3*m1 + 4*k2(vb)*q1L*q1tU*t2*P4*m1 + 4*k2(vb)*q1L*q1tU*t2*
     & P7*m1 + 4*k2(vb)*q2U*q2L*q1tq2t*P3*m1 + 4*k2(vb)*q2U*q2L*
     & q1tq2t*P4*m1 + 4*k2(vb)*q2U*q2L*q1tq2t*P7*m1 + 4*k2(vb)*q2U
     & *q1tq2t*q1tL*P3*m1 - 4*k2(vb)*q2U*q1tq2t*q1tL*P4*m1 - 4*
     & k2(vb)*q2U*q1tq2t*q1tL*P7*m1 - 4*k2(vb)*q2U*q1tq2t*q2tL*P3
     & *m1 + 4*k2(vb)*q2U*q1tq2t*q2tL*P4*m1 - 4*k2(vb)*q2U*q1tq2t*
     & q2tL*P7*m1 + 8*k2(vb)*q2U*q1tk1t*q2tL*P3*m1 + 4*k2(vb)*q2U*
     & q1tk1t*q2tL*P3*m2 - 4*k2(vb)*q2U*q1tk1t*q2tL*P4*m2 + 8*
     & k2(vb)*q2U*q1tk1t*q2tL*P7*m1
      AS1(4,vb) = AS1(4,vb) - 4*k2(vb)*q2U*q1tk2t*q2tL*P3*m1 - 4*
     & k2(vb)*q2U*q1tk2t*q2tL*P4*m1 - 8*k2(vb)*q2U*q1tL*q2tk1t*P3
     & *m1 - 4*k2(vb)*q2U*q1tL*q2tk1t*P3*m2 + 4*k2(vb)*q2U*q1tL*
     & q2tk1t*P4*m2 + 4*k2(vb)*q2U*q1tL*q2tk2t*P3*m1 + 4*k2(vb)*
     & q2U*q1tL*q2tk2t*P4*m1 + 8*k2(vb)*q2U*q1tL*q2tk2t*P7*m1 - 4
     & *k2(vb)*q2U*q1tL*t2*P3*m1 + 4*k2(vb)*q2U*q1tL*t2*P4*m1 + 4*
     & k2(vb)*q2U*q1tL*t2*P7*m1 + 4*k2(vb)*q2U*q2tL*t1*P3*m1 - 4*
     & k2(vb)*q2U*q2tL*t1*P4*m1 + 4*k2(vb)*q2U*q2tL*t1*P7*m1 + 4*
     & k2(vb)*q2L*q1tq2t*q1tU*P3*m1 + 4*k2(vb)*q2L*q1tq2t*q1tU*P4
     & *m1 + 4*k2(vb)*q2L*q1tq2t*q1tU*P7*m1 - 8*k2(vb)*q2L*q1tq2t*
     & q2tU*P4*m1 - 8*k2(vb)*q2L*q1tk1t*q2tU*P3*m1 + 4*k2(vb)*q2L*
     & q1tk1t*q2tU*P4*m2 - 8*k2(vb)*q2L*q1tk1t*q2tU*P7*m1 + 4*
     & k2(vb)*q2L*q1tk2t*q2tU*P4*m1 - 4*k2(vb)*q2L*q1tU*q2tk1t*P4
     & *m2 + 4*k2(vb)*q2L*q1tU*q2tk2t*P4*m1 - 4*k2(vb)*q2L*q2tU*t1
     & *P3*m1 + 4*k2(vb)*q2L*q2tU*t1*P4*m1 - 4*k2(vb)*q2L*q2tU*t1*
     & P7*m1
      AS1(4,vb) = AS1(4,vb) + 4*k2(vb)*q1tU*q2tL*t1*P3*m1 + 4*
     & k2(vb)*q1tU*q2tL*t2*P4*m1 + 4*k2(vb)*q1tU*q2tL*P3*m1*q1q2 - 
     & 8*k2(vb)*q1tU*q2tL*P3*m1*q1k1 + 4*k2(vb)*q1tU*q2tL*P3*m1*
     & q1k2 - 4*k2(vb)*q1tU*q2tL*P3*m2*q1k1 + 4*k2(vb)*q1tU*q2tL*P4
     & *m1*q1q2 - 4*k2(vb)*q1tU*q2tL*P4*m1*q2k2 - 4*k2(vb)*q1tU*
     & q2tL*P4*m2*q2k1 + 4*k2(vb)*q1tU*q2tL*P7*m1*q1q2 - 4*k2(vb)*
     & q1tL*q2tU*t1*P3*m1 - 4*k2(vb)*q1tL*q2tU*t2*P4*m1 - 4*k2(vb)*
     & q1tL*q2tU*P3*m1*q1q2 + 8*k2(vb)*q1tL*q2tU*P3*m1*q1k1 - 4*
     & k2(vb)*q1tL*q2tU*P3*m1*q1k2 + 4*k2(vb)*q1tL*q2tU*P3*m2*q1k1
     &  - 4*k2(vb)*q1tL*q2tU*P4*m1*q1q2 + 4*k2(vb)*q1tL*q2tU*P4*m1*
     & q2k2 + 4*k2(vb)*q1tL*q2tU*P4*m2*q2k1 - 4*k2(vb)*q1tL*q2tU*P7
     & *m1*q1q2 + 4*kw(vb)*q1L*q1tk1t*q2tU*P3*m2 - 4*kw(vb)*q1L*
     & q1tk2t*q2tU*P3*m1 - 4*kw(vb)*q1L*q1tU*q2tk1t*P3*m2 + 4*
     & kw(vb)*q1L*q1tU*q2tk2t*P3*m1 + 4*kw(vb)*q2U*q1tk1t*q2tL*P3
     & *m2 - 4*kw(vb)*q2U*q1tk1t*q2tL*P4*m2 - 4*kw(vb)*q2U*q1tk2t*
     & q2tL*P3*m1
      AS1(4,vb) = AS1(4,vb) + 4*kw(vb)*q2U*q1tk2t*q2tL*P4*m1 - 4*
     & kw(vb)*q2U*q1tL*q2tk1t*P3*m2 + 4*kw(vb)*q2U*q1tL*q2tk1t*P4
     & *m2 + 4*kw(vb)*q2U*q1tL*q2tk2t*P3*m1 - 4*kw(vb)*q2U*q1tL*
     & q2tk2t*P4*m1 + 4*kw(vb)*q2L*q1tk1t*q2tU*P4*m2 - 4*kw(vb)*
     & q2L*q1tk2t*q2tU*P4*m1 - 4*kw(vb)*q2L*q1tU*q2tk1t*P4*m2 + 4
     & *kw(vb)*q2L*q1tU*q2tk2t*P4*m1 + 4*kw(vb)*q1tU*q2tL*P3*m1*
     & q1k2 - 4*kw(vb)*q1tU*q2tL*P3*m2*q1k1 + 4*kw(vb)*q1tU*q2tL*P4
     & *m1*q2k2 - 4*kw(vb)*q1tU*q2tL*P4*m2*q2k1 - 4*kw(vb)*q1tL*
     & q2tU*P3*m1*q1k2 + 4*kw(vb)*q1tL*q2tU*P3*m2*q1k1 - 4*kw(vb)*
     & q1tL*q2tU*P4*m1*q2k2 + 4*kw(vb)*q1tL*q2tU*P4*m2*q2k1 - 4*
     & q1(vb)*q1L*q1tk1t*q2tU*P3*m2 + 4*q1(vb)*q1L*q1tk2t*q2tU*P3
     & *m1 + 4*q1(vb)*q1L*q1tU*q2tk1t*P3*m2 - 4*q1(vb)*q1L*q1tU*
     & q2tk2t*P3*m1 + 4*q1(vb)*q2U*q1tk1t*q2tL*P7*m2 - 4*q1(vb)*
     & q2U*q1tk2t*q2tL*P7*m1 - 4*q1(vb)*q2U*q1tL*q2tk1t*P7*m2 + 4
     & *q1(vb)*q2U*q1tL*q2tk2t*P7*m1 - 4*q1(vb)*q2L*q1tk1t*q2tU*
     & P3*m2
      AS1(4,vb) = AS1(4,vb) - 4*q1(vb)*q2L*q1tk1t*q2tU*P7*m2 + 4*
     & q1(vb)*q2L*q1tk2t*q2tU*P3*m1 + 4*q1(vb)*q2L*q1tk2t*q2tU*P7
     & *m1 + 4*q1(vb)*q2L*q1tU*q2tk1t*P3*m2 + 4*q1(vb)*q2L*q1tU*
     & q2tk1t*P7*m2 - 4*q1(vb)*q2L*q1tU*q2tk2t*P3*m1 - 4*q1(vb)*
     & q2L*q1tU*q2tk2t*P7*m1 - 4*q1(vb)*q1tU*q2tL*P3*m1*q1k2 - 4*
     & q1(vb)*q1tU*q2tL*P3*m1*q2k2 + 4*q1(vb)*q1tU*q2tL*P3*m2*q1k1
     &  + 4*q1(vb)*q1tU*q2tL*P3*m2*q2k1 - 4*q1(vb)*q1tU*q2tL*P7*m1*
     & q2k2 + 4*q1(vb)*q1tU*q2tL*P7*m2*q2k1 + 4*q1(vb)*q1tL*q2tU*P3
     & *m1*q1k2 + 4*q1(vb)*q1tL*q2tU*P3*m1*q2k2 - 4*q1(vb)*q1tL*
     & q2tU*P3*m2*q1k1 - 4*q1(vb)*q1tL*q2tU*P3*m2*q2k1 + 4*q1(vb)*
     & q1tL*q2tU*P7*m1*q2k2 - 4*q1(vb)*q1tL*q2tU*P7*m2*q2k1 - 4*
     & q2(vb)*q1L*q1tk1t*q2tU*P4*m2 - 4*q2(vb)*q1L*q1tk1t*q2tU*P7
     & *m2 + 4*q2(vb)*q1L*q1tk2t*q2tU*P4*m1 + 4*q2(vb)*q1L*q1tk2t*
     & q2tU*P7*m1 + 4*q2(vb)*q1L*q1tU*q2tk1t*P4*m2 + 4*q2(vb)*q1L*
     & q1tU*q2tk1t*P7*m2 - 4*q2(vb)*q1L*q1tU*q2tk2t*P4*m1 - 4*
     & q2(vb)*q1L*q1tU*q2tk2t*P7*m1
      AS1(4,vb) = AS1(4,vb) - 4*q2(vb)*q2U*q1tk1t*q2tL*P7*m2 + 4*
     & q2(vb)*q2U*q1tk2t*q2tL*P7*m1 + 4*q2(vb)*q2U*q1tL*q2tk1t*P7
     & *m2 - 4*q2(vb)*q2U*q1tL*q2tk2t*P7*m1 - 4*q2(vb)*q2L*q1tk1t*
     & q2tU*P4*m2 + 4*q2(vb)*q2L*q1tk2t*q2tU*P4*m1 + 4*q2(vb)*q2L*
     & q1tU*q2tk1t*P4*m2 - 4*q2(vb)*q2L*q1tU*q2tk2t*P4*m1 - 4*
     & q2(vb)*q1tU*q2tL*P4*m1*q1k2 - 4*q2(vb)*q1tU*q2tL*P4*m1*q2k2
     &  + 4*q2(vb)*q1tU*q2tL*P4*m2*q1k1 + 4*q2(vb)*q1tU*q2tL*P4*m2*
     & q2k1 - 4*q2(vb)*q1tU*q2tL*P7*m1*q1k2 + 4*q2(vb)*q1tU*q2tL*P7
     & *m2*q1k1 + 4*q2(vb)*q1tL*q2tU*P4*m1*q1k2 + 4*q2(vb)*q1tL*
     & q2tU*P4*m1*q2k2 - 4*q2(vb)*q1tL*q2tU*P4*m2*q1k1 - 4*q2(vb)*
     & q1tL*q2tU*P4*m2*q2k1 + 4*q2(vb)*q1tL*q2tU*P7*m1*q1k2 - 4*
     & q2(vb)*q1tL*q2tU*P7*m2*q1k1 + 4*q1t(vb)*q1L*q2U*q2tk1t*P3*
     & m2 - 4*q1t(vb)*q1L*q2U*q2tk1t*P4*m2 + 4*q1t(vb)*q1L*q2U*
     & q2tk1t*P7*m2 - 4*q1t(vb)*q1L*q2U*q2tk2t*P3*m1 + 4*q1t(vb)*
     & q1L*q2U*q2tk2t*P4*m1 - 4*q1t(vb)*q1L*q2U*q2tk2t*P7*m1 - 4*
     & q1t(vb)*q1L*q2tU*P3*m1*q2k2
      AS1(4,vb) = AS1(4,vb) + 4*q1t(vb)*q1L*q2tU*P3*m2*q2k1 + 4*
     & q1t(vb)*q1L*q2tU*P4*m1*q2k2 - 4*q1t(vb)*q1L*q2tU*P4*m2*q2k1
     &  - 4*q1t(vb)*q1L*q2tU*P7*m1*q2k2 + 4*q1t(vb)*q1L*q2tU*P7*m2*
     & q2k1 + 4*q1t(vb)*q2U*q2L*q2tk1t*P3*m2 - 4*q1t(vb)*q2U*q2L*
     & q2tk1t*P4*m2 + 4*q1t(vb)*q2U*q2L*q2tk1t*P7*m2 - 4*q1t(vb)*
     & q2U*q2L*q2tk2t*P3*m1 + 4*q1t(vb)*q2U*q2L*q2tk2t*P4*m1 - 4*
     & q1t(vb)*q2U*q2L*q2tk2t*P7*m1 - 4*q1t(vb)*q2U*q2tL*P3*m1*
     & q1k2 - 4*q1t(vb)*q2U*q2tL*P3*m1*q2k2 + 4*q1t(vb)*q2U*q2tL*P3
     & *m2*q1k1 + 4*q1t(vb)*q2U*q2tL*P3*m2*q2k1 + 4*q1t(vb)*q2U*
     & q2tL*P4*m1*q1k2 + 4*q1t(vb)*q2U*q2tL*P4*m1*q2k2 - 4*q1t(vb)*
     & q2U*q2tL*P4*m2*q1k1 - 4*q1t(vb)*q2U*q2tL*P4*m2*q2k1 - 4*
     & q1t(vb)*q2U*q2tL*P7*m1*q1k2 - 4*q1t(vb)*q2U*q2tL*P7*m1*q2k2
     &  + 4*q1t(vb)*q2U*q2tL*P7*m2*q1k1 + 4*q1t(vb)*q2U*q2tL*P7*m2*
     & q2k1 + 4*q1t(vb)*q2L*q2tU*P3*m1*q1k2 - 4*q1t(vb)*q2L*q2tU*P3
     & *m2*q1k1 - 4*q1t(vb)*q2L*q2tU*P4*m1*q1k2 + 4*q1t(vb)*q2L*
     & q2tU*P4*m2*q1k1
      AS1(4,vb) = AS1(4,vb) + 4*q1t(vb)*q2L*q2tU*P7*m1*q1k2 - 4*
     & q1t(vb)*q2L*q2tU*P7*m2*q1k1 - 4*q2t(vb)*q1L*q2U*q1tk1t*P3*
     & m2 + 4*q2t(vb)*q1L*q2U*q1tk1t*P4*m2 + 4*q2t(vb)*q1L*q2U*
     & q1tk1t*P7*m2 + 4*q2t(vb)*q1L*q2U*q1tk2t*P3*m1 - 4*q2t(vb)*
     & q1L*q2U*q1tk2t*P4*m1 - 4*q2t(vb)*q1L*q2U*q1tk2t*P7*m1 + 4*
     & q2t(vb)*q1L*q1tU*P3*m1*q2k2 - 4*q2t(vb)*q1L*q1tU*P3*m2*q2k1
     &  - 4*q2t(vb)*q1L*q1tU*P4*m1*q2k2 + 4*q2t(vb)*q1L*q1tU*P4*m2*
     & q2k1 - 4*q2t(vb)*q1L*q1tU*P7*m1*q2k2 + 4*q2t(vb)*q1L*q1tU*P7
     & *m2*q2k1 - 4*q2t(vb)*q2U*q2L*q1tk1t*P3*m2 + 4*q2t(vb)*q2U*
     & q2L*q1tk1t*P4*m2 + 4*q2t(vb)*q2U*q2L*q1tk1t*P7*m2 + 4*
     & q2t(vb)*q2U*q2L*q1tk2t*P3*m1 - 4*q2t(vb)*q2U*q2L*q1tk2t*P4
     & *m1 - 4*q2t(vb)*q2U*q2L*q1tk2t*P7*m1 + 4*q2t(vb)*q2U*q1tL*
     & P3*m1*q1k2 + 4*q2t(vb)*q2U*q1tL*P3*m1*q2k2 - 4*q2t(vb)*q2U*
     & q1tL*P3*m2*q1k1 - 4*q2t(vb)*q2U*q1tL*P3*m2*q2k1 - 4*q2t(vb)*
     & q2U*q1tL*P4*m1*q1k2 - 4*q2t(vb)*q2U*q1tL*P4*m1*q2k2 + 4*
     & q2t(vb)*q2U*q1tL*P4*m2*q1k1
      AS1(4,vb) = AS1(4,vb) + 4*q2t(vb)*q2U*q1tL*P4*m2*q2k1 - 4*
     & q2t(vb)*q2U*q1tL*P7*m1*q1k2 - 4*q2t(vb)*q2U*q1tL*P7*m1*q2k2
     &  + 4*q2t(vb)*q2U*q1tL*P7*m2*q1k1 + 4*q2t(vb)*q2U*q1tL*P7*m2*
     & q2k1 - 4*q2t(vb)*q2L*q1tU*P3*m1*q1k2 + 4*q2t(vb)*q2L*q1tU*P3
     & *m2*q1k1 + 4*q2t(vb)*q2L*q1tU*P4*m1*q1k2 - 4*q2t(vb)*q2L*
     & q1tU*P4*m2*q1k1 + 4*q2t(vb)*q2L*q1tU*P7*m1*q1k2 - 4*q2t(vb)*
     & q2L*q1tU*P7*m2*q1k1 - 8*U(vb)*q1L*q1tq2t*q1tk1t*P3*m2 + 8*
     & U(vb)*q1L*q1tq2t*q1tk2t*P3*m1 + 4*U(vb)*q1L*q1tq2t*q2tk1t*
     & P3*m2 + 4*U(vb)*q1L*q1tq2t*q2tk1t*P4*m2 + 4*U(vb)*q1L*
     & q1tq2t*q2tk1t*P7*m2 - 4*U(vb)*q1L*q1tq2t*q2tk2t*P3*m1 - 4*
     & U(vb)*q1L*q1tq2t*q2tk2t*P4*m1 - 4*U(vb)*q1L*q1tq2t*q2tk2t*
     & P7*m1 - 4*U(vb)*q1L*q1tq2t*P3*m1*q2k2 + 4*U(vb)*q1L*q1tq2t*
     & P3*m2*q2k1 - 4*U(vb)*q1L*q1tq2t*P4*m1*q2k2 + 4*U(vb)*q1L*
     & q1tq2t*P4*m2*q2k1 - 4*U(vb)*q1L*q1tq2t*P7*m1*q2k2 + 4*U(vb)*
     & q1L*q1tq2t*P7*m2*q2k1 + 8*U(vb)*q1L*q1tk1t*q2tk1t*P3*m2 - 8
     & *U(vb)*q1L*q1tk1t*q2tk2t*P4*m2
      AS1(4,vb) = AS1(4,vb) - 8*U(vb)*q1L*q1tk1t*q2tk2t*P7*m2 + 4*
     & U(vb)*q1L*q1tk1t*t2*P3*m2 - 4*U(vb)*q1L*q1tk1t*t2*P4*m2 - 4*
     & U(vb)*q1L*q1tk1t*t2*P7*m2 - 8*U(vb)*q1L*q1tk2t*q2tk1t*P3*m1
     &  + 8*U(vb)*q1L*q1tk2t*q2tk2t*P4*m1 + 8*U(vb)*q1L*q1tk2t*
     & q2tk2t*P7*m1 - 4*U(vb)*q1L*q1tk2t*t2*P3*m1 + 4*U(vb)*q1L*
     & q1tk2t*t2*P4*m1 + 4*U(vb)*q1L*q1tk2t*t2*P7*m1 - 4*U(vb)*q2L*
     & q1tq2t*q1tk1t*P3*m2 - 4*U(vb)*q2L*q1tq2t*q1tk1t*P4*m2 - 4*
     & U(vb)*q2L*q1tq2t*q1tk1t*P7*m2 + 4*U(vb)*q2L*q1tq2t*q1tk2t*
     & P3*m1 + 4*U(vb)*q2L*q1tq2t*q1tk2t*P4*m1 + 4*U(vb)*q2L*
     & q1tq2t*q1tk2t*P7*m1 + 8*U(vb)*q2L*q1tq2t*q2tk1t*P4*m2 - 8*
     & U(vb)*q2L*q1tq2t*q2tk2t*P4*m1 + 4*U(vb)*q2L*q1tq2t*P3*m1*
     & q1k2 - 4*U(vb)*q2L*q1tq2t*P3*m2*q1k1 + 4*U(vb)*q2L*q1tq2t*P4
     & *m1*q1k2 - 4*U(vb)*q2L*q1tq2t*P4*m2*q1k1 + 4*U(vb)*q2L*
     & q1tq2t*P7*m1*q1k2 - 4*U(vb)*q2L*q1tq2t*P7*m2*q1k1 + 8*U(vb)*
     & q2L*q1tk1t*q2tk1t*P3*m2 + 8*U(vb)*q2L*q1tk1t*q2tk1t*P7*m2
     &  - 8*U(vb)*q2L*q1tk1t*q2tk2t*P3*m1
      AS1(4,vb) = AS1(4,vb) - 8*U(vb)*q2L*q1tk1t*q2tk2t*P7*m1 - 8*
     & U(vb)*q2L*q1tk2t*q2tk1t*P4*m2 + 8*U(vb)*q2L*q1tk2t*q2tk2t*
     & P4*m1 + 4*U(vb)*q2L*q2tk1t*t1*P3*m2 - 4*U(vb)*q2L*q2tk1t*t1*
     & P4*m2 + 4*U(vb)*q2L*q2tk1t*t1*P7*m2 - 4*U(vb)*q2L*q2tk2t*t1*
     & P3*m1 + 4*U(vb)*q2L*q2tk2t*t1*P4*m1 - 4*U(vb)*q2L*q2tk2t*t1*
     & P7*m1 + 8*U(vb)*q1tq2t*q1tL*P3*m1*q1k2 + 4*U(vb)*q1tq2t*q1tL
     & *P3*m1*q2k2 - 8*U(vb)*q1tq2t*q1tL*P3*m2*q1k1 - 4*U(vb)*q1tq2t
     & *q1tL*P3*m2*q2k1 + 4*U(vb)*q1tq2t*q1tL*P4*m1*q2k2 - 4*U(vb)*
     & q1tq2t*q1tL*P4*m2*q2k1 + 4*U(vb)*q1tq2t*q1tL*P7*m1*q2k2 - 4*
     & U(vb)*q1tq2t*q1tL*P7*m2*q2k1 - 4*U(vb)*q1tq2t*q2tL*P3*m1*
     & q1k2 + 4*U(vb)*q1tq2t*q2tL*P3*m2*q1k1 - 4*U(vb)*q1tq2t*q2tL*
     & P4*m1*q1k2 - 8*U(vb)*q1tq2t*q2tL*P4*m1*q2k2 + 4*U(vb)*q1tq2t*
     & q2tL*P4*m2*q1k1 + 8*U(vb)*q1tq2t*q2tL*P4*m2*q2k1 - 4*U(vb)*
     & q1tq2t*q2tL*P7*m1*q1k2 + 4*U(vb)*q1tq2t*q2tL*P7*m2*q1k1 - 4*
     & U(vb)*q1tk1t*q2tL*t1*P3*m2 - 4*U(vb)*q1tk1t*q2tL*t2*P4*m2 - 
     & 8*U(vb)*q1tk1t*q2tL*P3*m1*q2k2
      AS1(4,vb) = AS1(4,vb) - 4*U(vb)*q1tk1t*q2tL*P3*m2*q1q2 + 8*
     & U(vb)*q1tk1t*q2tL*P3*m2*q1k1 + 8*U(vb)*q1tk1t*q2tL*P3*m2*
     & q2k1 - 4*U(vb)*q1tk1t*q2tL*P4*m2*q1q2 + 8*U(vb)*q1tk1t*q2tL*
     & P4*m2*q2k2 - 8*U(vb)*q1tk1t*q2tL*P7*m1*q2k2 - 4*U(vb)*q1tk1t*
     & q2tL*P7*m2*q1q2 + 8*U(vb)*q1tk1t*q2tL*P7*m2*q2k1 + 4*U(vb)*
     & q1tk2t*q2tL*t1*P3*m1 + 4*U(vb)*q1tk2t*q2tL*t2*P4*m1 + 4*
     & U(vb)*q1tk2t*q2tL*P3*m1*q1q2 - 8*U(vb)*q1tk2t*q2tL*P3*m1*
     & q1k1 + 4*U(vb)*q1tk2t*q2tL*P4*m1*q1q2 - 8*U(vb)*q1tk2t*q2tL*
     & P4*m2*q2k1 + 4*U(vb)*q1tk2t*q2tL*P7*m1*q1q2 + 4*U(vb)*q1tL*
     & q2tk1t*t1*P3*m2 + 4*U(vb)*q1tL*q2tk1t*t2*P4*m2 - 8*U(vb)*
     & q1tL*q2tk1t*P3*m1*q1k2 + 4*U(vb)*q1tL*q2tk1t*P3*m2*q1q2 + 4*
     & U(vb)*q1tL*q2tk1t*P4*m2*q1q2 - 8*U(vb)*q1tL*q2tk1t*P4*m2*
     & q2k2 + 4*U(vb)*q1tL*q2tk1t*P7*m2*q1q2 - 4*U(vb)*q1tL*q2tk2t*
     & t1*P3*m1 - 4*U(vb)*q1tL*q2tk2t*t2*P4*m1 - 4*U(vb)*q1tL*
     & q2tk2t*P3*m1*q1q2 + 8*U(vb)*q1tL*q2tk2t*P3*m1*q1k1 - 4*U(vb)*
     & q1tL*q2tk2t*P4*m1*q1q2
      AS1(4,vb) = AS1(4,vb) + 8*U(vb)*q1tL*q2tk2t*P4*m1*q1k2 + 8*
     & U(vb)*q1tL*q2tk2t*P4*m1*q2k2 - 8*U(vb)*q1tL*q2tk2t*P4*m2*
     & q1k1 - 4*U(vb)*q1tL*q2tk2t*P7*m1*q1q2 + 8*U(vb)*q1tL*q2tk2t*
     & P7*m1*q1k2 - 8*U(vb)*q1tL*q2tk2t*P7*m2*q1k1 - 4*U(vb)*q1tL*t2
     & *P3*m1*q1k2 + 4*U(vb)*q1tL*t2*P3*m2*q1k1 + 4*U(vb)*q1tL*t2*P4*
     & m1*q1k2 - 4*U(vb)*q1tL*t2*P4*m2*q1k1 + 4*U(vb)*q1tL*t2*P7*m1*
     & q1k2 - 4*U(vb)*q1tL*t2*P7*m2*q1k1 - 4*U(vb)*q2tL*t1*P3*m1*q2k2
     &  + 4*U(vb)*q2tL*t1*P3*m2*q2k1 + 4*U(vb)*q2tL*t1*P4*m1*q2k2 - 4
     & *U(vb)*q2tL*t1*P4*m2*q2k1 - 4*U(vb)*q2tL*t1*P7*m1*q2k2 + 4*
     & U(vb)*q2tL*t1*P7*m2*q2k1 - 4*L(vb)*q2U*q1tq2t*q1tk1t*P3*m2
     &  + 4*L(vb)*q2U*q1tq2t*q1tk1t*P4*m2 + 4*L(vb)*q2U*q1tq2t*
     & q1tk1t*P7*m2 + 4*L(vb)*q2U*q1tq2t*q1tk2t*P3*m1 - 4*L(vb)*
     & q2U*q1tq2t*q1tk2t*P4*m1 - 4*L(vb)*q2U*q1tq2t*q1tk2t*P7*m1
     &  + 4*L(vb)*q2U*q1tq2t*q2tk1t*P3*m2 - 4*L(vb)*q2U*q1tq2t*
     & q2tk1t*P4*m2 + 4*L(vb)*q2U*q1tq2t*q2tk1t*P7*m2 - 4*L(vb)*
     & q2U*q1tq2t*q2tk2t*P3*m1
      AS1(4,vb) = AS1(4,vb) + 4*L(vb)*q2U*q1tq2t*q2tk2t*P4*m1 - 4*
     & L(vb)*q2U*q1tq2t*q2tk2t*P7*m1 - 4*L(vb)*q2U*q1tq2t*P3*m1*
     & q1k2 - 4*L(vb)*q2U*q1tq2t*P3*m1*q2k2 + 4*L(vb)*q2U*q1tq2t*P3
     & *m2*q1k1 + 4*L(vb)*q2U*q1tq2t*P3*m2*q2k1 - 4*L(vb)*q2U*
     & q1tq2t*P4*m1*q1k2 - 4*L(vb)*q2U*q1tq2t*P4*m1*q2k2 + 4*L(vb)*
     & q2U*q1tq2t*P4*m2*q1k1 + 4*L(vb)*q2U*q1tq2t*P4*m2*q2k1 - 4*
     & L(vb)*q2U*q1tq2t*P7*m1*q1k2 - 4*L(vb)*q2U*q1tq2t*P7*m1*q2k2
     &  + 4*L(vb)*q2U*q1tq2t*P7*m2*q1k1 + 4*L(vb)*q2U*q1tq2t*P7*m2*
     & q2k1 - 8*L(vb)*q2U*q1tk1t*q2tk1t*P7*m2 + 8*L(vb)*q2U*q1tk1t
     & *q2tk2t*P3*m1 - 8*L(vb)*q2U*q1tk1t*q2tk2t*P4*m2 + 8*L(vb)*
     & q2U*q1tk1t*q2tk2t*P7*m1 - 8*L(vb)*q2U*q1tk1t*q2tk2t*P7*m2
     &  + 4*L(vb)*q2U*q1tk1t*t2*P3*m2 - 4*L(vb)*q2U*q1tk1t*t2*P4*m2
     &  - 4*L(vb)*q2U*q1tk1t*t2*P7*m2 - 8*L(vb)*q2U*q1tk2t*q2tk1t*
     & P3*m1 + 8*L(vb)*q2U*q1tk2t*q2tk1t*P4*m2 + 8*L(vb)*q2U*
     & q1tk2t*q2tk2t*P7*m1 - 4*L(vb)*q2U*q1tk2t*t2*P3*m1 + 4*L(vb)*
     & q2U*q1tk2t*t2*P4*m1
      AS1(4,vb) = AS1(4,vb) + 4*L(vb)*q2U*q1tk2t*t2*P7*m1 - 4*L(vb)
     & *q2U*q2tk1t*t1*P3*m2 + 4*L(vb)*q2U*q2tk1t*t1*P4*m2 - 4*L(vb)
     & *q2U*q2tk1t*t1*P7*m2 + 4*L(vb)*q2U*q2tk2t*t1*P3*m1 - 4*L(vb)
     & *q2U*q2tk2t*t1*P4*m1 + 4*L(vb)*q2U*q2tk2t*t1*P7*m1 - 8*L(vb)
     & *q1tq2t*q1tU*P3*m1*q1k2 - 4*L(vb)*q1tq2t*q1tU*P3*m1*q2k2 + 8
     & *L(vb)*q1tq2t*q1tU*P3*m2*q1k1 + 4*L(vb)*q1tq2t*q1tU*P3*m2*
     & q2k1 - 4*L(vb)*q1tq2t*q1tU*P4*m1*q2k2 + 4*L(vb)*q1tq2t*q1tU*
     & P4*m2*q2k1 - 4*L(vb)*q1tq2t*q1tU*P7*m1*q2k2 + 4*L(vb)*q1tq2t*
     & q1tU*P7*m2*q2k1 + 4*L(vb)*q1tq2t*q2tU*P3*m1*q1k2 - 4*L(vb)*
     & q1tq2t*q2tU*P3*m2*q1k1 + 4*L(vb)*q1tq2t*q2tU*P4*m1*q1k2 + 8*
     & L(vb)*q1tq2t*q2tU*P4*m1*q2k2 - 4*L(vb)*q1tq2t*q2tU*P4*m2*
     & q1k1 - 8*L(vb)*q1tq2t*q2tU*P4*m2*q2k1 + 4*L(vb)*q1tq2t*q2tU*
     & P7*m1*q1k2 - 4*L(vb)*q1tq2t*q2tU*P7*m2*q1k1 + 4*L(vb)*q1tk1t*
     & q2tU*t1*P3*m2 + 4*L(vb)*q1tk1t*q2tU*t2*P4*m2 + 8*L(vb)*
     & q1tk1t*q2tU*P3*m1*q2k2 + 4*L(vb)*q1tk1t*q2tU*P3*m2*q1q2 - 8*
     & L(vb)*q1tk1t*q2tU*P3*m2*q1k1
      AS1(4,vb) = AS1(4,vb) - 8*L(vb)*q1tk1t*q2tU*P3*m2*q2k1 + 4*
     & L(vb)*q1tk1t*q2tU*P4*m2*q1q2 - 8*L(vb)*q1tk1t*q2tU*P4*m2*
     & q2k2 + 8*L(vb)*q1tk1t*q2tU*P7*m1*q2k2 + 4*L(vb)*q1tk1t*q2tU*
     & P7*m2*q1q2 - 8*L(vb)*q1tk1t*q2tU*P7*m2*q2k1 - 4*L(vb)*q1tk2t*
     & q2tU*t1*P3*m1 - 4*L(vb)*q1tk2t*q2tU*t2*P4*m1 - 4*L(vb)*
     & q1tk2t*q2tU*P3*m1*q1q2 + 8*L(vb)*q1tk2t*q2tU*P3*m1*q1k1 - 4*
     & L(vb)*q1tk2t*q2tU*P4*m1*q1q2 + 8*L(vb)*q1tk2t*q2tU*P4*m2*
     & q2k1 - 4*L(vb)*q1tk2t*q2tU*P7*m1*q1q2 - 4*L(vb)*q1tU*q2tk1t*
     & t1*P3*m2 - 4*L(vb)*q1tU*q2tk1t*t2*P4*m2 + 8*L(vb)*q1tU*
     & q2tk1t*P3*m1*q1k2 - 4*L(vb)*q1tU*q2tk1t*P3*m2*q1q2 - 4*L(vb)*
     & q1tU*q2tk1t*P4*m2*q1q2 + 8*L(vb)*q1tU*q2tk1t*P4*m2*q2k2 - 4*
     & L(vb)*q1tU*q2tk1t*P7*m2*q1q2 + 4*L(vb)*q1tU*q2tk2t*t1*P3*m1
     &  + 4*L(vb)*q1tU*q2tk2t*t2*P4*m1 + 4*L(vb)*q1tU*q2tk2t*P3*m1*
     & q1q2 - 8*L(vb)*q1tU*q2tk2t*P3*m1*q1k1 + 4*L(vb)*q1tU*q2tk2t*
     & P4*m1*q1q2 - 8*L(vb)*q1tU*q2tk2t*P4*m1*q1k2 - 8*L(vb)*q1tU*
     & q2tk2t*P4*m1*q2k2
      AS1(4,vb) = AS1(4,vb) + 8*L(vb)*q1tU*q2tk2t*P4*m2*q1k1 + 4*
     & L(vb)*q1tU*q2tk2t*P7*m1*q1q2 - 8*L(vb)*q1tU*q2tk2t*P7*m1*
     & q1k2 + 8*L(vb)*q1tU*q2tk2t*P7*m2*q1k1 + 4*L(vb)*q1tU*t2*P3*m1
     & *q1k2 - 4*L(vb)*q1tU*t2*P3*m2*q1k1 - 4*L(vb)*q1tU*t2*P4*m1*
     & q1k2 + 4*L(vb)*q1tU*t2*P4*m2*q1k1 - 4*L(vb)*q1tU*t2*P7*m1*q1k2
     &  + 4*L(vb)*q1tU*t2*P7*m2*q1k1 + 4*L(vb)*q2tU*t1*P3*m1*q2k2 - 4
     & *L(vb)*q2tU*t1*P3*m2*q2k1 - 4*L(vb)*q2tU*t1*P4*m1*q2k2 + 4*
     & L(vb)*q2tU*t1*P4*m2*q2k1 + 4*L(vb)*q2tU*t1*P7*m1*q2k2 - 4*
     & L(vb)*q2tU*t1*P7*m2*q2k1

      AS2(1,vb) =
     & 4*k1(vb)*t1*t2*P5*m2 + 4*k1(vb)*t1*t2*P6*m2 - 4*k1(vb)*t1*t2*P8*
     & m2 + 4*k1(vb)*q1tq2t*q1tk1t*P6*m2 - 4*k1(vb)*q1tq2t*q1tk2t*
     & P6*m1 - 16*k1(vb)*q1tq2t*q1tk2t*P6*m2 - 12*k1(vb)*q1tq2t*
     & q2tk1t*P5*m2 - 4*k1(vb)*q1tq2t*q2tk2t*P5*m1 - 4*k1(vb)*
     & q1tq2t*t1*P6*m2 - 4*k1(vb)*q1tq2t*t2*P5*m2 - 4*k1(vb)*q1tq2t*
     & P5*m1*q2k2 + 4*k1(vb)*q1tq2t*P5*m2*q1q2 - 4*k1(vb)*q1tq2t*P5*
     & m2*q2k1 - 4*k1(vb)*q1tq2t*P6*m1*q1k2 + 4*k1(vb)*q1tq2t*P6*m2*
     & q1q2 + 4*k1(vb)*q1tq2t*P6*m2*q1k1 - 8*k1(vb)*q1tq2t*P6*m2*q1k2
     &  + 4*k1(vb)*q1tq2t*P8*m2*q1q2 + 4*k1(vb)*q1tq2t**2*P5*m2 + 4*
     & k1(vb)*q1tq2t**2*P6*m2 + 4*k1(vb)*q1tq2t**2*P8*m2 + 8*k1(vb)*
     & q1tk1t*q2tk1t*P5*m2 + 4*k1(vb)*q1tk1t*t2*P5*m2 + 8*k1(vb)*
     & q1tk2t*q2tk1t*P5*m1 - 8*k1(vb)*q1tk2t*q2tk1t*P6*m2 - 16*
     & k1(vb)*q1tk2t*q2tk1t*P8*m2 + 8*k1(vb)*q1tk2t*q2tk2t*P6*m1 + 
     & 16*k1(vb)*q1tk2t*q2tk2t*P6*m2 + 4*k1(vb)*q1tk2t*t2*P5*m1 + 8*
     & k1(vb)*q1tk2t*t2*P6*m2
      AS2(1,vb) = AS2(1,vb) - 8*k1(vb)*q1tk2t*t2*P8*m2 + 8*k1(vb)*
     & q2tk1t*t1*P5*m2 - 4*k1(vb)*q2tk1t*t1*P6*m2 - 8*k1(vb)*q2tk1t*
     & t1*P8*m2 + 4*k1(vb)*q2tk2t*t1*P6*m1 + 8*k1(vb)*q2tk2t*t1*P6*m2
     &  - 4*k2(vb)*t1*t2*P5*m1 - 4*k2(vb)*t1*t2*P6*m1 + 4*k2(vb)*t1*t2*
     & P8*m1 + 4*k2(vb)*q1tq2t*q1tk1t*P6*m2 + 12*k2(vb)*q1tq2t*
     & q1tk2t*P6*m1 + 16*k2(vb)*q1tq2t*q2tk1t*P5*m1 + 4*k2(vb)*
     & q1tq2t*q2tk1t*P5*m2 - 4*k2(vb)*q1tq2t*q2tk2t*P5*m1 + 4*
     & k2(vb)*q1tq2t*t1*P6*m1 + 4*k2(vb)*q1tq2t*t2*P5*m1 - 4*k2(vb)*
     & q1tq2t*P5*m1*q1q2 + 8*k2(vb)*q1tq2t*P5*m1*q2k1 - 4*k2(vb)*
     & q1tq2t*P5*m1*q2k2 + 4*k2(vb)*q1tq2t*P5*m2*q2k1 - 4*k2(vb)*
     & q1tq2t*P6*m1*q1q2 + 4*k2(vb)*q1tq2t*P6*m1*q1k2 + 4*k2(vb)*
     & q1tq2t*P6*m2*q1k1 - 4*k2(vb)*q1tq2t*P8*m1*q1q2 - 4*k2(vb)*
     & q1tq2t**2*P5*m1 - 4*k2(vb)*q1tq2t**2*P6*m1 - 4*k2(vb)*
     & q1tq2t**2*P8*m1 - 16*k2(vb)*q1tk1t*q2tk1t*P5*m1 - 8*k2(vb)*
     & q1tk1t*q2tk1t*P5*m2 - 8*k2(vb)*q1tk1t*t2*P5*m1 - 4*k2(vb)*
     & q1tk1t*t2*P5*m2
      AS2(1,vb) = AS2(1,vb) + 8*k2(vb)*q1tk2t*q2tk1t*P5*m1 - 8*
     & k2(vb)*q1tk2t*q2tk1t*P6*m2 + 16*k2(vb)*q1tk2t*q2tk1t*P8*m1
     &  - 8*k2(vb)*q1tk2t*q2tk2t*P6*m1 + 4*k2(vb)*q1tk2t*t2*P5*m1 - 
     & 8*k2(vb)*q1tk2t*t2*P6*m1 + 8*k2(vb)*q1tk2t*t2*P8*m1 - 8*k2(vb)
     & *q2tk1t*t1*P5*m1 - 4*k2(vb)*q2tk1t*t1*P6*m2 + 8*k2(vb)*q2tk1t
     & *t1*P8*m1 - 4*k2(vb)*q2tk2t*t1*P6*m1 + 4*kw(vb)*q1tq2t*q1tk1t
     & *P6*m2 - 4*kw(vb)*q1tq2t*q1tk2t*P6*m1 + 4*kw(vb)*q1tq2t*
     & q2tk1t*P5*m2 - 4*kw(vb)*q1tq2t*q2tk2t*P5*m1 - 4*kw(vb)*
     & q1tq2t*P5*m1*q2k2 + 4*kw(vb)*q1tq2t*P5*m2*q2k1 - 4*kw(vb)*
     & q1tq2t*P6*m1*q1k2 + 4*kw(vb)*q1tq2t*P6*m2*q1k1 - 8*kw(vb)*
     & q1tk1t*q2tk1t*P5*m2 - 4*kw(vb)*q1tk1t*t2*P5*m2 + 8*kw(vb)*
     & q1tk2t*q2tk1t*P5*m1 - 8*kw(vb)*q1tk2t*q2tk1t*P6*m2 + 8*
     & kw(vb)*q1tk2t*q2tk2t*P6*m1 + 4*kw(vb)*q1tk2t*t2*P5*m1 - 4*
     & kw(vb)*q2tk1t*t1*P6*m2 + 4*kw(vb)*q2tk2t*t1*P6*m1 + 4*q1(vb)*
     & q1tq2t*q1tk1t*P6*m2 - 4*q1(vb)*q1tq2t*q1tk2t*P6*m1 - 4*
     & q1(vb)*q1tq2t*q2tk1t*P6*m2
      AS2(1,vb) = AS2(1,vb) - 4*q1(vb)*q1tq2t*q2tk1t*P8*m2 + 4*
     & q1(vb)*q1tq2t*q2tk2t*P6*m1 + 4*q1(vb)*q1tq2t*q2tk2t*P8*m1 + 
     & 4*q1(vb)*q1tq2t*P6*m1*q1k2 + 4*q1(vb)*q1tq2t*P6*m1*q2k2 - 4*
     & q1(vb)*q1tq2t*P6*m2*q1k1 - 4*q1(vb)*q1tq2t*P6*m2*q2k1 + 4*
     & q1(vb)*q1tq2t*P8*m1*q2k2 - 4*q1(vb)*q1tq2t*P8*m2*q2k1 + 8*
     & q1(vb)*q1tk1t*q2tk1t*P8*m2 - 8*q1(vb)*q1tk1t*q2tk2t*P6*m2 - 
     & 4*q1(vb)*q1tk1t*t2*P6*m2 + 4*q1(vb)*q1tk1t*t2*P8*m2 + 8*q1(vb)
     & *q1tk2t*q2tk1t*P6*m2 - 8*q1(vb)*q1tk2t*q2tk1t*P8*m1 + 4*
     & q1(vb)*q1tk2t*t2*P6*m1 - 4*q1(vb)*q1tk2t*t2*P8*m1 + 4*q1(vb)*
     & q2tk1t*t1*P6*m2 - 4*q1(vb)*q2tk2t*t1*P6*m1 - 4*q2(vb)*q1tq2t*
     & q1tk1t*P5*m2 - 4*q2(vb)*q1tq2t*q1tk1t*P8*m2 + 4*q2(vb)*
     & q1tq2t*q1tk2t*P5*m1 + 4*q2(vb)*q1tq2t*q1tk2t*P8*m1 + 4*
     & q2(vb)*q1tq2t*q2tk1t*P5*m2 - 4*q2(vb)*q1tq2t*q2tk2t*P5*m1 + 
     & 4*q2(vb)*q1tq2t*P5*m1*q1k2 + 4*q2(vb)*q1tq2t*P5*m1*q2k2 - 4*
     & q2(vb)*q1tq2t*P5*m2*q1k1 - 4*q2(vb)*q1tq2t*P5*m2*q2k1 + 4*
     & q2(vb)*q1tq2t*P8*m1*q1k2
      AS2(1,vb) = AS2(1,vb) - 4*q2(vb)*q1tq2t*P8*m2*q1k1 + 8*q2(vb)*
     & q1tk1t*q2tk2t*P5*m1 + 4*q2(vb)*q1tk1t*t2*P5*m2 - 8*q2(vb)*
     & q1tk2t*q2tk1t*P5*m1 + 8*q2(vb)*q1tk2t*q2tk1t*P8*m2 - 8*
     & q2(vb)*q1tk2t*q2tk2t*P8*m1 - 4*q2(vb)*q1tk2t*t2*P5*m1 - 4*
     & q2(vb)*q2tk1t*t1*P5*m2 + 4*q2(vb)*q2tk1t*t1*P8*m2 + 4*q2(vb)*
     & q2tk2t*t1*P5*m1 - 4*q2(vb)*q2tk2t*t1*P8*m1 + 4*q1t(vb)*t2*P5*
     & m1*q1k2 - 4*q1t(vb)*t2*P5*m2*q1k1 + 4*q1t(vb)*t2*P6*m1*q1k2 - 4*
     & q1t(vb)*t2*P6*m2*q1k1 - 4*q1t(vb)*t2*P8*m1*q1k2 + 4*q1t(vb)*t2*
     & P8*m2*q1k1 + 4*q1t(vb)*q1tq2t*P5*m1*q2k2 - 4*q1t(vb)*q1tq2t*P5
     & *m2*q2k1 - 8*q1t(vb)*q1tq2t*P6*m1*q1k2 - 4*q1t(vb)*q1tq2t*P6*
     & m1*q2k2 + 8*q1t(vb)*q1tq2t*P6*m2*q1k1 + 4*q1t(vb)*q1tq2t*P6*m2
     & *q2k1 - 4*q1t(vb)*q1tq2t*P8*m1*q2k2 + 4*q1t(vb)*q1tq2t*P8*m2*
     & q2k1 - 4*q1t(vb)*q2tk1t*t1*P6*m2 + 4*q1t(vb)*q2tk1t*t2*P5*m2
     &  + 8*q1t(vb)*q2tk1t*P5*m1*q1k2 + 4*q1t(vb)*q2tk1t*P5*m2*q1q2
     &  - 8*q1t(vb)*q2tk1t*P5*m2*q1k1 - 8*q1t(vb)*q2tk1t*P5*m2*q2k1
     &  - 4*q1t(vb)*q2tk1t*P6*m2*q1q2
      AS2(1,vb) = AS2(1,vb) + 8*q1t(vb)*q2tk1t*P6*m2*q1k2 - 8*
     & q1t(vb)*q2tk1t*P8*m1*q1k2 - 4*q1t(vb)*q2tk1t*P8*m2*q1q2 + 8*
     & q1t(vb)*q2tk1t*P8*m2*q1k1 + 4*q1t(vb)*q2tk2t*t1*P6*m1 - 4*
     & q1t(vb)*q2tk2t*t2*P5*m1 - 4*q1t(vb)*q2tk2t*P5*m1*q1q2 + 8*
     & q1t(vb)*q2tk2t*P5*m1*q2k1 + 4*q1t(vb)*q2tk2t*P6*m1*q1q2 - 8*
     & q1t(vb)*q2tk2t*P6*m2*q1k1 + 4*q1t(vb)*q2tk2t*P8*m1*q1q2 + 4*
     & q2t(vb)*t1*P5*m1*q2k2 - 4*q2t(vb)*t1*P5*m2*q2k1 + 4*q2t(vb)*t1*
     & P6*m1*q2k2 - 4*q2t(vb)*t1*P6*m2*q2k1 - 4*q2t(vb)*t1*P8*m1*q2k2
     &  + 4*q2t(vb)*t1*P8*m2*q2k1 - 4*q2t(vb)*q1tq2t*P5*m1*q1k2 - 8*
     & q2t(vb)*q1tq2t*P5*m1*q2k2 + 4*q2t(vb)*q1tq2t*P5*m2*q1k1 + 8*
     & q2t(vb)*q1tq2t*P5*m2*q2k1 + 4*q2t(vb)*q1tq2t*P6*m1*q1k2 - 4*
     & q2t(vb)*q1tq2t*P6*m2*q1k1 - 4*q2t(vb)*q1tq2t*P8*m1*q1k2 + 4*
     & q2t(vb)*q1tq2t*P8*m2*q1k1 + 4*q2t(vb)*q1tk1t*t1*P6*m2 - 4*
     & q2t(vb)*q1tk1t*t2*P5*m2 + 8*q2t(vb)*q1tk1t*P5*m1*q2k2 - 4*
     & q2t(vb)*q1tk1t*P5*m2*q1q2 + 4*q2t(vb)*q1tk1t*P6*m2*q1q2 - 8*
     & q2t(vb)*q1tk1t*P6*m2*q1k2
      AS2(1,vb) = AS2(1,vb) - 4*q2t(vb)*q1tk1t*P8*m2*q1q2 - 4*
     & q2t(vb)*q1tk2t*t1*P6*m1 + 4*q2t(vb)*q1tk2t*t2*P5*m1 + 4*
     & q2t(vb)*q1tk2t*P5*m1*q1q2 - 8*q2t(vb)*q1tk2t*P5*m1*q2k1 - 4*
     & q2t(vb)*q1tk2t*P6*m1*q1q2 + 8*q2t(vb)*q1tk2t*P6*m1*q1k2 + 8*
     & q2t(vb)*q1tk2t*P6*m1*q2k2 - 8*q2t(vb)*q1tk2t*P6*m2*q2k1 + 4*
     & q2t(vb)*q1tk2t*P8*m1*q1q2 - 8*q2t(vb)*q1tk2t*P8*m1*q2k2 + 8*
     & q2t(vb)*q1tk2t*P8*m2*q2k1

      AS2(2,vb) =
     & 4*k1(vb)*q2U*q1tq2t*q1tk2t*P5 - 4*k1(vb)*q2U*q1tq2t*q1tk2t
     & *P6 + 4*k1(vb)*q2U*q1tq2t*q1tk2t*P8 - 4*k1(vb)*q2U*q1tq2t*
     & q2tk2t*P5 + 4*k1(vb)*q2U*q1tq2t*q2tk2t*P6 + 4*k1(vb)*q2U*
     & q1tq2t*q2tk2t*P8 - 4*k1(vb)*q2U*q1tq2t*P5*m1*m2 + 4*k1(vb)*
     & q2U*q1tq2t*P5*q1k2 + 4*k1(vb)*q2U*q1tq2t*P5*q2k2 - 4*k1(vb)*
     & q2U*q1tq2t*P5*k1k2 + 4*k1(vb)*q2U*q1tq2t*P6*m1*m2 + 4*k1(vb)
     & *q2U*q1tq2t*P6*q1k2 + 4*k1(vb)*q2U*q1tq2t*P6*q2k2 + 4*k1(vb)
     & *q2U*q1tq2t*P6*k1k2 + 4*k1(vb)*q2U*q1tq2t*P8*q1k2 + 4*k1(vb)
     & *q2U*q1tq2t*P8*q2k2 + 4*k1(vb)*q2U*q1tk1t*q2tk2t*P5 + 4*
     & k1(vb)*q2U*q1tk1t*q2tk2t*P6 - 4*k1(vb)*q2U*q1tk2t*q2tk1t*
     & P5 - 4*k1(vb)*q2U*q1tk2t*q2tk1t*P6 - 8*k1(vb)*q2U*q1tk2t*
     & q2tk1t*P8 - 8*k1(vb)*q2U*q1tk2t*q2tk2t*P8 - 4*k1(vb)*q2U*
     & q1tk2t*t2*P5 + 4*k1(vb)*q2U*q1tk2t*t2*P6 - 4*k1(vb)*q2U*
     & q1tk2t*t2*P8 + 4*k1(vb)*q2U*q2tk2t*t1*P5 - 4*k1(vb)*q2U*
     & q2tk2t*t1*P6
      AS2(2,vb) = AS2(2,vb) - 4*k1(vb)*q2U*q2tk2t*t1*P8 + 4*k1(vb)*
     & q1tq2t*q1tU*P5*q2k2 + 4*k1(vb)*q1tq2t*q1tU*P6*m1*m2 + 8*
     & k1(vb)*q1tq2t*q1tU*P6*q1k2 + 4*k1(vb)*q1tq2t*q1tU*P6*q2k2 + 
     & 4*k1(vb)*q1tq2t*q1tU*P6*k1k2 + 4*k1(vb)*q1tq2t*q1tU*P8*q2k2
     &  + 4*k1(vb)*q1tq2t*q2tU*P5*m1*m2 - 4*k1(vb)*q1tq2t*q2tU*P5*
     & q1k2 - 8*k1(vb)*q1tq2t*q2tU*P5*q2k2 + 4*k1(vb)*q1tq2t*q2tU*
     & P5*k1k2 - 4*k1(vb)*q1tq2t*q2tU*P6*q1k2 - 4*k1(vb)*q1tq2t*
     & q2tU*P8*q1k2 + 4*k1(vb)*q1tk1t*q2tU*P5*q2k2 - 4*k1(vb)*
     & q1tk1t*q2tU*P6*q1k2 + 4*k1(vb)*q1tk2t*q2tU*t1*P6 + 4*k1(vb)*
     & q1tk2t*q2tU*t2*P5 + 4*k1(vb)*q1tk2t*q2tU*P5*q1q2 - 4*k1(vb)*
     & q1tk2t*q2tU*P5*q2k1 - 8*k1(vb)*q1tk2t*q2tU*P6*m1*m2 + 4*
     & k1(vb)*q1tk2t*q2tU*P6*q1q2 + 4*k1(vb)*q1tk2t*q2tU*P6*q1k1 - 
     & 8*k1(vb)*q1tk2t*q2tU*P6*q1k2 - 8*k1(vb)*q1tk2t*q2tU*P6*q2k2
     &  - 8*k1(vb)*q1tk2t*q2tU*P6*k1k2 + 4*k1(vb)*q1tk2t*q2tU*P8*
     & q1q2 - 8*k1(vb)*q1tk2t*q2tU*P8*q2k2 - 8*k1(vb)*q1tU*q2tk1t*
     & P5*m1*m2
      AS2(2,vb) = AS2(2,vb) + 8*k1(vb)*q1tU*q2tk1t*P5*q1k2 + 4*
     & k1(vb)*q1tU*q2tk1t*P5*q2k2 - 8*k1(vb)*q1tU*q2tk1t*P5*k1k2 + 
     & 4*k1(vb)*q1tU*q2tk1t*P6*q1k2 + 8*k1(vb)*q1tU*q2tk1t*P8*q1k2
     &  - 4*k1(vb)*q1tU*q2tk2t*t1*P6 - 4*k1(vb)*q1tU*q2tk2t*t2*P5
     &  - 4*k1(vb)*q1tU*q2tk2t*P5*q1q2 + 4*k1(vb)*q1tU*q2tk2t*P5*
     & q2k1 - 4*k1(vb)*q1tU*q2tk2t*P6*q1q2 - 4*k1(vb)*q1tU*q2tk2t*
     & P6*q1k1 - 4*k1(vb)*q1tU*q2tk2t*P8*q1q2 - 4*k1(vb)*q1tU*t2*P5*
     & m1*m2 + 4*k1(vb)*q1tU*t2*P5*q1k2 - 4*k1(vb)*q1tU*t2*P5*k1k2 - 
     & 4*k1(vb)*q1tU*t2*P6*q1k2 + 4*k1(vb)*q1tU*t2*P8*q1k2 + 4*k1(vb)
     & *q2tU*t1*P5*q2k2 - 4*k1(vb)*q2tU*t1*P6*m1*m2 - 4*k1(vb)*q2tU*
     & t1*P6*q2k2 - 4*k1(vb)*q2tU*t1*P6*k1k2 - 4*k1(vb)*q2tU*t1*P8*
     & q2k2 - 4*k2(vb)*q2U*q1tq2t*q1tk1t*P5 + 4*k2(vb)*q2U*q1tq2t*
     & q1tk1t*P6 - 4*k2(vb)*q2U*q1tq2t*q1tk1t*P8 + 4*k2(vb)*q2U*
     & q1tq2t*q2tk1t*P5 - 4*k2(vb)*q2U*q1tq2t*q2tk1t*P6 - 4*k2(vb)
     & *q2U*q1tq2t*q2tk1t*P8 - 4*k2(vb)*q2U*q1tq2t*P5*m1*m2 - 4*
     & k2(vb)*q2U*q1tq2t*P5*q1k1
      AS2(2,vb) = AS2(2,vb) - 4*k2(vb)*q2U*q1tq2t*P5*q2k1 - 4*
     & k2(vb)*q2U*q1tq2t*P5*k1k2 + 4*k2(vb)*q2U*q1tq2t*P6*m1*m2 - 4
     & *k2(vb)*q2U*q1tq2t*P6*q1k1 - 4*k2(vb)*q2U*q1tq2t*P6*q2k1 + 4
     & *k2(vb)*q2U*q1tq2t*P6*k1k2 - 4*k2(vb)*q2U*q1tq2t*P8*q1k1 - 4
     & *k2(vb)*q2U*q1tq2t*P8*q2k1 + 8*k2(vb)*q2U*q1tk1t*q2tk1t*P8
     &  - 4*k2(vb)*q2U*q1tk1t*q2tk2t*P5 - 4*k2(vb)*q2U*q1tk1t*
     & q2tk2t*P6 + 4*k2(vb)*q2U*q1tk1t*t2*P5 - 4*k2(vb)*q2U*q1tk1t
     & *t2*P6 + 4*k2(vb)*q2U*q1tk1t*t2*P8 + 4*k2(vb)*q2U*q1tk2t*
     & q2tk1t*P5 + 4*k2(vb)*q2U*q1tk2t*q2tk1t*P6 + 8*k2(vb)*q2U*
     & q1tk2t*q2tk1t*P8 - 4*k2(vb)*q2U*q2tk1t*t1*P5 + 4*k2(vb)*q2U
     & *q2tk1t*t1*P6 + 4*k2(vb)*q2U*q2tk1t*t1*P8 - 4*k2(vb)*q1tq2t*
     & q1tU*P5*q2k1 + 4*k2(vb)*q1tq2t*q1tU*P6*m1*m2 - 8*k2(vb)*
     & q1tq2t*q1tU*P6*q1k1 - 4*k2(vb)*q1tq2t*q1tU*P6*q2k1 + 4*
     & k2(vb)*q1tq2t*q1tU*P6*k1k2 - 4*k2(vb)*q1tq2t*q1tU*P8*q2k1 + 
     & 4*k2(vb)*q1tq2t*q2tU*P5*m1*m2 + 4*k2(vb)*q1tq2t*q2tU*P5*q1k1
     &  + 8*k2(vb)*q1tq2t*q2tU*P5*q2k1
      AS2(2,vb) = AS2(2,vb) + 4*k2(vb)*q1tq2t*q2tU*P5*k1k2 + 4*
     & k2(vb)*q1tq2t*q2tU*P6*q1k1 + 4*k2(vb)*q1tq2t*q2tU*P8*q1k1 - 
     & 4*k2(vb)*q1tk1t*q2tU*t1*P6 - 4*k2(vb)*q1tk1t*q2tU*t2*P5 - 4*
     & k2(vb)*q1tk1t*q2tU*P5*q1q2 - 4*k2(vb)*q1tk1t*q2tU*P5*q2k2 - 
     & 4*k2(vb)*q1tk1t*q2tU*P6*q1q2 + 4*k2(vb)*q1tk1t*q2tU*P6*q1k2
     &  - 4*k2(vb)*q1tk1t*q2tU*P8*q1q2 + 4*k2(vb)*q1tk2t*q2tU*P5*
     & q2k1 - 8*k2(vb)*q1tk2t*q2tU*P6*m1*m2 + 4*k2(vb)*q1tk2t*q2tU*
     & P6*q1k1 + 8*k2(vb)*q1tk2t*q2tU*P6*q2k1 - 8*k2(vb)*q1tk2t*
     & q2tU*P6*k1k2 + 8*k2(vb)*q1tk2t*q2tU*P8*q2k1 + 4*k2(vb)*q1tU*
     & q2tk1t*t1*P6 + 4*k2(vb)*q1tU*q2tk1t*t2*P5 - 8*k2(vb)*q1tU*
     & q2tk1t*P5*m1*m2 + 4*k2(vb)*q1tU*q2tk1t*P5*q1q2 - 8*k2(vb)*
     & q1tU*q2tk1t*P5*q1k1 - 8*k2(vb)*q1tU*q2tk1t*P5*q2k1 + 4*
     & k2(vb)*q1tU*q2tk1t*P5*q2k2 - 8*k2(vb)*q1tU*q2tk1t*P5*k1k2 + 
     & 4*k2(vb)*q1tU*q2tk1t*P6*q1q2 - 4*k2(vb)*q1tU*q2tk1t*P6*q1k2
     &  + 4*k2(vb)*q1tU*q2tk1t*P8*q1q2 - 8*k2(vb)*q1tU*q2tk1t*P8*
     & q1k1
      AS2(2,vb) = AS2(2,vb) - 4*k2(vb)*q1tU*q2tk2t*P5*q2k1 + 4*
     & k2(vb)*q1tU*q2tk2t*P6*q1k1 - 4*k2(vb)*q1tU*t2*P5*m1*m2 - 4*
     & k2(vb)*q1tU*t2*P5*q1k1 - 4*k2(vb)*q1tU*t2*P5*k1k2 + 4*k2(vb)*
     & q1tU*t2*P6*q1k1 - 4*k2(vb)*q1tU*t2*P8*q1k1 - 4*k2(vb)*q2tU*t1
     & *P5*q2k1 - 4*k2(vb)*q2tU*t1*P6*m1*m2 + 4*k2(vb)*q2tU*t1*P6*
     & q2k1 - 4*k2(vb)*q2tU*t1*P6*k1k2 + 4*k2(vb)*q2tU*t1*P8*q2k1 - 4
     & *kw(vb)*q2U*q1tq2t*P5*m1*m2 - 4*kw(vb)*q2U*q1tq2t*P5*k1k2 + 
     & 4*kw(vb)*q2U*q1tq2t*P6*m1*m2 + 4*kw(vb)*q2U*q1tq2t*P6*k1k2
     &  - 4*kw(vb)*q2U*q1tk1t*q2tk2t*P5 + 4*kw(vb)*q2U*q1tk1t*
     & q2tk2t*P6 + 4*kw(vb)*q2U*q1tk2t*q2tk1t*P5 - 4*kw(vb)*q2U*
     & q1tk2t*q2tk1t*P6 + 4*kw(vb)*q1tq2t*q1tU*P6*m1*m2 + 4*kw(vb)*
     & q1tq2t*q1tU*P6*k1k2 + 4*kw(vb)*q1tq2t*q2tU*P5*m1*m2 + 4*
     & kw(vb)*q1tq2t*q2tU*P5*k1k2 - 4*kw(vb)*q1tk1t*q2tU*P5*q2k2 - 
     & 4*kw(vb)*q1tk1t*q2tU*P6*q1k2 + 4*kw(vb)*q1tk2t*q2tU*P5*q2k1
     &  - 8*kw(vb)*q1tk2t*q2tU*P6*m1*m2 + 4*kw(vb)*q1tk2t*q2tU*P6*
     & q1k1
      AS2(2,vb) = AS2(2,vb) - 8*kw(vb)*q1tk2t*q2tU*P6*k1k2 - 8*
     & kw(vb)*q1tU*q2tk1t*P5*m1*m2 + 4*kw(vb)*q1tU*q2tk1t*P5*q2k2
     &  - 8*kw(vb)*q1tU*q2tk1t*P5*k1k2 + 4*kw(vb)*q1tU*q2tk1t*P6*
     & q1k2 - 4*kw(vb)*q1tU*q2tk2t*P5*q2k1 - 4*kw(vb)*q1tU*q2tk2t*
     & P6*q1k1 - 4*kw(vb)*q1tU*t2*P5*m1*m2 - 4*kw(vb)*q1tU*t2*P5*k1k2
     &  - 4*kw(vb)*q2tU*t1*P6*m1*m2 - 4*kw(vb)*q2tU*t1*P6*k1k2 + 4*
     & q1(vb)*q2U*q1tq2t*P8*m1*m2 + 4*q1(vb)*q2U*q1tq2t*P8*k1k2 + 4
     & *q1(vb)*q2U*q1tk1t*q2tk2t*P8 - 4*q1(vb)*q2U*q1tk2t*q2tk1t*
     & P8 + 4*q1(vb)*q1tq2t*q1tU*P6*m1*m2 + 4*q1(vb)*q1tq2t*q1tU*P6
     & *k1k2 - 4*q1(vb)*q1tq2t*q2tU*P6*m1*m2 - 4*q1(vb)*q1tq2t*q2tU
     & *P6*k1k2 - 4*q1(vb)*q1tq2t*q2tU*P8*m1*m2 - 4*q1(vb)*q1tq2t*
     & q2tU*P8*k1k2 + 4*q1(vb)*q1tk1t*q2tU*P6*q1k2 + 4*q1(vb)*
     & q1tk1t*q2tU*P6*q2k2 + 4*q1(vb)*q1tk1t*q2tU*P8*q2k2 + 8*
     & q1(vb)*q1tk2t*q2tU*P6*m1*m2 - 4*q1(vb)*q1tk2t*q2tU*P6*q1k1
     &  - 4*q1(vb)*q1tk2t*q2tU*P6*q2k1 + 8*q1(vb)*q1tk2t*q2tU*P6*
     & k1k2
      AS2(2,vb) = AS2(2,vb) - 4*q1(vb)*q1tk2t*q2tU*P8*q2k1 - 4*
     & q1(vb)*q1tU*q2tk1t*P6*q1k2 - 4*q1(vb)*q1tU*q2tk1t*P6*q2k2 + 
     & 8*q1(vb)*q1tU*q2tk1t*P8*m1*m2 - 4*q1(vb)*q1tU*q2tk1t*P8*q2k2
     &  + 8*q1(vb)*q1tU*q2tk1t*P8*k1k2 - 8*q1(vb)*q1tU*q2tk2t*P6*m1
     & *m2 + 4*q1(vb)*q1tU*q2tk2t*P6*q1k1 + 4*q1(vb)*q1tU*q2tk2t*P6
     & *q2k1 - 8*q1(vb)*q1tU*q2tk2t*P6*k1k2 + 4*q1(vb)*q1tU*q2tk2t*
     & P8*q2k1 - 4*q1(vb)*q1tU*t2*P6*m1*m2 - 4*q1(vb)*q1tU*t2*P6*k1k2
     &  + 4*q1(vb)*q1tU*t2*P8*m1*m2 + 4*q1(vb)*q1tU*t2*P8*k1k2 + 4*
     & q1(vb)*q2tU*t1*P6*m1*m2 + 4*q1(vb)*q2tU*t1*P6*k1k2 - 4*q2(vb)*
     & q2U*q1tq2t*P8*m1*m2 - 4*q2(vb)*q2U*q1tq2t*P8*k1k2 - 4*q2(vb)
     & *q2U*q1tk1t*q2tk2t*P8 + 4*q2(vb)*q2U*q1tk2t*q2tk1t*P8 - 4*
     & q2(vb)*q1tq2t*q1tU*P5*m1*m2 - 4*q2(vb)*q1tq2t*q1tU*P5*k1k2
     &  - 4*q2(vb)*q1tq2t*q1tU*P8*m1*m2 - 4*q2(vb)*q1tq2t*q1tU*P8*
     & k1k2 + 4*q2(vb)*q1tq2t*q2tU*P5*m1*m2 + 4*q2(vb)*q1tq2t*q2tU*
     & P5*k1k2 - 8*q2(vb)*q1tk1t*q2tU*P5*m1*m2 + 4*q2(vb)*q1tk1t*
     & q2tU*P5*q1k2
      AS2(2,vb) = AS2(2,vb) + 4*q2(vb)*q1tk1t*q2tU*P5*q2k2 - 8*
     & q2(vb)*q1tk1t*q2tU*P5*k1k2 + 4*q2(vb)*q1tk1t*q2tU*P8*q1k2 - 
     & 4*q2(vb)*q1tk2t*q2tU*P5*q1k1 - 4*q2(vb)*q1tk2t*q2tU*P5*q2k1
     &  + 8*q2(vb)*q1tk2t*q2tU*P8*m1*m2 - 4*q2(vb)*q1tk2t*q2tU*P8*
     & q1k1 + 8*q2(vb)*q1tk2t*q2tU*P8*k1k2 + 8*q2(vb)*q1tU*q2tk1t*
     & P5*m1*m2 - 4*q2(vb)*q1tU*q2tk1t*P5*q1k2 - 4*q2(vb)*q1tU*
     & q2tk1t*P5*q2k2 + 8*q2(vb)*q1tU*q2tk1t*P5*k1k2 - 4*q2(vb)*
     & q1tU*q2tk1t*P8*q1k2 + 4*q2(vb)*q1tU*q2tk2t*P5*q1k1 + 4*
     & q2(vb)*q1tU*q2tk2t*P5*q2k1 + 4*q2(vb)*q1tU*q2tk2t*P8*q1k1 + 
     & 4*q2(vb)*q1tU*t2*P5*m1*m2 + 4*q2(vb)*q1tU*t2*P5*k1k2 - 4*
     & q2(vb)*q2tU*t1*P5*m1*m2 - 4*q2(vb)*q2tU*t1*P5*k1k2 + 4*q2(vb)*
     & q2tU*t1*P8*m1*m2 + 4*q2(vb)*q2tU*t1*P8*k1k2 + 4*q1t(vb)*q2U*
     & q1tq2t*P5*m1*m2 + 4*q1t(vb)*q2U*q1tq2t*P5*k1k2 + 4*q1t(vb)*
     & q2U*q1tq2t*P6*m1*m2 + 4*q1t(vb)*q2U*q1tq2t*P6*k1k2 - 4*
     & q1t(vb)*q2U*q1tq2t*P8*m1*m2 - 4*q1t(vb)*q2U*q1tq2t*P8*k1k2
     &  - 8*q1t(vb)*q2U*q2tk1t*P5*m1*m2
      AS2(2,vb) = AS2(2,vb) + 4*q1t(vb)*q2U*q2tk1t*P5*q1k2 + 4*
     & q1t(vb)*q2U*q2tk1t*P5*q2k2 - 8*q1t(vb)*q2U*q2tk1t*P5*k1k2 - 
     & 4*q1t(vb)*q2U*q2tk1t*P6*q1k2 - 4*q1t(vb)*q2U*q2tk1t*P6*q2k2
     &  + 8*q1t(vb)*q2U*q2tk1t*P8*m1*m2 - 4*q1t(vb)*q2U*q2tk1t*P8*
     & q1k2 - 4*q1t(vb)*q2U*q2tk1t*P8*q2k2 + 8*q1t(vb)*q2U*q2tk1t*
     & P8*k1k2 - 4*q1t(vb)*q2U*q2tk2t*P5*q1k1 - 4*q1t(vb)*q2U*
     & q2tk2t*P5*q2k1 - 8*q1t(vb)*q2U*q2tk2t*P6*m1*m2 + 4*q1t(vb)*
     & q2U*q2tk2t*P6*q1k1 + 4*q1t(vb)*q2U*q2tk2t*P6*q2k1 - 8*
     & q1t(vb)*q2U*q2tk2t*P6*k1k2 + 4*q1t(vb)*q2U*q2tk2t*P8*q1k1 + 
     & 4*q1t(vb)*q2U*q2tk2t*P8*q2k1 - 4*q1t(vb)*q2U*t2*P5*m1*m2 - 4*
     & q1t(vb)*q2U*t2*P5*k1k2 - 4*q1t(vb)*q2U*t2*P6*m1*m2 - 4*q1t(vb)
     & *q2U*t2*P6*k1k2 + 4*q1t(vb)*q2U*t2*P8*m1*m2 + 4*q1t(vb)*q2U*
     & t2*P8*k1k2 - 4*q1t(vb)*q2tU*t1*P6*m1*m2 - 4*q1t(vb)*q2tU*t1*P6
     & *k1k2 + 4*q1t(vb)*q2tU*t2*P5*m1*m2 + 4*q1t(vb)*q2tU*t2*P5*k1k2
     &  + 4*q1t(vb)*q2tU*P5*m1*m2*q1q2 - 8*q1t(vb)*q2tU*P5*m1*m2*q2k1
     &  + 4*q1t(vb)*q2tU*P5*q1q2*k1k2
      AS2(2,vb) = AS2(2,vb) - 4*q1t(vb)*q2tU*P5*q1k1*q2k2 + 4*
     & q1t(vb)*q2tU*P5*q1k2*q2k1 - 8*q1t(vb)*q2tU*P5*q2k1*k1k2 - 4*
     & q1t(vb)*q2tU*P6*m1*m2*q1q2 + 8*q1t(vb)*q2tU*P6*m1*m2*q1k2 - 4*
     & q1t(vb)*q2tU*P6*q1q2*k1k2 + 4*q1t(vb)*q2tU*P6*q1k1*q2k2 - 4*
     & q1t(vb)*q2tU*P6*q1k2*q2k1 + 8*q1t(vb)*q2tU*P6*q1k2*k1k2 - 4*
     & q1t(vb)*q2tU*P8*m1*m2*q1q2 - 4*q1t(vb)*q2tU*P8*q1q2*k1k2 + 4*
     & q1t(vb)*q2tU*P8*q1k1*q2k2 - 4*q1t(vb)*q2tU*P8*q1k2*q2k1 - 4*
     & q2t(vb)*q2U*q1tq2t*P5*m1*m2 - 4*q2t(vb)*q2U*q1tq2t*P5*k1k2
     &  - 4*q2t(vb)*q2U*q1tq2t*P6*m1*m2 - 4*q2t(vb)*q2U*q1tq2t*P6*
     & k1k2 + 4*q2t(vb)*q2U*q1tq2t*P8*m1*m2 + 4*q2t(vb)*q2U*q1tq2t*
     & P8*k1k2 + 8*q2t(vb)*q2U*q1tk1t*P5*m1*m2 - 4*q2t(vb)*q2U*
     & q1tk1t*P5*q1k2 - 4*q2t(vb)*q2U*q1tk1t*P5*q2k2 + 8*q2t(vb)*
     & q2U*q1tk1t*P5*k1k2 + 4*q2t(vb)*q2U*q1tk1t*P6*q1k2 + 4*
     & q2t(vb)*q2U*q1tk1t*P6*q2k2 - 4*q2t(vb)*q2U*q1tk1t*P8*q1k2 - 
     & 4*q2t(vb)*q2U*q1tk1t*P8*q2k2 + 4*q2t(vb)*q2U*q1tk2t*P5*q1k1
     &  + 4*q2t(vb)*q2U*q1tk2t*P5*q2k1
      AS2(2,vb) = AS2(2,vb) + 8*q2t(vb)*q2U*q1tk2t*P6*m1*m2 - 4*
     & q2t(vb)*q2U*q1tk2t*P6*q1k1 - 4*q2t(vb)*q2U*q1tk2t*P6*q2k1 + 
     & 8*q2t(vb)*q2U*q1tk2t*P6*k1k2 - 8*q2t(vb)*q2U*q1tk2t*P8*m1*m2
     &  + 4*q2t(vb)*q2U*q1tk2t*P8*q1k1 + 4*q2t(vb)*q2U*q1tk2t*P8*
     & q2k1 - 8*q2t(vb)*q2U*q1tk2t*P8*k1k2 + 4*q2t(vb)*q2U*t1*P5*m1*
     & m2 + 4*q2t(vb)*q2U*t1*P5*k1k2 + 4*q2t(vb)*q2U*t1*P6*m1*m2 + 4*
     & q2t(vb)*q2U*t1*P6*k1k2 - 4*q2t(vb)*q2U*t1*P8*m1*m2 - 4*q2t(vb)
     & *q2U*t1*P8*k1k2 + 4*q2t(vb)*q1tU*t1*P6*m1*m2 + 4*q2t(vb)*q1tU
     & *t1*P6*k1k2 - 4*q2t(vb)*q1tU*t2*P5*m1*m2 - 4*q2t(vb)*q1tU*t2*
     & P5*k1k2 - 4*q2t(vb)*q1tU*P5*m1*m2*q1q2 + 8*q2t(vb)*q1tU*P5*m1*
     & m2*q2k1 - 4*q2t(vb)*q1tU*P5*q1q2*k1k2 + 4*q2t(vb)*q1tU*P5*q1k1
     & *q2k2 - 4*q2t(vb)*q1tU*P5*q1k2*q2k1 + 8*q2t(vb)*q1tU*P5*q2k1*
     & k1k2 + 4*q2t(vb)*q1tU*P6*m1*m2*q1q2 - 8*q2t(vb)*q1tU*P6*m1*m2*
     & q1k2 + 4*q2t(vb)*q1tU*P6*q1q2*k1k2 - 4*q2t(vb)*q1tU*P6*q1k1*
     & q2k2 + 4*q2t(vb)*q1tU*P6*q1k2*q2k1 - 8*q2t(vb)*q1tU*P6*q1k2*
     & k1k2
      AS2(2,vb) = AS2(2,vb) - 4*q2t(vb)*q1tU*P8*m1*m2*q1q2 - 4*
     & q2t(vb)*q1tU*P8*q1q2*k1k2 + 4*q2t(vb)*q1tU*P8*q1k1*q2k2 - 4*
     & q2t(vb)*q1tU*P8*q1k2*q2k1 - 4*U(vb)*t1*t2*P5*m1*m2 - 4*U(vb)*t1
     & *t2*P5*k1k2 - 4*U(vb)*t1*t2*P6*m1*m2 - 4*U(vb)*t1*t2*P6*k1k2 + 4
     & *U(vb)*t1*t2*P8*m1*m2 + 4*U(vb)*t1*t2*P8*k1k2 + 4*U(vb)*q1tq2t*
     & q1tk1t*P5*q2k2 + 8*U(vb)*q1tq2t*q1tk1t*P6*q1k2 + 4*U(vb)*
     & q1tq2t*q1tk1t*P6*q2k2 + 4*U(vb)*q1tq2t*q1tk1t*P8*q2k2 - 4*
     & U(vb)*q1tq2t*q1tk2t*P5*q2k1 + 16*U(vb)*q1tq2t*q1tk2t*P6*m1*
     & m2 - 8*U(vb)*q1tq2t*q1tk2t*P6*q1k1 - 4*U(vb)*q1tq2t*q1tk2t*
     & P6*q2k1 + 16*U(vb)*q1tq2t*q1tk2t*P6*k1k2 - 4*U(vb)*q1tq2t*
     & q1tk2t*P8*q2k1 + 16*U(vb)*q1tq2t*q2tk1t*P5*m1*m2 - 4*U(vb)*
     & q1tq2t*q2tk1t*P5*q1k2 - 8*U(vb)*q1tq2t*q2tk1t*P5*q2k2 + 16*
     & U(vb)*q1tq2t*q2tk1t*P5*k1k2 - 4*U(vb)*q1tq2t*q2tk1t*P6*q1k2
     &  - 4*U(vb)*q1tq2t*q2tk1t*P8*q1k2 + 4*U(vb)*q1tq2t*q2tk2t*P5*
     & q1k1 + 8*U(vb)*q1tq2t*q2tk2t*P5*q2k1 + 4*U(vb)*q1tq2t*q2tk2t
     & *P6*q1k1
      AS2(2,vb) = AS2(2,vb) + 4*U(vb)*q1tq2t*q2tk2t*P8*q1k1 + 4*
     & U(vb)*q1tq2t*t1*P6*m1*m2 + 4*U(vb)*q1tq2t*t1*P6*k1k2 + 4*U(vb)
     & *q1tq2t*t2*P5*m1*m2 + 4*U(vb)*q1tq2t*t2*P5*k1k2 - 4*U(vb)*
     & q1tq2t*P5*m1*m2*q1q2 + 8*U(vb)*q1tq2t*P5*m1*m2*q2k1 - 4*U(vb)*
     & q1tq2t*P5*q1q2*k1k2 + 4*U(vb)*q1tq2t*P5*q1k1*q2k2 - 4*U(vb)*
     & q1tq2t*P5*q1k2*q2k1 + 8*U(vb)*q1tq2t*P5*q2k1*k1k2 - 4*U(vb)*
     & q1tq2t*P6*m1*m2*q1q2 + 8*U(vb)*q1tq2t*P6*m1*m2*q1k2 - 4*U(vb)*
     & q1tq2t*P6*q1q2*k1k2 + 4*U(vb)*q1tq2t*P6*q1k1*q2k2 - 4*U(vb)*
     & q1tq2t*P6*q1k2*q2k1 + 8*U(vb)*q1tq2t*P6*q1k2*k1k2 - 4*U(vb)*
     & q1tq2t*P8*m1*m2*q1q2 - 4*U(vb)*q1tq2t*P8*q1q2*k1k2 + 4*U(vb)*
     & q1tq2t*P8*q1k1*q2k2 - 4*U(vb)*q1tq2t*P8*q1k2*q2k1 - 4*U(vb)*
     & q1tq2t**2*P5*m1*m2 - 4*U(vb)*q1tq2t**2*P5*k1k2 - 4*U(vb)*
     & q1tq2t**2*P6*m1*m2 - 4*U(vb)*q1tq2t**2*P6*k1k2 - 4*U(vb)*
     & q1tq2t**2*P8*m1*m2 - 4*U(vb)*q1tq2t**2*P8*k1k2 - 16*U(vb)*
     & q1tk1t*q2tk1t*P5*m1*m2 + 8*U(vb)*q1tk1t*q2tk1t*P5*q1k2 + 8*
     & U(vb)*q1tk1t*q2tk1t*P5*q2k2
      AS2(2,vb) = AS2(2,vb) - 16*U(vb)*q1tk1t*q2tk1t*P5*k1k2 + 8*
     & U(vb)*q1tk1t*q2tk1t*P8*q1k2 - 4*U(vb)*q1tk1t*q2tk2t*t1*P6 - 
     & 4*U(vb)*q1tk1t*q2tk2t*t2*P5 - 4*U(vb)*q1tk1t*q2tk2t*P5*q1q2
     &  - 4*U(vb)*q1tk1t*q2tk2t*P6*q1q2 - 4*U(vb)*q1tk1t*q2tk2t*P8*
     & q1q2 - 8*U(vb)*q1tk1t*t2*P5*m1*m2 + 4*U(vb)*q1tk1t*t2*P5*q1k2
     &  - 8*U(vb)*q1tk1t*t2*P5*k1k2 - 4*U(vb)*q1tk1t*t2*P6*q1k2 + 4*
     & U(vb)*q1tk1t*t2*P8*q1k2 + 4*U(vb)*q1tk2t*q2tk1t*t1*P6 + 4*
     & U(vb)*q1tk2t*q2tk1t*t2*P5 + 4*U(vb)*q1tk2t*q2tk1t*P5*q1q2 - 
     & 8*U(vb)*q1tk2t*q2tk1t*P5*q1k1 - 8*U(vb)*q1tk2t*q2tk1t*P5*
     & q2k1 + 4*U(vb)*q1tk2t*q2tk1t*P6*q1q2 - 8*U(vb)*q1tk2t*q2tk1t
     & *P6*q1k2 - 8*U(vb)*q1tk2t*q2tk1t*P6*q2k2 + 16*U(vb)*q1tk2t*
     & q2tk1t*P8*m1*m2 + 4*U(vb)*q1tk2t*q2tk1t*P8*q1q2 - 8*U(vb)*
     & q1tk2t*q2tk1t*P8*q1k1 - 8*U(vb)*q1tk2t*q2tk1t*P8*q2k2 + 16*
     & U(vb)*q1tk2t*q2tk1t*P8*k1k2 - 16*U(vb)*q1tk2t*q2tk2t*P6*m1*
     & m2 + 8*U(vb)*q1tk2t*q2tk2t*P6*q1k1 + 8*U(vb)*q1tk2t*q2tk2t*
     & P6*q2k1
      AS2(2,vb) = AS2(2,vb) - 16*U(vb)*q1tk2t*q2tk2t*P6*k1k2 + 8*
     & U(vb)*q1tk2t*q2tk2t*P8*q2k1 - 4*U(vb)*q1tk2t*t2*P5*q1k1 - 8*
     & U(vb)*q1tk2t*t2*P6*m1*m2 + 4*U(vb)*q1tk2t*t2*P6*q1k1 - 8*U(vb)
     & *q1tk2t*t2*P6*k1k2 + 8*U(vb)*q1tk2t*t2*P8*m1*m2 - 4*U(vb)*
     & q1tk2t*t2*P8*q1k1 + 8*U(vb)*q1tk2t*t2*P8*k1k2 - 8*U(vb)*
     & q2tk1t*t1*P5*m1*m2 + 4*U(vb)*q2tk1t*t1*P5*q2k2 - 8*U(vb)*
     & q2tk1t*t1*P5*k1k2 - 4*U(vb)*q2tk1t*t1*P6*q2k2 + 8*U(vb)*
     & q2tk1t*t1*P8*m1*m2 - 4*U(vb)*q2tk1t*t1*P8*q2k2 + 8*U(vb)*
     & q2tk1t*t1*P8*k1k2 - 4*U(vb)*q2tk2t*t1*P5*q2k1 - 8*U(vb)*
     & q2tk2t*t1*P6*m1*m2 + 4*U(vb)*q2tk2t*t1*P6*q2k1 - 8*U(vb)*
     & q2tk2t*t1*P6*k1k2 + 4*U(vb)*q2tk2t*t1*P8*q2k1

      AS2(3,vb) =
     & 8*k1(vb)*q1L*q1tq2t*q1tk2t*P6 - 4*k1(vb)*q1L*q1tq2t*q2tk2t
     & *P5 - 4*k1(vb)*q1L*q1tq2t*q2tk2t*P6 - 4*k1(vb)*q1L*q1tq2t*
     & q2tk2t*P8 - 4*k1(vb)*q1L*q1tq2t*P5*q2k2 - 4*k1(vb)*q1L*
     & q1tq2t*P6*m1*m2 - 4*k1(vb)*q1L*q1tq2t*P6*q2k2 - 4*k1(vb)*q1L
     & *q1tq2t*P6*k1k2 - 4*k1(vb)*q1L*q1tq2t*P8*q2k2 - 4*k1(vb)*q1L
     & *q1tk1t*q2tk2t*P6 + 8*k1(vb)*q1L*q1tk2t*q2tk1t*P5 + 4*
     & k1(vb)*q1L*q1tk2t*q2tk1t*P6 + 8*k1(vb)*q1L*q1tk2t*q2tk1t*
     & P8 - 8*k1(vb)*q1L*q1tk2t*q2tk2t*P6 + 4*k1(vb)*q1L*q1tk2t*t2
     & *P5 - 4*k1(vb)*q1L*q1tk2t*t2*P6 + 4*k1(vb)*q1L*q1tk2t*t2*P8
     &  + 4*k1(vb)*q2L*q1tq2t*q1tk2t*P5 + 4*k1(vb)*q2L*q1tq2t*
     & q1tk2t*P6 + 4*k1(vb)*q2L*q1tq2t*q1tk2t*P8 - 8*k1(vb)*q2L*
     & q1tq2t*q2tk2t*P5 - 4*k1(vb)*q2L*q1tq2t*P5*m1*m2 + 4*k1(vb)*
     & q2L*q1tq2t*P5*q1k2 - 4*k1(vb)*q2L*q1tq2t*P5*k1k2 + 4*k1(vb)*
     & q2L*q1tq2t*P6*q1k2 + 4*k1(vb)*q2L*q1tq2t*P8*q1k2 + 4*k1(vb)*
     & q2L*q1tk1t*q2tk2t*P5
      AS2(3,vb) = AS2(3,vb) + 4*k1(vb)*q2L*q1tk2t*q2tk1t*P5 - 8*
     & k1(vb)*q2L*q1tk2t*q2tk2t*P6 - 8*k1(vb)*q2L*q1tk2t*q2tk2t*
     & P8 + 4*k1(vb)*q2L*q2tk2t*t1*P5 - 4*k1(vb)*q2L*q2tk2t*t1*P6
     &  - 4*k1(vb)*q2L*q2tk2t*t1*P8 + 4*k1(vb)*q1tq2t*q1tL*P5*q2k2
     &  + 4*k1(vb)*q1tq2t*q1tL*P6*m1*m2 + 8*k1(vb)*q1tq2t*q1tL*P6*
     & q1k2 + 4*k1(vb)*q1tq2t*q1tL*P6*q2k2 + 4*k1(vb)*q1tq2t*q1tL*
     & P6*k1k2 + 4*k1(vb)*q1tq2t*q1tL*P8*q2k2 + 4*k1(vb)*q1tq2t*
     & q2tL*P5*m1*m2 - 4*k1(vb)*q1tq2t*q2tL*P5*q1k2 - 8*k1(vb)*
     & q1tq2t*q2tL*P5*q2k2 + 4*k1(vb)*q1tq2t*q2tL*P5*k1k2 - 4*
     & k1(vb)*q1tq2t*q2tL*P6*q1k2 - 4*k1(vb)*q1tq2t*q2tL*P8*q1k2 + 
     & 4*k1(vb)*q1tk1t*q2tL*P5*q2k2 - 4*k1(vb)*q1tk1t*q2tL*P6*q1k2
     &  + 4*k1(vb)*q1tk2t*q2tL*t1*P6 + 4*k1(vb)*q1tk2t*q2tL*t2*P5
     &  + 4*k1(vb)*q1tk2t*q2tL*P5*q1q2 - 4*k1(vb)*q1tk2t*q2tL*P5*
     & q2k1 - 8*k1(vb)*q1tk2t*q2tL*P6*m1*m2 + 4*k1(vb)*q1tk2t*q2tL*
     & P6*q1q2 + 4*k1(vb)*q1tk2t*q2tL*P6*q1k1 - 8*k1(vb)*q1tk2t*
     & q2tL*P6*q1k2
      AS2(3,vb) = AS2(3,vb) - 8*k1(vb)*q1tk2t*q2tL*P6*q2k2 - 8*
     & k1(vb)*q1tk2t*q2tL*P6*k1k2 + 4*k1(vb)*q1tk2t*q2tL*P8*q1q2 - 
     & 8*k1(vb)*q1tk2t*q2tL*P8*q2k2 - 8*k1(vb)*q1tL*q2tk1t*P5*m1*m2
     &  + 8*k1(vb)*q1tL*q2tk1t*P5*q1k2 + 4*k1(vb)*q1tL*q2tk1t*P5*
     & q2k2 - 8*k1(vb)*q1tL*q2tk1t*P5*k1k2 + 4*k1(vb)*q1tL*q2tk1t*
     & P6*q1k2 + 8*k1(vb)*q1tL*q2tk1t*P8*q1k2 - 4*k1(vb)*q1tL*
     & q2tk2t*t1*P6 - 4*k1(vb)*q1tL*q2tk2t*t2*P5 - 4*k1(vb)*q1tL*
     & q2tk2t*P5*q1q2 + 4*k1(vb)*q1tL*q2tk2t*P5*q2k1 - 4*k1(vb)*
     & q1tL*q2tk2t*P6*q1q2 - 4*k1(vb)*q1tL*q2tk2t*P6*q1k1 - 4*
     & k1(vb)*q1tL*q2tk2t*P8*q1q2 - 4*k1(vb)*q1tL*t2*P5*m1*m2 + 4*
     & k1(vb)*q1tL*t2*P5*q1k2 - 4*k1(vb)*q1tL*t2*P5*k1k2 - 4*k1(vb)*
     & q1tL*t2*P6*q1k2 + 4*k1(vb)*q1tL*t2*P8*q1k2 + 4*k1(vb)*q2tL*t1
     & *P5*q2k2 - 4*k1(vb)*q2tL*t1*P6*m1*m2 - 4*k1(vb)*q2tL*t1*P6*
     & q2k2 - 4*k1(vb)*q2tL*t1*P6*k1k2 - 4*k1(vb)*q2tL*t1*P8*q2k2 - 8
     & *k2(vb)*q1L*q1tq2t*q1tk1t*P6 + 4*k2(vb)*q1L*q1tq2t*q2tk1t*
     & P5
      AS2(3,vb) = AS2(3,vb) + 4*k2(vb)*q1L*q1tq2t*q2tk1t*P6 + 4*
     & k2(vb)*q1L*q1tq2t*q2tk1t*P8 + 4*k2(vb)*q1L*q1tq2t*P5*q2k1
     &  - 4*k2(vb)*q1L*q1tq2t*P6*m1*m2 + 4*k2(vb)*q1L*q1tq2t*P6*
     & q2k1 - 4*k2(vb)*q1L*q1tq2t*P6*k1k2 + 4*k2(vb)*q1L*q1tq2t*P8*
     & q2k1 - 8*k2(vb)*q1L*q1tk1t*q2tk1t*P5 - 8*k2(vb)*q1L*q1tk1t*
     & q2tk1t*P8 + 4*k2(vb)*q1L*q1tk1t*q2tk2t*P6 - 4*k2(vb)*q1L*
     & q1tk1t*t2*P5 + 4*k2(vb)*q1L*q1tk1t*t2*P6 - 4*k2(vb)*q1L*
     & q1tk1t*t2*P8 + 4*k2(vb)*q1L*q1tk2t*q2tk1t*P6 - 4*k2(vb)*q2L
     & *q1tq2t*q1tk1t*P5 - 4*k2(vb)*q2L*q1tq2t*q1tk1t*P6 - 4*
     & k2(vb)*q2L*q1tq2t*q1tk1t*P8 + 8*k2(vb)*q2L*q1tq2t*q2tk1t*
     & P5 - 4*k2(vb)*q2L*q1tq2t*P5*m1*m2 - 4*k2(vb)*q2L*q1tq2t*P5*
     & q1k1 - 4*k2(vb)*q2L*q1tq2t*P5*k1k2 - 4*k2(vb)*q2L*q1tq2t*P6*
     & q1k1 - 4*k2(vb)*q2L*q1tq2t*P8*q1k1 - 8*k2(vb)*q2L*q1tk1t*
     & q2tk1t*P5 - 4*k2(vb)*q2L*q1tk1t*q2tk2t*P5 + 4*k2(vb)*q2L*
     & q1tk2t*q2tk1t*P5 + 8*k2(vb)*q2L*q1tk2t*q2tk1t*P6 + 8*k2(vb)
     & *q2L*q1tk2t*q2tk1t*P8
      AS2(3,vb) = AS2(3,vb) - 4*k2(vb)*q2L*q2tk1t*t1*P5 + 4*k2(vb)*
     & q2L*q2tk1t*t1*P6 + 4*k2(vb)*q2L*q2tk1t*t1*P8 - 4*k2(vb)*
     & q1tq2t*q1tL*P5*q2k1 + 4*k2(vb)*q1tq2t*q1tL*P6*m1*m2 - 8*
     & k2(vb)*q1tq2t*q1tL*P6*q1k1 - 4*k2(vb)*q1tq2t*q1tL*P6*q2k1 + 
     & 4*k2(vb)*q1tq2t*q1tL*P6*k1k2 - 4*k2(vb)*q1tq2t*q1tL*P8*q2k1
     &  + 4*k2(vb)*q1tq2t*q2tL*P5*m1*m2 + 4*k2(vb)*q1tq2t*q2tL*P5*
     & q1k1 + 8*k2(vb)*q1tq2t*q2tL*P5*q2k1 + 4*k2(vb)*q1tq2t*q2tL*
     & P5*k1k2 + 4*k2(vb)*q1tq2t*q2tL*P6*q1k1 + 4*k2(vb)*q1tq2t*
     & q2tL*P8*q1k1 - 4*k2(vb)*q1tk1t*q2tL*t1*P6 - 4*k2(vb)*q1tk1t*
     & q2tL*t2*P5 - 4*k2(vb)*q1tk1t*q2tL*P5*q1q2 - 4*k2(vb)*q1tk1t*
     & q2tL*P5*q2k2 - 4*k2(vb)*q1tk1t*q2tL*P6*q1q2 + 4*k2(vb)*
     & q1tk1t*q2tL*P6*q1k2 - 4*k2(vb)*q1tk1t*q2tL*P8*q1q2 + 4*
     & k2(vb)*q1tk2t*q2tL*P5*q2k1 - 8*k2(vb)*q1tk2t*q2tL*P6*m1*m2
     &  + 4*k2(vb)*q1tk2t*q2tL*P6*q1k1 + 8*k2(vb)*q1tk2t*q2tL*P6*
     & q2k1 - 8*k2(vb)*q1tk2t*q2tL*P6*k1k2 + 8*k2(vb)*q1tk2t*q2tL*
     & P8*q2k1
      AS2(3,vb) = AS2(3,vb) + 4*k2(vb)*q1tL*q2tk1t*t1*P6 + 4*k2(vb)
     & *q1tL*q2tk1t*t2*P5 - 8*k2(vb)*q1tL*q2tk1t*P5*m1*m2 + 4*
     & k2(vb)*q1tL*q2tk1t*P5*q1q2 - 8*k2(vb)*q1tL*q2tk1t*P5*q1k1 - 
     & 8*k2(vb)*q1tL*q2tk1t*P5*q2k1 + 4*k2(vb)*q1tL*q2tk1t*P5*q2k2
     &  - 8*k2(vb)*q1tL*q2tk1t*P5*k1k2 + 4*k2(vb)*q1tL*q2tk1t*P6*
     & q1q2 - 4*k2(vb)*q1tL*q2tk1t*P6*q1k2 + 4*k2(vb)*q1tL*q2tk1t*
     & P8*q1q2 - 8*k2(vb)*q1tL*q2tk1t*P8*q1k1 - 4*k2(vb)*q1tL*
     & q2tk2t*P5*q2k1 + 4*k2(vb)*q1tL*q2tk2t*P6*q1k1 - 4*k2(vb)*
     & q1tL*t2*P5*m1*m2 - 4*k2(vb)*q1tL*t2*P5*q1k1 - 4*k2(vb)*q1tL*
     & t2*P5*k1k2 + 4*k2(vb)*q1tL*t2*P6*q1k1 - 4*k2(vb)*q1tL*t2*P8*
     & q1k1 - 4*k2(vb)*q2tL*t1*P5*q2k1 - 4*k2(vb)*q2tL*t1*P6*m1*m2 + 
     & 4*k2(vb)*q2tL*t1*P6*q2k1 - 4*k2(vb)*q2tL*t1*P6*k1k2 + 4*k2(vb)
     & *q2tL*t1*P8*q2k1 - 4*kw(vb)*q1L*q1tq2t*P6*m1*m2 - 4*kw(vb)*
     & q1L*q1tq2t*P6*k1k2 - 4*kw(vb)*q1L*q1tk1t*q2tk2t*P6 + 4*
     & kw(vb)*q1L*q1tk2t*q2tk1t*P6 - 4*kw(vb)*q2L*q1tq2t*P5*m1*m2
     &  - 4*kw(vb)*q2L*q1tq2t*P5*k1k2
      AS2(3,vb) = AS2(3,vb) - 4*kw(vb)*q2L*q1tk1t*q2tk2t*P5 + 4*
     & kw(vb)*q2L*q1tk2t*q2tk1t*P5 + 4*kw(vb)*q1tq2t*q1tL*P6*m1*m2
     &  + 4*kw(vb)*q1tq2t*q1tL*P6*k1k2 + 4*kw(vb)*q1tq2t*q2tL*P5*m1
     & *m2 + 4*kw(vb)*q1tq2t*q2tL*P5*k1k2 - 4*kw(vb)*q1tk1t*q2tL*P5
     & *q2k2 - 4*kw(vb)*q1tk1t*q2tL*P6*q1k2 + 4*kw(vb)*q1tk2t*q2tL*
     & P5*q2k1 - 8*kw(vb)*q1tk2t*q2tL*P6*m1*m2 + 4*kw(vb)*q1tk2t*
     & q2tL*P6*q1k1 - 8*kw(vb)*q1tk2t*q2tL*P6*k1k2 - 8*kw(vb)*q1tL*
     & q2tk1t*P5*m1*m2 + 4*kw(vb)*q1tL*q2tk1t*P5*q2k2 - 8*kw(vb)*
     & q1tL*q2tk1t*P5*k1k2 + 4*kw(vb)*q1tL*q2tk1t*P6*q1k2 - 4*
     & kw(vb)*q1tL*q2tk2t*P5*q2k1 - 4*kw(vb)*q1tL*q2tk2t*P6*q1k1 - 
     & 4*kw(vb)*q1tL*t2*P5*m1*m2 - 4*kw(vb)*q1tL*t2*P5*k1k2 - 4*
     & kw(vb)*q2tL*t1*P6*m1*m2 - 4*kw(vb)*q2tL*t1*P6*k1k2 + 4*q1(vb)*
     & q1L*q1tq2t*P6*m1*m2 + 4*q1(vb)*q1L*q1tq2t*P6*k1k2 + 4*q1(vb)
     & *q1L*q1tk1t*q2tk2t*P6 - 4*q1(vb)*q1L*q1tk2t*q2tk1t*P6 + 4*
     & q1(vb)*q2L*q1tq2t*P6*m1*m2 + 4*q1(vb)*q2L*q1tq2t*P6*k1k2 + 4
     & *q1(vb)*q2L*q1tq2t*P8*m1*m2
      AS2(3,vb) = AS2(3,vb) + 4*q1(vb)*q2L*q1tq2t*P8*k1k2 + 4*
     & q1(vb)*q2L*q1tk1t*q2tk2t*P6 + 4*q1(vb)*q2L*q1tk1t*q2tk2t*
     & P8 - 4*q1(vb)*q2L*q1tk2t*q2tk1t*P6 - 4*q1(vb)*q2L*q1tk2t*
     & q2tk1t*P8 + 4*q1(vb)*q1tq2t*q1tL*P6*m1*m2 + 4*q1(vb)*q1tq2t*
     & q1tL*P6*k1k2 - 4*q1(vb)*q1tq2t*q2tL*P6*m1*m2 - 4*q1(vb)*
     & q1tq2t*q2tL*P6*k1k2 - 4*q1(vb)*q1tq2t*q2tL*P8*m1*m2 - 4*
     & q1(vb)*q1tq2t*q2tL*P8*k1k2 + 4*q1(vb)*q1tk1t*q2tL*P6*q1k2 + 
     & 4*q1(vb)*q1tk1t*q2tL*P6*q2k2 + 4*q1(vb)*q1tk1t*q2tL*P8*q2k2
     &  + 8*q1(vb)*q1tk2t*q2tL*P6*m1*m2 - 4*q1(vb)*q1tk2t*q2tL*P6*
     & q1k1 - 4*q1(vb)*q1tk2t*q2tL*P6*q2k1 + 8*q1(vb)*q1tk2t*q2tL*
     & P6*k1k2 - 4*q1(vb)*q1tk2t*q2tL*P8*q2k1 - 4*q1(vb)*q1tL*
     & q2tk1t*P6*q1k2 - 4*q1(vb)*q1tL*q2tk1t*P6*q2k2 + 8*q1(vb)*
     & q1tL*q2tk1t*P8*m1*m2 - 4*q1(vb)*q1tL*q2tk1t*P8*q2k2 + 8*
     & q1(vb)*q1tL*q2tk1t*P8*k1k2 - 8*q1(vb)*q1tL*q2tk2t*P6*m1*m2
     &  + 4*q1(vb)*q1tL*q2tk2t*P6*q1k1 + 4*q1(vb)*q1tL*q2tk2t*P6*
     & q2k1
      AS2(3,vb) = AS2(3,vb) - 8*q1(vb)*q1tL*q2tk2t*P6*k1k2 + 4*
     & q1(vb)*q1tL*q2tk2t*P8*q2k1 - 4*q1(vb)*q1tL*t2*P6*m1*m2 - 4*
     & q1(vb)*q1tL*t2*P6*k1k2 + 4*q1(vb)*q1tL*t2*P8*m1*m2 + 4*q1(vb)*
     & q1tL*t2*P8*k1k2 + 4*q1(vb)*q2tL*t1*P6*m1*m2 + 4*q1(vb)*q2tL*
     & t1*P6*k1k2 + 4*q2(vb)*q1L*q1tq2t*P5*m1*m2 + 4*q2(vb)*q1L*
     & q1tq2t*P5*k1k2 + 4*q2(vb)*q1L*q1tq2t*P8*m1*m2 + 4*q2(vb)*q1L
     & *q1tq2t*P8*k1k2 + 4*q2(vb)*q1L*q1tk1t*q2tk2t*P5 + 4*q2(vb)*
     & q1L*q1tk1t*q2tk2t*P8 - 4*q2(vb)*q1L*q1tk2t*q2tk1t*P5 - 4*
     & q2(vb)*q1L*q1tk2t*q2tk1t*P8 + 4*q2(vb)*q2L*q1tq2t*P5*m1*m2
     &  + 4*q2(vb)*q2L*q1tq2t*P5*k1k2 + 4*q2(vb)*q2L*q1tk1t*q2tk2t
     & *P5 - 4*q2(vb)*q2L*q1tk2t*q2tk1t*P5 - 4*q2(vb)*q1tq2t*q1tL*
     & P5*m1*m2 - 4*q2(vb)*q1tq2t*q1tL*P5*k1k2 - 4*q2(vb)*q1tq2t*
     & q1tL*P8*m1*m2 - 4*q2(vb)*q1tq2t*q1tL*P8*k1k2 + 4*q2(vb)*
     & q1tq2t*q2tL*P5*m1*m2 + 4*q2(vb)*q1tq2t*q2tL*P5*k1k2 - 8*
     & q2(vb)*q1tk1t*q2tL*P5*m1*m2 + 4*q2(vb)*q1tk1t*q2tL*P5*q1k2
     &  + 4*q2(vb)*q1tk1t*q2tL*P5*q2k2
      AS2(3,vb) = AS2(3,vb) - 8*q2(vb)*q1tk1t*q2tL*P5*k1k2 + 4*
     & q2(vb)*q1tk1t*q2tL*P8*q1k2 - 4*q2(vb)*q1tk2t*q2tL*P5*q1k1 - 
     & 4*q2(vb)*q1tk2t*q2tL*P5*q2k1 + 8*q2(vb)*q1tk2t*q2tL*P8*m1*m2
     &  - 4*q2(vb)*q1tk2t*q2tL*P8*q1k1 + 8*q2(vb)*q1tk2t*q2tL*P8*
     & k1k2 + 8*q2(vb)*q1tL*q2tk1t*P5*m1*m2 - 4*q2(vb)*q1tL*q2tk1t*
     & P5*q1k2 - 4*q2(vb)*q1tL*q2tk1t*P5*q2k2 + 8*q2(vb)*q1tL*
     & q2tk1t*P5*k1k2 - 4*q2(vb)*q1tL*q2tk1t*P8*q1k2 + 4*q2(vb)*
     & q1tL*q2tk2t*P5*q1k1 + 4*q2(vb)*q1tL*q2tk2t*P5*q2k1 + 4*
     & q2(vb)*q1tL*q2tk2t*P8*q1k1 + 4*q2(vb)*q1tL*t2*P5*m1*m2 + 4*
     & q2(vb)*q1tL*t2*P5*k1k2 - 4*q2(vb)*q2tL*t1*P5*m1*m2 - 4*q2(vb)*
     & q2tL*t1*P5*k1k2 + 4*q2(vb)*q2tL*t1*P8*m1*m2 + 4*q2(vb)*q2tL*
     & t1*P8*k1k2 - 8*q1t(vb)*q1L*q1tq2t*P6*m1*m2 - 8*q1t(vb)*q1L*
     & q1tq2t*P6*k1k2 + 8*q1t(vb)*q1L*q2tk1t*P5*m1*m2 - 4*q1t(vb)*
     & q1L*q2tk1t*P5*q2k2 + 8*q1t(vb)*q1L*q2tk1t*P5*k1k2 + 4*
     & q1t(vb)*q1L*q2tk1t*P6*q2k2 - 8*q1t(vb)*q1L*q2tk1t*P8*m1*m2
     &  + 4*q1t(vb)*q1L*q2tk1t*P8*q2k2
      AS2(3,vb) = AS2(3,vb) - 8*q1t(vb)*q1L*q2tk1t*P8*k1k2 + 4*
     & q1t(vb)*q1L*q2tk2t*P5*q2k1 + 8*q1t(vb)*q1L*q2tk2t*P6*m1*m2
     &  - 4*q1t(vb)*q1L*q2tk2t*P6*q2k1 + 8*q1t(vb)*q1L*q2tk2t*P6*
     & k1k2 - 4*q1t(vb)*q1L*q2tk2t*P8*q2k1 + 4*q1t(vb)*q1L*t2*P5*m1*
     & m2 + 4*q1t(vb)*q1L*t2*P5*k1k2 + 4*q1t(vb)*q1L*t2*P6*m1*m2 + 4*
     & q1t(vb)*q1L*t2*P6*k1k2 - 4*q1t(vb)*q1L*t2*P8*m1*m2 - 4*q1t(vb)
     & *q1L*t2*P8*k1k2 + 4*q1t(vb)*q2L*q1tq2t*P5*m1*m2 + 4*q1t(vb)*
     & q2L*q1tq2t*P5*k1k2 - 4*q1t(vb)*q2L*q1tq2t*P6*m1*m2 - 4*
     & q1t(vb)*q2L*q1tq2t*P6*k1k2 - 4*q1t(vb)*q2L*q1tq2t*P8*m1*m2
     &  - 4*q1t(vb)*q2L*q1tq2t*P8*k1k2 + 4*q1t(vb)*q2L*q2tk1t*P5*
     & q1k2 - 4*q1t(vb)*q2L*q2tk1t*P6*q1k2 - 4*q1t(vb)*q2L*q2tk1t*
     & P8*q1k2 - 4*q1t(vb)*q2L*q2tk2t*P5*q1k1 + 4*q1t(vb)*q2L*
     & q2tk2t*P6*q1k1 + 4*q1t(vb)*q2L*q2tk2t*P8*q1k1 - 4*q1t(vb)*
     & q2tL*t1*P6*m1*m2 - 4*q1t(vb)*q2tL*t1*P6*k1k2 + 4*q1t(vb)*q2tL
     & *t2*P5*m1*m2 + 4*q1t(vb)*q2tL*t2*P5*k1k2 + 4*q1t(vb)*q2tL*P5*
     & m1*m2*q1q2
      AS2(3,vb) = AS2(3,vb) - 8*q1t(vb)*q2tL*P5*m1*m2*q2k1 + 4*
     & q1t(vb)*q2tL*P5*q1q2*k1k2 - 4*q1t(vb)*q2tL*P5*q1k1*q2k2 + 4*
     & q1t(vb)*q2tL*P5*q1k2*q2k1 - 8*q1t(vb)*q2tL*P5*q2k1*k1k2 - 4*
     & q1t(vb)*q2tL*P6*m1*m2*q1q2 + 8*q1t(vb)*q2tL*P6*m1*m2*q1k2 - 4*
     & q1t(vb)*q2tL*P6*q1q2*k1k2 + 4*q1t(vb)*q2tL*P6*q1k1*q2k2 - 4*
     & q1t(vb)*q2tL*P6*q1k2*q2k1 + 8*q1t(vb)*q2tL*P6*q1k2*k1k2 - 4*
     & q1t(vb)*q2tL*P8*m1*m2*q1q2 - 4*q1t(vb)*q2tL*P8*q1q2*k1k2 + 4*
     & q1t(vb)*q2tL*P8*q1k1*q2k2 - 4*q1t(vb)*q2tL*P8*q1k2*q2k1 - 4*
     & q2t(vb)*q1L*q1tq2t*P5*m1*m2 - 4*q2t(vb)*q1L*q1tq2t*P5*k1k2
     &  + 4*q2t(vb)*q1L*q1tq2t*P6*m1*m2 + 4*q2t(vb)*q1L*q1tq2t*P6*
     & k1k2 - 4*q2t(vb)*q1L*q1tq2t*P8*m1*m2 - 4*q2t(vb)*q1L*q1tq2t*
     & P8*k1k2 + 4*q2t(vb)*q1L*q1tk1t*P5*q2k2 - 4*q2t(vb)*q1L*
     & q1tk1t*P6*q2k2 + 4*q2t(vb)*q1L*q1tk1t*P8*q2k2 - 4*q2t(vb)*
     & q1L*q1tk2t*P5*q2k1 + 4*q2t(vb)*q1L*q1tk2t*P6*q2k1 - 4*
     & q2t(vb)*q1L*q1tk2t*P8*q2k1 - 8*q2t(vb)*q2L*q1tq2t*P5*m1*m2
     &  - 8*q2t(vb)*q2L*q1tq2t*P5*k1k2
      AS2(3,vb) = AS2(3,vb) + 8*q2t(vb)*q2L*q1tk1t*P5*m1*m2 - 4*
     & q2t(vb)*q2L*q1tk1t*P5*q1k2 + 8*q2t(vb)*q2L*q1tk1t*P5*k1k2 + 
     & 4*q2t(vb)*q2L*q1tk1t*P6*q1k2 - 4*q2t(vb)*q2L*q1tk1t*P8*q1k2
     &  + 4*q2t(vb)*q2L*q1tk2t*P5*q1k1 + 8*q2t(vb)*q2L*q1tk2t*P6*m1
     & *m2 - 4*q2t(vb)*q2L*q1tk2t*P6*q1k1 + 8*q2t(vb)*q2L*q1tk2t*P6
     & *k1k2 - 8*q2t(vb)*q2L*q1tk2t*P8*m1*m2 + 4*q2t(vb)*q2L*q1tk2t
     & *P8*q1k1 - 8*q2t(vb)*q2L*q1tk2t*P8*k1k2 + 4*q2t(vb)*q2L*t1*P5
     & *m1*m2 + 4*q2t(vb)*q2L*t1*P5*k1k2 + 4*q2t(vb)*q2L*t1*P6*m1*m2
     &  + 4*q2t(vb)*q2L*t1*P6*k1k2 - 4*q2t(vb)*q2L*t1*P8*m1*m2 - 4*
     & q2t(vb)*q2L*t1*P8*k1k2 + 4*q2t(vb)*q1tL*t1*P6*m1*m2 + 4*
     & q2t(vb)*q1tL*t1*P6*k1k2 - 4*q2t(vb)*q1tL*t2*P5*m1*m2 - 4*
     & q2t(vb)*q1tL*t2*P5*k1k2 - 4*q2t(vb)*q1tL*P5*m1*m2*q1q2 + 8*
     & q2t(vb)*q1tL*P5*m1*m2*q2k1 - 4*q2t(vb)*q1tL*P5*q1q2*k1k2 + 4*
     & q2t(vb)*q1tL*P5*q1k1*q2k2 - 4*q2t(vb)*q1tL*P5*q1k2*q2k1 + 8*
     & q2t(vb)*q1tL*P5*q2k1*k1k2 + 4*q2t(vb)*q1tL*P6*m1*m2*q1q2 - 8*
     & q2t(vb)*q1tL*P6*m1*m2*q1k2
      AS2(3,vb) = AS2(3,vb) + 4*q2t(vb)*q1tL*P6*q1q2*k1k2 - 4*
     & q2t(vb)*q1tL*P6*q1k1*q2k2 + 4*q2t(vb)*q1tL*P6*q1k2*q2k1 - 8*
     & q2t(vb)*q1tL*P6*q1k2*k1k2 - 4*q2t(vb)*q1tL*P8*m1*m2*q1q2 - 4*
     & q2t(vb)*q1tL*P8*q1q2*k1k2 + 4*q2t(vb)*q1tL*P8*q1k1*q2k2 - 4*
     & q2t(vb)*q1tL*P8*q1k2*q2k1 - 4*L(vb)*t1*t2*P5*m1*m2 - 4*L(vb)*t1
     & *t2*P5*k1k2 - 4*L(vb)*t1*t2*P6*m1*m2 - 4*L(vb)*t1*t2*P6*k1k2 + 4
     & *L(vb)*t1*t2*P8*m1*m2 + 4*L(vb)*t1*t2*P8*k1k2 + 4*L(vb)*q1tq2t*
     & q1tk1t*P5*q2k2 + 8*L(vb)*q1tq2t*q1tk1t*P6*q1k2 + 4*L(vb)*
     & q1tq2t*q1tk1t*P6*q2k2 + 4*L(vb)*q1tq2t*q1tk1t*P8*q2k2 - 4*
     & L(vb)*q1tq2t*q1tk2t*P5*q2k1 + 16*L(vb)*q1tq2t*q1tk2t*P6*m1*
     & m2 - 8*L(vb)*q1tq2t*q1tk2t*P6*q1k1 - 4*L(vb)*q1tq2t*q1tk2t*
     & P6*q2k1 + 16*L(vb)*q1tq2t*q1tk2t*P6*k1k2 - 4*L(vb)*q1tq2t*
     & q1tk2t*P8*q2k1 + 16*L(vb)*q1tq2t*q2tk1t*P5*m1*m2 - 4*L(vb)*
     & q1tq2t*q2tk1t*P5*q1k2 - 8*L(vb)*q1tq2t*q2tk1t*P5*q2k2 + 16*
     & L(vb)*q1tq2t*q2tk1t*P5*k1k2 - 4*L(vb)*q1tq2t*q2tk1t*P6*q1k2
     &  - 4*L(vb)*q1tq2t*q2tk1t*P8*q1k2
      AS2(3,vb) = AS2(3,vb) + 4*L(vb)*q1tq2t*q2tk2t*P5*q1k1 + 8*
     & L(vb)*q1tq2t*q2tk2t*P5*q2k1 + 4*L(vb)*q1tq2t*q2tk2t*P6*q1k1
     &  + 4*L(vb)*q1tq2t*q2tk2t*P8*q1k1 + 4*L(vb)*q1tq2t*t1*P6*m1*m2
     &  + 4*L(vb)*q1tq2t*t1*P6*k1k2 + 4*L(vb)*q1tq2t*t2*P5*m1*m2 + 4*
     & L(vb)*q1tq2t*t2*P5*k1k2 - 4*L(vb)*q1tq2t*P5*m1*m2*q1q2 + 8*
     & L(vb)*q1tq2t*P5*m1*m2*q2k1 - 4*L(vb)*q1tq2t*P5*q1q2*k1k2 + 4*
     & L(vb)*q1tq2t*P5*q1k1*q2k2 - 4*L(vb)*q1tq2t*P5*q1k2*q2k1 + 8*
     & L(vb)*q1tq2t*P5*q2k1*k1k2 - 4*L(vb)*q1tq2t*P6*m1*m2*q1q2 + 8*
     & L(vb)*q1tq2t*P6*m1*m2*q1k2 - 4*L(vb)*q1tq2t*P6*q1q2*k1k2 + 4*
     & L(vb)*q1tq2t*P6*q1k1*q2k2 - 4*L(vb)*q1tq2t*P6*q1k2*q2k1 + 8*
     & L(vb)*q1tq2t*P6*q1k2*k1k2 - 4*L(vb)*q1tq2t*P8*m1*m2*q1q2 - 4*
     & L(vb)*q1tq2t*P8*q1q2*k1k2 + 4*L(vb)*q1tq2t*P8*q1k1*q2k2 - 4*
     & L(vb)*q1tq2t*P8*q1k2*q2k1 - 4*L(vb)*q1tq2t**2*P5*m1*m2 - 4*
     & L(vb)*q1tq2t**2*P5*k1k2 - 4*L(vb)*q1tq2t**2*P6*m1*m2 - 4*L(vb)
     & *q1tq2t**2*P6*k1k2 - 4*L(vb)*q1tq2t**2*P8*m1*m2 - 4*L(vb)*
     & q1tq2t**2*P8*k1k2
      AS2(3,vb) = AS2(3,vb) - 16*L(vb)*q1tk1t*q2tk1t*P5*m1*m2 + 8*
     & L(vb)*q1tk1t*q2tk1t*P5*q1k2 + 8*L(vb)*q1tk1t*q2tk1t*P5*q2k2
     &  - 16*L(vb)*q1tk1t*q2tk1t*P5*k1k2 + 8*L(vb)*q1tk1t*q2tk1t*P8
     & *q1k2 - 4*L(vb)*q1tk1t*q2tk2t*t1*P6 - 4*L(vb)*q1tk1t*q2tk2t*
     & t2*P5 - 4*L(vb)*q1tk1t*q2tk2t*P5*q1q2 - 4*L(vb)*q1tk1t*
     & q2tk2t*P6*q1q2 - 4*L(vb)*q1tk1t*q2tk2t*P8*q1q2 - 8*L(vb)*
     & q1tk1t*t2*P5*m1*m2 + 4*L(vb)*q1tk1t*t2*P5*q1k2 - 8*L(vb)*
     & q1tk1t*t2*P5*k1k2 - 4*L(vb)*q1tk1t*t2*P6*q1k2 + 4*L(vb)*
     & q1tk1t*t2*P8*q1k2 + 4*L(vb)*q1tk2t*q2tk1t*t1*P6 + 4*L(vb)*
     & q1tk2t*q2tk1t*t2*P5 + 4*L(vb)*q1tk2t*q2tk1t*P5*q1q2 - 8*
     & L(vb)*q1tk2t*q2tk1t*P5*q1k1 - 8*L(vb)*q1tk2t*q2tk1t*P5*q2k1
     &  + 4*L(vb)*q1tk2t*q2tk1t*P6*q1q2 - 8*L(vb)*q1tk2t*q2tk1t*P6*
     & q1k2 - 8*L(vb)*q1tk2t*q2tk1t*P6*q2k2 + 16*L(vb)*q1tk2t*
     & q2tk1t*P8*m1*m2 + 4*L(vb)*q1tk2t*q2tk1t*P8*q1q2 - 8*L(vb)*
     & q1tk2t*q2tk1t*P8*q1k1 - 8*L(vb)*q1tk2t*q2tk1t*P8*q2k2 + 16*
     & L(vb)*q1tk2t*q2tk1t*P8*k1k2
      AS2(3,vb) = AS2(3,vb) - 16*L(vb)*q1tk2t*q2tk2t*P6*m1*m2 + 8*
     & L(vb)*q1tk2t*q2tk2t*P6*q1k1 + 8*L(vb)*q1tk2t*q2tk2t*P6*q2k1
     &  - 16*L(vb)*q1tk2t*q2tk2t*P6*k1k2 + 8*L(vb)*q1tk2t*q2tk2t*P8
     & *q2k1 - 4*L(vb)*q1tk2t*t2*P5*q1k1 - 8*L(vb)*q1tk2t*t2*P6*m1*m2
     &  + 4*L(vb)*q1tk2t*t2*P6*q1k1 - 8*L(vb)*q1tk2t*t2*P6*k1k2 + 8*
     & L(vb)*q1tk2t*t2*P8*m1*m2 - 4*L(vb)*q1tk2t*t2*P8*q1k1 + 8*L(vb)
     & *q1tk2t*t2*P8*k1k2 - 8*L(vb)*q2tk1t*t1*P5*m1*m2 + 4*L(vb)*
     & q2tk1t*t1*P5*q2k2 - 8*L(vb)*q2tk1t*t1*P5*k1k2 - 4*L(vb)*
     & q2tk1t*t1*P6*q2k2 + 8*L(vb)*q2tk1t*t1*P8*m1*m2 - 4*L(vb)*
     & q2tk1t*t1*P8*q2k2 + 8*L(vb)*q2tk1t*t1*P8*k1k2 - 4*L(vb)*
     & q2tk2t*t1*P5*q2k1 - 8*L(vb)*q2tk2t*t1*P6*m1*m2 + 4*L(vb)*
     & q2tk2t*t1*P6*q2k1 - 8*L(vb)*q2tk2t*t1*P6*k1k2 + 4*L(vb)*
     & q2tk2t*t1*P8*q2k1

      AS2(4,vb) =
     & 4*k1(vb)*q1L*q2U*q1tq2t*P5*m2 + 4*k1(vb)*q1L*q2U*q1tq2t*P6
     & *m2 + 4*k1(vb)*q1L*q2U*q1tq2t*P8*m2 + 8*k1(vb)*q1L*q1tq2t*
     & q1tU*P6*m2 - 4*k1(vb)*q1L*q1tq2t*q2tU*P5*m2 - 4*k1(vb)*q1L*
     & q1tq2t*q2tU*P6*m2 - 4*k1(vb)*q1L*q1tq2t*q2tU*P8*m2 - 4*
     & k1(vb)*q1L*q1tk1t*q2tU*P6*m2 + 4*k1(vb)*q1L*q1tk2t*q2tU*P6
     & *m1 + 8*k1(vb)*q1L*q1tU*q2tk1t*P5*m2 + 4*k1(vb)*q1L*q1tU*
     & q2tk1t*P6*m2 + 8*k1(vb)*q1L*q1tU*q2tk1t*P8*m2 - 4*k1(vb)*
     & q1L*q1tU*q2tk2t*P6*m1 - 8*k1(vb)*q1L*q1tU*q2tk2t*P6*m2 + 4
     & *k1(vb)*q1L*q1tU*t2*P5*m2 - 4*k1(vb)*q1L*q1tU*t2*P6*m2 + 4*
     & k1(vb)*q1L*q1tU*t2*P8*m2 + 4*k1(vb)*q2U*q2L*q1tq2t*P5*m2 + 
     & 4*k1(vb)*q2U*q2L*q1tq2t*P6*m2 + 4*k1(vb)*q2U*q2L*q1tq2t*P8
     & *m2 - 4*k1(vb)*q2U*q1tq2t*q1tL*P5*m2 + 4*k1(vb)*q2U*q1tq2t*
     & q1tL*P6*m2 - 4*k1(vb)*q2U*q1tq2t*q1tL*P8*m2 + 4*k1(vb)*q2U*
     & q1tq2t*q2tL*P5*m2 - 4*k1(vb)*q2U*q1tq2t*q2tL*P6*m2 - 4*
     & k1(vb)*q2U*q1tq2t*q2tL*P8*m2
      AS2(4,vb) = AS2(4,vb) - 4*k1(vb)*q2U*q1tk1t*q2tL*P5*m2 - 4*
     & k1(vb)*q2U*q1tk1t*q2tL*P6*m2 - 4*k1(vb)*q2U*q1tk2t*q2tL*P5
     & *m1 + 4*k1(vb)*q2U*q1tk2t*q2tL*P6*m1 + 8*k1(vb)*q2U*q1tk2t*
     & q2tL*P6*m2 + 8*k1(vb)*q2U*q1tk2t*q2tL*P8*m2 + 4*k1(vb)*q2U*
     & q1tL*q2tk1t*P5*m2 + 4*k1(vb)*q2U*q1tL*q2tk1t*P6*m2 + 8*
     & k1(vb)*q2U*q1tL*q2tk1t*P8*m2 + 4*k1(vb)*q2U*q1tL*q2tk2t*P5
     & *m1 - 4*k1(vb)*q2U*q1tL*q2tk2t*P6*m1 - 8*k1(vb)*q2U*q1tL*
     & q2tk2t*P6*m2 + 4*k1(vb)*q2U*q1tL*t2*P5*m2 - 4*k1(vb)*q2U*
     & q1tL*t2*P6*m2 + 4*k1(vb)*q2U*q1tL*t2*P8*m2 - 4*k1(vb)*q2U*
     & q2tL*t1*P5*m2 + 4*k1(vb)*q2U*q2tL*t1*P6*m2 + 4*k1(vb)*q2U*
     & q2tL*t1*P8*m2 + 4*k1(vb)*q2L*q1tq2t*q1tU*P5*m2 + 4*k1(vb)*
     & q2L*q1tq2t*q1tU*P6*m2 + 4*k1(vb)*q2L*q1tq2t*q1tU*P8*m2 - 8
     & *k1(vb)*q2L*q1tq2t*q2tU*P5*m2 + 4*k1(vb)*q2L*q1tk1t*q2tU*
     & P5*m2 + 4*k1(vb)*q2L*q1tk2t*q2tU*P5*m1 - 8*k1(vb)*q2L*
     & q1tk2t*q2tU*P6*m2 - 8*k1(vb)*q2L*q1tk2t*q2tU*P8*m2 + 4*
     & k1(vb)*q2L*q1tU*q2tk1t*P5*m2
      AS2(4,vb) = AS2(4,vb) - 4*k1(vb)*q2L*q1tU*q2tk2t*P5*m1 + 4*
     & k1(vb)*q2L*q2tU*t1*P5*m2 - 4*k1(vb)*q2L*q2tU*t1*P6*m2 - 4*
     & k1(vb)*q2L*q2tU*t1*P8*m2 + 4*k1(vb)*q1tU*q2tL*t1*P6*m2 + 4*
     & k1(vb)*q1tU*q2tL*t2*P5*m2 - 4*k1(vb)*q1tU*q2tL*P5*m1*q2k2 +
     & 4*k1(vb)*q1tU*q2tL*P5*m2*q1q2 - 4*k1(vb)*q1tU*q2tL*P5*m2*
     & q2k1 - 4*k1(vb)*q1tU*q2tL*P6*m1*q1k2 + 4*k1(vb)*q1tU*q2tL*P6
     & *m2*q1q2 + 4*k1(vb)*q1tU*q2tL*P6*m2*q1k1 - 8*k1(vb)*q1tU*
     & q2tL*P6*m2*q1k2 + 4*k1(vb)*q1tU*q2tL*P8*m2*q1q2 - 4*k1(vb)*
     & q1tL*q2tU*t1*P6*m2 - 4*k1(vb)*q1tL*q2tU*t2*P5*m2 + 4*k1(vb)*
     & q1tL*q2tU*P5*m1*q2k2 - 4*k1(vb)*q1tL*q2tU*P5*m2*q1q2 + 4*
     & k1(vb)*q1tL*q2tU*P5*m2*q2k1 + 4*k1(vb)*q1tL*q2tU*P6*m1*q1k2
     &  - 4*k1(vb)*q1tL*q2tU*P6*m2*q1q2 - 4*k1(vb)*q1tL*q2tU*P6*m2*
     & q1k1 + 8*k1(vb)*q1tL*q2tU*P6*m2*q1k2 - 4*k1(vb)*q1tL*q2tU*P8
     & *m2*q1q2 - 4*k2(vb)*q1L*q2U*q1tq2t*P5*m1 - 4*k2(vb)*q1L*q2U
     & *q1tq2t*P6*m1 - 4*k2(vb)*q1L*q2U*q1tq2t*P8*m1 - 8*k2(vb)*
     & q1L*q1tq2t*q1tU*P6*m1
      AS2(4,vb) = AS2(4,vb) + 4*k2(vb)*q1L*q1tq2t*q2tU*P5*m1 + 4*
     & k2(vb)*q1L*q1tq2t*q2tU*P6*m1 + 4*k2(vb)*q1L*q1tq2t*q2tU*P8
     & *m1 - 4*k2(vb)*q1L*q1tk1t*q2tU*P6*m2 + 4*k2(vb)*q1L*q1tk2t*
     & q2tU*P6*m1 - 8*k2(vb)*q1L*q1tU*q2tk1t*P5*m1 + 4*k2(vb)*q1L*
     & q1tU*q2tk1t*P6*m2 - 8*k2(vb)*q1L*q1tU*q2tk1t*P8*m1 + 4*
     & k2(vb)*q1L*q1tU*q2tk2t*P6*m1 - 4*k2(vb)*q1L*q1tU*t2*P5*m1
     &  + 4*k2(vb)*q1L*q1tU*t2*P6*m1 - 4*k2(vb)*q1L*q1tU*t2*P8*m1
     &  - 4*k2(vb)*q2U*q2L*q1tq2t*P5*m1 - 4*k2(vb)*q2U*q2L*q1tq2t
     & *P6*m1 - 4*k2(vb)*q2U*q2L*q1tq2t*P8*m1 + 4*k2(vb)*q2U*
     & q1tq2t*q1tL*P5*m1 - 4*k2(vb)*q2U*q1tq2t*q1tL*P6*m1 + 4*
     & k2(vb)*q2U*q1tq2t*q1tL*P8*m1 - 4*k2(vb)*q2U*q1tq2t*q2tL*P5
     & *m1 + 4*k2(vb)*q2U*q1tq2t*q2tL*P6*m1 + 4*k2(vb)*q2U*q1tq2t*
     & q2tL*P8*m1 + 8*k2(vb)*q2U*q1tk1t*q2tL*P5*m1 + 4*k2(vb)*q2U*
     & q1tk1t*q2tL*P5*m2 - 4*k2(vb)*q2U*q1tk1t*q2tL*P6*m2 - 4*
     & k2(vb)*q2U*q1tk2t*q2tL*P5*m1 - 4*k2(vb)*q2U*q1tk2t*q2tL*P6
     & *m1
      AS2(4,vb) = AS2(4,vb) - 8*k2(vb)*q2U*q1tk2t*q2tL*P8*m1 - 8*
     & k2(vb)*q2U*q1tL*q2tk1t*P5*m1 - 4*k2(vb)*q2U*q1tL*q2tk1t*P5
     & *m2 + 4*k2(vb)*q2U*q1tL*q2tk1t*P6*m2 - 8*k2(vb)*q2U*q1tL*
     & q2tk1t*P8*m1 + 4*k2(vb)*q2U*q1tL*q2tk2t*P5*m1 + 4*k2(vb)*
     & q2U*q1tL*q2tk2t*P6*m1 - 4*k2(vb)*q2U*q1tL*t2*P5*m1 + 4*
     & k2(vb)*q2U*q1tL*t2*P6*m1 - 4*k2(vb)*q2U*q1tL*t2*P8*m1 + 4*
     & k2(vb)*q2U*q2tL*t1*P5*m1 - 4*k2(vb)*q2U*q2tL*t1*P6*m1 - 4*
     & k2(vb)*q2U*q2tL*t1*P8*m1 - 4*k2(vb)*q2L*q1tq2t*q1tU*P5*m1
     &  - 4*k2(vb)*q2L*q1tq2t*q1tU*P6*m1 - 4*k2(vb)*q2L*q1tq2t*
     & q1tU*P8*m1 + 8*k2(vb)*q2L*q1tq2t*q2tU*P5*m1 - 8*k2(vb)*q2L*
     & q1tk1t*q2tU*P5*m1 - 4*k2(vb)*q2L*q1tk1t*q2tU*P5*m2 + 4*
     & k2(vb)*q2L*q1tk2t*q2tU*P5*m1 + 8*k2(vb)*q2L*q1tk2t*q2tU*P6
     & *m1 + 8*k2(vb)*q2L*q1tk2t*q2tU*P8*m1 + 4*k2(vb)*q2L*q1tU*
     & q2tk1t*P5*m2 - 4*k2(vb)*q2L*q1tU*q2tk2t*P5*m1 - 4*k2(vb)*
     & q2L*q2tU*t1*P5*m1 + 4*k2(vb)*q2L*q2tU*t1*P6*m1 + 4*k2(vb)*
     & q2L*q2tU*t1*P8*m1
      AS2(4,vb) = AS2(4,vb) - 4*k2(vb)*q1tU*q2tL*t1*P6*m1 - 4*
     & k2(vb)*q1tU*q2tL*t2*P5*m1 - 4*k2(vb)*q1tU*q2tL*P5*m1*q1q2 +
     & 8*k2(vb)*q1tU*q2tL*P5*m1*q2k1 - 4*k2(vb)*q1tU*q2tL*P5*m1*
     & q2k2 + 4*k2(vb)*q1tU*q2tL*P5*m2*q2k1 - 4*k2(vb)*q1tU*q2tL*P6
     & *m1*q1q2 + 4*k2(vb)*q1tU*q2tL*P6*m1*q1k2 + 4*k2(vb)*q1tU*
     & q2tL*P6*m2*q1k1 - 4*k2(vb)*q1tU*q2tL*P8*m1*q1q2 + 4*k2(vb)*
     & q1tL*q2tU*t1*P6*m1 + 4*k2(vb)*q1tL*q2tU*t2*P5*m1 + 4*k2(vb)*
     & q1tL*q2tU*P5*m1*q1q2 - 8*k2(vb)*q1tL*q2tU*P5*m1*q2k1 + 4*
     & k2(vb)*q1tL*q2tU*P5*m1*q2k2 - 4*k2(vb)*q1tL*q2tU*P5*m2*q2k1
     &  + 4*k2(vb)*q1tL*q2tU*P6*m1*q1q2 - 4*k2(vb)*q1tL*q2tU*P6*m1*
     & q1k2 - 4*k2(vb)*q1tL*q2tU*P6*m2*q1k1 + 4*k2(vb)*q1tL*q2tU*P8
     & *m1*q1q2 - 4*kw(vb)*q1L*q1tk1t*q2tU*P6*m2 + 4*kw(vb)*q1L*
     & q1tk2t*q2tU*P6*m1 + 4*kw(vb)*q1L*q1tU*q2tk1t*P6*m2 - 4*
     & kw(vb)*q1L*q1tU*q2tk2t*P6*m1 + 4*kw(vb)*q2U*q1tk1t*q2tL*P5
     & *m2 - 4*kw(vb)*q2U*q1tk1t*q2tL*P6*m2 - 4*kw(vb)*q2U*q1tk2t*
     & q2tL*P5*m1
      AS2(4,vb) = AS2(4,vb) + 4*kw(vb)*q2U*q1tk2t*q2tL*P6*m1 - 4*
     & kw(vb)*q2U*q1tL*q2tk1t*P5*m2 + 4*kw(vb)*q2U*q1tL*q2tk1t*P6
     & *m2 + 4*kw(vb)*q2U*q1tL*q2tk2t*P5*m1 - 4*kw(vb)*q2U*q1tL*
     & q2tk2t*P6*m1 - 4*kw(vb)*q2L*q1tk1t*q2tU*P5*m2 + 4*kw(vb)*
     & q2L*q1tk2t*q2tU*P5*m1 + 4*kw(vb)*q2L*q1tU*q2tk1t*P5*m2 - 4
     & *kw(vb)*q2L*q1tU*q2tk2t*P5*m1 - 4*kw(vb)*q1tU*q2tL*P5*m1*
     & q2k2 + 4*kw(vb)*q1tU*q2tL*P5*m2*q2k1 - 4*kw(vb)*q1tU*q2tL*P6
     & *m1*q1k2 + 4*kw(vb)*q1tU*q2tL*P6*m2*q1k1 + 4*kw(vb)*q1tL*
     & q2tU*P5*m1*q2k2 - 4*kw(vb)*q1tL*q2tU*P5*m2*q2k1 + 4*kw(vb)*
     & q1tL*q2tU*P6*m1*q1k2 - 4*kw(vb)*q1tL*q2tU*P6*m2*q1k1 + 4*
     & q1(vb)*q1L*q1tk1t*q2tU*P6*m2 - 4*q1(vb)*q1L*q1tk2t*q2tU*P6
     & *m1 - 4*q1(vb)*q1L*q1tU*q2tk1t*P6*m2 + 4*q1(vb)*q1L*q1tU*
     & q2tk2t*P6*m1 - 4*q1(vb)*q2U*q1tk1t*q2tL*P8*m2 + 4*q1(vb)*
     & q2U*q1tk2t*q2tL*P8*m1 + 4*q1(vb)*q2U*q1tL*q2tk1t*P8*m2 - 4
     & *q1(vb)*q2U*q1tL*q2tk2t*P8*m1 + 4*q1(vb)*q2L*q1tk1t*q2tU*
     & P6*m2
      AS2(4,vb) = AS2(4,vb) + 4*q1(vb)*q2L*q1tk1t*q2tU*P8*m2 - 4*
     & q1(vb)*q2L*q1tk2t*q2tU*P6*m1 - 4*q1(vb)*q2L*q1tk2t*q2tU*P8
     & *m1 - 4*q1(vb)*q2L*q1tU*q2tk1t*P6*m2 - 4*q1(vb)*q2L*q1tU*
     & q2tk1t*P8*m2 + 4*q1(vb)*q2L*q1tU*q2tk2t*P6*m1 + 4*q1(vb)*
     & q2L*q1tU*q2tk2t*P8*m1 + 4*q1(vb)*q1tU*q2tL*P6*m1*q1k2 + 4*
     & q1(vb)*q1tU*q2tL*P6*m1*q2k2 - 4*q1(vb)*q1tU*q2tL*P6*m2*q1k1
     &  - 4*q1(vb)*q1tU*q2tL*P6*m2*q2k1 + 4*q1(vb)*q1tU*q2tL*P8*m1*
     & q2k2 - 4*q1(vb)*q1tU*q2tL*P8*m2*q2k1 - 4*q1(vb)*q1tL*q2tU*P6
     & *m1*q1k2 - 4*q1(vb)*q1tL*q2tU*P6*m1*q2k2 + 4*q1(vb)*q1tL*
     & q2tU*P6*m2*q1k1 + 4*q1(vb)*q1tL*q2tU*P6*m2*q2k1 - 4*q1(vb)*
     & q1tL*q2tU*P8*m1*q2k2 + 4*q1(vb)*q1tL*q2tU*P8*m2*q2k1 + 4*
     & q2(vb)*q1L*q1tk1t*q2tU*P5*m2 + 4*q2(vb)*q1L*q1tk1t*q2tU*P8
     & *m2 - 4*q2(vb)*q1L*q1tk2t*q2tU*P5*m1 - 4*q2(vb)*q1L*q1tk2t*
     & q2tU*P8*m1 - 4*q2(vb)*q1L*q1tU*q2tk1t*P5*m2 - 4*q2(vb)*q1L*
     & q1tU*q2tk1t*P8*m2 + 4*q2(vb)*q1L*q1tU*q2tk2t*P5*m1 + 4*
     & q2(vb)*q1L*q1tU*q2tk2t*P8*m1
      AS2(4,vb) = AS2(4,vb) + 4*q2(vb)*q2U*q1tk1t*q2tL*P8*m2 - 4*
     & q2(vb)*q2U*q1tk2t*q2tL*P8*m1 - 4*q2(vb)*q2U*q1tL*q2tk1t*P8
     & *m2 + 4*q2(vb)*q2U*q1tL*q2tk2t*P8*m1 + 4*q2(vb)*q2L*q1tk1t*
     & q2tU*P5*m2 - 4*q2(vb)*q2L*q1tk2t*q2tU*P5*m1 - 4*q2(vb)*q2L*
     & q1tU*q2tk1t*P5*m2 + 4*q2(vb)*q2L*q1tU*q2tk2t*P5*m1 + 4*
     & q2(vb)*q1tU*q2tL*P5*m1*q1k2 + 4*q2(vb)*q1tU*q2tL*P5*m1*q2k2
     &  - 4*q2(vb)*q1tU*q2tL*P5*m2*q1k1 - 4*q2(vb)*q1tU*q2tL*P5*m2*
     & q2k1 + 4*q2(vb)*q1tU*q2tL*P8*m1*q1k2 - 4*q2(vb)*q1tU*q2tL*P8
     & *m2*q1k1 - 4*q2(vb)*q1tL*q2tU*P5*m1*q1k2 - 4*q2(vb)*q1tL*
     & q2tU*P5*m1*q2k2 + 4*q2(vb)*q1tL*q2tU*P5*m2*q1k1 + 4*q2(vb)*
     & q1tL*q2tU*P5*m2*q2k1 - 4*q2(vb)*q1tL*q2tU*P8*m1*q1k2 + 4*
     & q2(vb)*q1tL*q2tU*P8*m2*q1k1 + 4*q1t(vb)*q1L*q2U*q2tk1t*P5*
     & m2 - 4*q1t(vb)*q1L*q2U*q2tk1t*P6*m2 - 4*q1t(vb)*q1L*q2U*
     & q2tk1t*P8*m2 - 4*q1t(vb)*q1L*q2U*q2tk2t*P5*m1 + 4*q1t(vb)*
     & q1L*q2U*q2tk2t*P6*m1 + 4*q1t(vb)*q1L*q2U*q2tk2t*P8*m1 - 4*
     & q1t(vb)*q1L*q2tU*P5*m1*q2k2
      AS2(4,vb) = AS2(4,vb) + 4*q1t(vb)*q1L*q2tU*P5*m2*q2k1 + 4*
     & q1t(vb)*q1L*q2tU*P6*m1*q2k2 - 4*q1t(vb)*q1L*q2tU*P6*m2*q2k1
     &  + 4*q1t(vb)*q1L*q2tU*P8*m1*q2k2 - 4*q1t(vb)*q1L*q2tU*P8*m2*
     & q2k1 + 4*q1t(vb)*q2U*q2L*q2tk1t*P5*m2 - 4*q1t(vb)*q2U*q2L*
     & q2tk1t*P6*m2 - 4*q1t(vb)*q2U*q2L*q2tk1t*P8*m2 - 4*q1t(vb)*
     & q2U*q2L*q2tk2t*P5*m1 + 4*q1t(vb)*q2U*q2L*q2tk2t*P6*m1 + 4*
     & q1t(vb)*q2U*q2L*q2tk2t*P8*m1 - 4*q1t(vb)*q2U*q2tL*P5*m1*
     & q1k2 - 4*q1t(vb)*q2U*q2tL*P5*m1*q2k2 + 4*q1t(vb)*q2U*q2tL*P5
     & *m2*q1k1 + 4*q1t(vb)*q2U*q2tL*P5*m2*q2k1 + 4*q1t(vb)*q2U*
     & q2tL*P6*m1*q1k2 + 4*q1t(vb)*q2U*q2tL*P6*m1*q2k2 - 4*q1t(vb)*
     & q2U*q2tL*P6*m2*q1k1 - 4*q1t(vb)*q2U*q2tL*P6*m2*q2k1 + 4*
     & q1t(vb)*q2U*q2tL*P8*m1*q1k2 + 4*q1t(vb)*q2U*q2tL*P8*m1*q2k2
     &  - 4*q1t(vb)*q2U*q2tL*P8*m2*q1k1 - 4*q1t(vb)*q2U*q2tL*P8*m2*
     & q2k1 + 4*q1t(vb)*q2L*q2tU*P5*m1*q1k2 - 4*q1t(vb)*q2L*q2tU*P5
     & *m2*q1k1 - 4*q1t(vb)*q2L*q2tU*P6*m1*q1k2 + 4*q1t(vb)*q2L*
     & q2tU*P6*m2*q1k1
      AS2(4,vb) = AS2(4,vb) - 4*q1t(vb)*q2L*q2tU*P8*m1*q1k2 + 4*
     & q1t(vb)*q2L*q2tU*P8*m2*q1k1 - 4*q2t(vb)*q1L*q2U*q1tk1t*P5*
     & m2 + 4*q2t(vb)*q1L*q2U*q1tk1t*P6*m2 - 4*q2t(vb)*q1L*q2U*
     & q1tk1t*P8*m2 + 4*q2t(vb)*q1L*q2U*q1tk2t*P5*m1 - 4*q2t(vb)*
     & q1L*q2U*q1tk2t*P6*m1 + 4*q2t(vb)*q1L*q2U*q1tk2t*P8*m1 + 4*
     & q2t(vb)*q1L*q1tU*P5*m1*q2k2 - 4*q2t(vb)*q1L*q1tU*P5*m2*q2k1
     &  - 4*q2t(vb)*q1L*q1tU*P6*m1*q2k2 + 4*q2t(vb)*q1L*q1tU*P6*m2*
     & q2k1 + 4*q2t(vb)*q1L*q1tU*P8*m1*q2k2 - 4*q2t(vb)*q1L*q1tU*P8
     & *m2*q2k1 - 4*q2t(vb)*q2U*q2L*q1tk1t*P5*m2 + 4*q2t(vb)*q2U*
     & q2L*q1tk1t*P6*m2 - 4*q2t(vb)*q2U*q2L*q1tk1t*P8*m2 + 4*
     & q2t(vb)*q2U*q2L*q1tk2t*P5*m1 - 4*q2t(vb)*q2U*q2L*q1tk2t*P6
     & *m1 + 4*q2t(vb)*q2U*q2L*q1tk2t*P8*m1 + 4*q2t(vb)*q2U*q1tL*
     & P5*m1*q1k2 + 4*q2t(vb)*q2U*q1tL*P5*m1*q2k2 - 4*q2t(vb)*q2U*
     & q1tL*P5*m2*q1k1 - 4*q2t(vb)*q2U*q1tL*P5*m2*q2k1 - 4*q2t(vb)*
     & q2U*q1tL*P6*m1*q1k2 - 4*q2t(vb)*q2U*q1tL*P6*m1*q2k2 + 4*
     & q2t(vb)*q2U*q1tL*P6*m2*q1k1
      AS2(4,vb) = AS2(4,vb) + 4*q2t(vb)*q2U*q1tL*P6*m2*q2k1 + 4*
     & q2t(vb)*q2U*q1tL*P8*m1*q1k2 + 4*q2t(vb)*q2U*q1tL*P8*m1*q2k2
     &  - 4*q2t(vb)*q2U*q1tL*P8*m2*q1k1 - 4*q2t(vb)*q2U*q1tL*P8*m2*
     & q2k1 - 4*q2t(vb)*q2L*q1tU*P5*m1*q1k2 + 4*q2t(vb)*q2L*q1tU*P5
     & *m2*q1k1 + 4*q2t(vb)*q2L*q1tU*P6*m1*q1k2 - 4*q2t(vb)*q2L*
     & q1tU*P6*m2*q1k1 - 4*q2t(vb)*q2L*q1tU*P8*m1*q1k2 + 4*q2t(vb)*
     & q2L*q1tU*P8*m2*q1k1 + 8*U(vb)*q1L*q1tq2t*q1tk1t*P6*m2 - 8*
     & U(vb)*q1L*q1tq2t*q1tk2t*P6*m1 - 4*U(vb)*q1L*q1tq2t*q2tk1t*
     & P5*m2 - 4*U(vb)*q1L*q1tq2t*q2tk1t*P6*m2 - 4*U(vb)*q1L*
     & q1tq2t*q2tk1t*P8*m2 + 4*U(vb)*q1L*q1tq2t*q2tk2t*P5*m1 + 4*
     & U(vb)*q1L*q1tq2t*q2tk2t*P6*m1 + 4*U(vb)*q1L*q1tq2t*q2tk2t*
     & P8*m1 + 4*U(vb)*q1L*q1tq2t*P5*m1*q2k2 - 4*U(vb)*q1L*q1tq2t*
     & P5*m2*q2k1 + 4*U(vb)*q1L*q1tq2t*P6*m1*q2k2 - 4*U(vb)*q1L*
     & q1tq2t*P6*m2*q2k1 + 4*U(vb)*q1L*q1tq2t*P8*m1*q2k2 - 4*U(vb)*
     & q1L*q1tq2t*P8*m2*q2k1 + 8*U(vb)*q1L*q1tk1t*q2tk1t*P5*m2 + 8
     & *U(vb)*q1L*q1tk1t*q2tk1t*P8*m2
      AS2(4,vb) = AS2(4,vb) - 8*U(vb)*q1L*q1tk1t*q2tk2t*P6*m2 + 4*
     & U(vb)*q1L*q1tk1t*t2*P5*m2 - 4*U(vb)*q1L*q1tk1t*t2*P6*m2 + 4*
     & U(vb)*q1L*q1tk1t*t2*P8*m2 - 8*U(vb)*q1L*q1tk2t*q2tk1t*P5*m1
     &  - 8*U(vb)*q1L*q1tk2t*q2tk1t*P8*m1 + 8*U(vb)*q1L*q1tk2t*
     & q2tk2t*P6*m1 - 4*U(vb)*q1L*q1tk2t*t2*P5*m1 + 4*U(vb)*q1L*
     & q1tk2t*t2*P6*m1 - 4*U(vb)*q1L*q1tk2t*t2*P8*m1 + 4*U(vb)*q2L*
     & q1tq2t*q1tk1t*P5*m2 + 4*U(vb)*q2L*q1tq2t*q1tk1t*P6*m2 + 4*
     & U(vb)*q2L*q1tq2t*q1tk1t*P8*m2 - 4*U(vb)*q2L*q1tq2t*q1tk2t*
     & P5*m1 - 4*U(vb)*q2L*q1tq2t*q1tk2t*P6*m1 - 4*U(vb)*q2L*
     & q1tq2t*q1tk2t*P8*m1 - 8*U(vb)*q2L*q1tq2t*q2tk1t*P5*m2 + 8*
     & U(vb)*q2L*q1tq2t*q2tk2t*P5*m1 - 4*U(vb)*q2L*q1tq2t*P5*m1*
     & q1k2 + 4*U(vb)*q2L*q1tq2t*P5*m2*q1k1 - 4*U(vb)*q2L*q1tq2t*P6
     & *m1*q1k2 + 4*U(vb)*q2L*q1tq2t*P6*m2*q1k1 - 4*U(vb)*q2L*
     & q1tq2t*P8*m1*q1k2 + 4*U(vb)*q2L*q1tq2t*P8*m2*q1k1 + 8*U(vb)*
     & q2L*q1tk1t*q2tk1t*P5*m2 - 8*U(vb)*q2L*q1tk1t*q2tk2t*P5*m1
     &  - 8*U(vb)*q2L*q1tk2t*q2tk1t*P6*m2
      AS2(4,vb) = AS2(4,vb) - 8*U(vb)*q2L*q1tk2t*q2tk1t*P8*m2 + 8*
     & U(vb)*q2L*q1tk2t*q2tk2t*P6*m1 + 8*U(vb)*q2L*q1tk2t*q2tk2t*
     & P8*m1 + 4*U(vb)*q2L*q2tk1t*t1*P5*m2 - 4*U(vb)*q2L*q2tk1t*t1*
     & P6*m2 - 4*U(vb)*q2L*q2tk1t*t1*P8*m2 - 4*U(vb)*q2L*q2tk2t*t1*
     & P5*m1 + 4*U(vb)*q2L*q2tk2t*t1*P6*m1 + 4*U(vb)*q2L*q2tk2t*t1*
     & P8*m1 - 4*U(vb)*q1tq2t*q1tL*P5*m1*q2k2 + 4*U(vb)*q1tq2t*q1tL
     & *P5*m2*q2k1 - 8*U(vb)*q1tq2t*q1tL*P6*m1*q1k2 - 4*U(vb)*q1tq2t
     & *q1tL*P6*m1*q2k2 + 8*U(vb)*q1tq2t*q1tL*P6*m2*q1k1 + 4*U(vb)*
     & q1tq2t*q1tL*P6*m2*q2k1 - 4*U(vb)*q1tq2t*q1tL*P8*m1*q2k2 + 4*
     & U(vb)*q1tq2t*q1tL*P8*m2*q2k1 + 4*U(vb)*q1tq2t*q2tL*P5*m1*
     & q1k2 + 8*U(vb)*q1tq2t*q2tL*P5*m1*q2k2 - 4*U(vb)*q1tq2t*q2tL*
     & P5*m2*q1k1 - 8*U(vb)*q1tq2t*q2tL*P5*m2*q2k1 + 4*U(vb)*q1tq2t*
     & q2tL*P6*m1*q1k2 - 4*U(vb)*q1tq2t*q2tL*P6*m2*q1k1 + 4*U(vb)*
     & q1tq2t*q2tL*P8*m1*q1k2 - 4*U(vb)*q1tq2t*q2tL*P8*m2*q1k1 + 4*
     & U(vb)*q1tk1t*q2tL*t1*P6*m2 + 4*U(vb)*q1tk1t*q2tL*t2*P5*m2 -
     & 8*U(vb)*q1tk1t*q2tL*P5*m1*q2k2
      AS2(4,vb) = AS2(4,vb) + 4*U(vb)*q1tk1t*q2tL*P5*m2*q1q2 + 4*
     & U(vb)*q1tk1t*q2tL*P6*m2*q1q2 - 8*U(vb)*q1tk1t*q2tL*P6*m2*
     & q1k2 + 4*U(vb)*q1tk1t*q2tL*P8*m2*q1q2 - 4*U(vb)*q1tk2t*q2tL*
     & t1*P6*m1 - 4*U(vb)*q1tk2t*q2tL*t2*P5*m1 - 4*U(vb)*q1tk2t*
     & q2tL*P5*m1*q1q2 + 8*U(vb)*q1tk2t*q2tL*P5*m1*q2k1 - 4*U(vb)*
     & q1tk2t*q2tL*P6*m1*q1q2 + 8*U(vb)*q1tk2t*q2tL*P6*m1*q1k2 + 8*
     & U(vb)*q1tk2t*q2tL*P6*m1*q2k2 - 8*U(vb)*q1tk2t*q2tL*P6*m2*
     & q2k1 - 4*U(vb)*q1tk2t*q2tL*P8*m1*q1q2 + 8*U(vb)*q1tk2t*q2tL*
     & P8*m1*q2k2 - 8*U(vb)*q1tk2t*q2tL*P8*m2*q2k1 - 4*U(vb)*q1tL*
     & q2tk1t*t1*P6*m2 - 4*U(vb)*q1tL*q2tk1t*t2*P5*m2 - 8*U(vb)*
     & q1tL*q2tk1t*P5*m1*q1k2 - 4*U(vb)*q1tL*q2tk1t*P5*m2*q1q2 + 8*
     & U(vb)*q1tL*q2tk1t*P5*m2*q1k1 + 8*U(vb)*q1tL*q2tk1t*P5*m2*
     & q2k1 - 4*U(vb)*q1tL*q2tk1t*P6*m2*q1q2 + 8*U(vb)*q1tL*q2tk1t*
     & P6*m2*q1k2 - 8*U(vb)*q1tL*q2tk1t*P8*m1*q1k2 - 4*U(vb)*q1tL*
     & q2tk1t*P8*m2*q1q2 + 8*U(vb)*q1tL*q2tk1t*P8*m2*q1k1 + 4*U(vb)*
     & q1tL*q2tk2t*t1*P6*m1
      AS2(4,vb) = AS2(4,vb) + 4*U(vb)*q1tL*q2tk2t*t2*P5*m1 + 4*
     & U(vb)*q1tL*q2tk2t*P5*m1*q1q2 - 8*U(vb)*q1tL*q2tk2t*P5*m1*
     & q2k1 + 4*U(vb)*q1tL*q2tk2t*P6*m1*q1q2 - 8*U(vb)*q1tL*q2tk2t*
     & P6*m2*q1k1 + 4*U(vb)*q1tL*q2tk2t*P8*m1*q1q2 - 4*U(vb)*q1tL*t2
     & *P5*m1*q1k2 + 4*U(vb)*q1tL*t2*P5*m2*q1k1 + 4*U(vb)*q1tL*t2*P6*
     & m1*q1k2 - 4*U(vb)*q1tL*t2*P6*m2*q1k1 - 4*U(vb)*q1tL*t2*P8*m1*
     & q1k2 + 4*U(vb)*q1tL*t2*P8*m2*q1k1 - 4*U(vb)*q2tL*t1*P5*m1*q2k2
     &  + 4*U(vb)*q2tL*t1*P5*m2*q2k1 + 4*U(vb)*q2tL*t1*P6*m1*q2k2 - 4
     & *U(vb)*q2tL*t1*P6*m2*q2k1 + 4*U(vb)*q2tL*t1*P8*m1*q2k2 - 4*
     & U(vb)*q2tL*t1*P8*m2*q2k1 - 4*L(vb)*q2U*q1tq2t*q1tk1t*P5*m2
     &  + 4*L(vb)*q2U*q1tq2t*q1tk1t*P6*m2 - 4*L(vb)*q2U*q1tq2t*
     & q1tk1t*P8*m2 + 4*L(vb)*q2U*q1tq2t*q1tk2t*P5*m1 - 4*L(vb)*
     & q2U*q1tq2t*q1tk2t*P6*m1 + 4*L(vb)*q2U*q1tq2t*q1tk2t*P8*m1
     &  + 4*L(vb)*q2U*q1tq2t*q2tk1t*P5*m2 - 4*L(vb)*q2U*q1tq2t*
     & q2tk1t*P6*m2 - 4*L(vb)*q2U*q1tq2t*q2tk1t*P8*m2 - 4*L(vb)*
     & q2U*q1tq2t*q2tk2t*P5*m1
      AS2(4,vb) = AS2(4,vb) + 4*L(vb)*q2U*q1tq2t*q2tk2t*P6*m1 + 4*
     & L(vb)*q2U*q1tq2t*q2tk2t*P8*m1 + 4*L(vb)*q2U*q1tq2t*P5*m1*
     & q1k2 + 4*L(vb)*q2U*q1tq2t*P5*m1*q2k2 - 4*L(vb)*q2U*q1tq2t*P5
     & *m2*q1k1 - 4*L(vb)*q2U*q1tq2t*P5*m2*q2k1 + 4*L(vb)*q2U*
     & q1tq2t*P6*m1*q1k2 + 4*L(vb)*q2U*q1tq2t*P6*m1*q2k2 - 4*L(vb)*
     & q2U*q1tq2t*P6*m2*q1k1 - 4*L(vb)*q2U*q1tq2t*P6*m2*q2k1 + 4*
     & L(vb)*q2U*q1tq2t*P8*m1*q1k2 + 4*L(vb)*q2U*q1tq2t*P8*m1*q2k2
     &  - 4*L(vb)*q2U*q1tq2t*P8*m2*q1k1 - 4*L(vb)*q2U*q1tq2t*P8*m2*
     & q2k1 + 8*L(vb)*q2U*q1tk1t*q2tk1t*P8*m2 + 8*L(vb)*q2U*q1tk1t
     & *q2tk2t*P5*m1 - 8*L(vb)*q2U*q1tk1t*q2tk2t*P6*m2 + 4*L(vb)*
     & q2U*q1tk1t*t2*P5*m2 - 4*L(vb)*q2U*q1tk1t*t2*P6*m2 + 4*L(vb)*
     & q2U*q1tk1t*t2*P8*m2 - 8*L(vb)*q2U*q1tk2t*q2tk1t*P5*m1 + 8*
     & L(vb)*q2U*q1tk2t*q2tk1t*P6*m2 - 8*L(vb)*q2U*q1tk2t*q2tk1t*
     & P8*m1 + 8*L(vb)*q2U*q1tk2t*q2tk1t*P8*m2 - 8*L(vb)*q2U*
     & q1tk2t*q2tk2t*P8*m1 - 4*L(vb)*q2U*q1tk2t*t2*P5*m1 + 4*L(vb)*
     & q2U*q1tk2t*t2*P6*m1
      AS2(4,vb) = AS2(4,vb) - 4*L(vb)*q2U*q1tk2t*t2*P8*m1 - 4*L(vb)
     & *q2U*q2tk1t*t1*P5*m2 + 4*L(vb)*q2U*q2tk1t*t1*P6*m2 + 4*L(vb)
     & *q2U*q2tk1t*t1*P8*m2 + 4*L(vb)*q2U*q2tk2t*t1*P5*m1 - 4*L(vb)
     & *q2U*q2tk2t*t1*P6*m1 - 4*L(vb)*q2U*q2tk2t*t1*P8*m1 + 4*L(vb)
     & *q1tq2t*q1tU*P5*m1*q2k2 - 4*L(vb)*q1tq2t*q1tU*P5*m2*q2k1 + 8
     & *L(vb)*q1tq2t*q1tU*P6*m1*q1k2 + 4*L(vb)*q1tq2t*q1tU*P6*m1*
     & q2k2 - 8*L(vb)*q1tq2t*q1tU*P6*m2*q1k1 - 4*L(vb)*q1tq2t*q1tU*
     & P6*m2*q2k1 + 4*L(vb)*q1tq2t*q1tU*P8*m1*q2k2 - 4*L(vb)*q1tq2t*
     & q1tU*P8*m2*q2k1 - 4*L(vb)*q1tq2t*q2tU*P5*m1*q1k2 - 8*L(vb)*
     & q1tq2t*q2tU*P5*m1*q2k2 + 4*L(vb)*q1tq2t*q2tU*P5*m2*q1k1 + 8*
     & L(vb)*q1tq2t*q2tU*P5*m2*q2k1 - 4*L(vb)*q1tq2t*q2tU*P6*m1*
     & q1k2 + 4*L(vb)*q1tq2t*q2tU*P6*m2*q1k1 - 4*L(vb)*q1tq2t*q2tU*
     & P8*m1*q1k2 + 4*L(vb)*q1tq2t*q2tU*P8*m2*q1k1 - 4*L(vb)*q1tk1t*
     & q2tU*t1*P6*m2 - 4*L(vb)*q1tk1t*q2tU*t2*P5*m2 + 8*L(vb)*
     & q1tk1t*q2tU*P5*m1*q2k2 - 4*L(vb)*q1tk1t*q2tU*P5*m2*q1q2 - 4*
     & L(vb)*q1tk1t*q2tU*P6*m2*q1q2
      AS2(4,vb) = AS2(4,vb) + 8*L(vb)*q1tk1t*q2tU*P6*m2*q1k2 - 4*
     & L(vb)*q1tk1t*q2tU*P8*m2*q1q2 + 4*L(vb)*q1tk2t*q2tU*t1*P6*m1
     &  + 4*L(vb)*q1tk2t*q2tU*t2*P5*m1 + 4*L(vb)*q1tk2t*q2tU*P5*m1*
     & q1q2 - 8*L(vb)*q1tk2t*q2tU*P5*m1*q2k1 + 4*L(vb)*q1tk2t*q2tU*
     & P6*m1*q1q2 - 8*L(vb)*q1tk2t*q2tU*P6*m1*q1k2 - 8*L(vb)*q1tk2t*
     & q2tU*P6*m1*q2k2 + 8*L(vb)*q1tk2t*q2tU*P6*m2*q2k1 + 4*L(vb)*
     & q1tk2t*q2tU*P8*m1*q1q2 - 8*L(vb)*q1tk2t*q2tU*P8*m1*q2k2 + 8*
     & L(vb)*q1tk2t*q2tU*P8*m2*q2k1 + 4*L(vb)*q1tU*q2tk1t*t1*P6*m2
     &  + 4*L(vb)*q1tU*q2tk1t*t2*P5*m2 + 8*L(vb)*q1tU*q2tk1t*P5*m1*
     & q1k2 + 4*L(vb)*q1tU*q2tk1t*P5*m2*q1q2 - 8*L(vb)*q1tU*q2tk1t*
     & P5*m2*q1k1 - 8*L(vb)*q1tU*q2tk1t*P5*m2*q2k1 + 4*L(vb)*q1tU*
     & q2tk1t*P6*m2*q1q2 - 8*L(vb)*q1tU*q2tk1t*P6*m2*q1k2 + 8*L(vb)*
     & q1tU*q2tk1t*P8*m1*q1k2 + 4*L(vb)*q1tU*q2tk1t*P8*m2*q1q2 - 8*
     & L(vb)*q1tU*q2tk1t*P8*m2*q1k1 - 4*L(vb)*q1tU*q2tk2t*t1*P6*m1
     &  - 4*L(vb)*q1tU*q2tk2t*t2*P5*m1 - 4*L(vb)*q1tU*q2tk2t*P5*m1*
     & q1q2
      AS2(4,vb) = AS2(4,vb) + 8*L(vb)*q1tU*q2tk2t*P5*m1*q2k1 - 4*
     & L(vb)*q1tU*q2tk2t*P6*m1*q1q2 + 8*L(vb)*q1tU*q2tk2t*P6*m2*
     & q1k1 - 4*L(vb)*q1tU*q2tk2t*P8*m1*q1q2 + 4*L(vb)*q1tU*t2*P5*m1
     & *q1k2 - 4*L(vb)*q1tU*t2*P5*m2*q1k1 - 4*L(vb)*q1tU*t2*P6*m1*
     & q1k2 + 4*L(vb)*q1tU*t2*P6*m2*q1k1 + 4*L(vb)*q1tU*t2*P8*m1*q1k2
     &  - 4*L(vb)*q1tU*t2*P8*m2*q1k1 + 4*L(vb)*q2tU*t1*P5*m1*q2k2 - 4
     & *L(vb)*q2tU*t1*P5*m2*q2k1 - 4*L(vb)*q2tU*t1*P6*m1*q2k2 + 4*
     & L(vb)*q2tU*t1*P6*m2*q2k1 - 4*L(vb)*q2tU*t1*P8*m1*q2k2 + 4*
     & L(vb)*q2tU*t1*P8*m2*q2k1

       A58(1,vb)=0.
       A58(2,vb)=0.
       A58(3,vb)=0.
       A58(4,vb)=0.

       A51(1,vb)=0.
       A51(2,vb)=0.
       A51(3,vb)=0.
       A51(4,vb)=0.

       A52(1,vb)=0.
       A52(2,vb)=0.
       A52(3,vb)=0.
       A52(4,vb)=0.

      A58(1,vb) =
     &  - 8*k1(vb)*q1tq2t*q1tk1t*P1*m2 - 8*k1(vb)*q1tq2t*q1tk1t*P2*
     & m2 + 8*k1(vb)*q1tq2t*q1tk2t*P1*m1 - 8*k1(vb)*q1tq2t*q1tk2t*
     & P2*m1 + 16*k1(vb)*q1tq2t*q1tk2t*P2*m2 + 8*k1(vb)*q1tq2t*
     & q2tk1t*P1*m2 + 8*k1(vb)*q1tq2t*q2tk1t*P2*m2 - 8*k1(vb)*
     & q1tq2t*q2tk2t*P1*m1 + 8*k1(vb)*q1tq2t*q2tk2t*P2*m1 - 16*
     & k1(vb)*q1tq2t*q2tk2t*P2*m2 + 4*k1(vb)*q1tq2t*P1*m1*q1k2 - 4*
     & k1(vb)*q1tq2t*P1*m1*q2k2 - 4*k1(vb)*q1tq2t*P1*m2*q1k1 + 4*
     & k1(vb)*q1tq2t*P1*m2*q2k1 - 4*k1(vb)*q1tq2t*P2*m1*q1k2 + 4*
     & k1(vb)*q1tq2t*P2*m1*q2k2 - 4*k1(vb)*q1tq2t*P2*m2*q1k1 + 8*
     & k1(vb)*q1tq2t*P2*m2*q1k2 + 4*k1(vb)*q1tq2t*P2*m2*q2k1 - 8*
     & k1(vb)*q1tq2t*P2*m2*q2k2 + 4*k1(vb)*q1tk1t*t2*P1*m2 + 4*k1(vb)
     & *q1tk1t*t2*P2*m2 - 4*k1(vb)*q1tk2t*t2*P1*m1 + 4*k1(vb)*q1tk2t
     & *t2*P2*m1 - 8*k1(vb)*q1tk2t*t2*P2*m2 - 4*k1(vb)*q2tk1t*t1*P1*
     & m2 - 4*k1(vb)*q2tk1t*t1*P2*m2 + 4*k1(vb)*q2tk2t*t1*P1*m1 - 4*
     & k1(vb)*q2tk2t*t1*P2*m1
      A58(1,vb) = A58(1,vb) + 8*k1(vb)*q2tk2t*t1*P2*m2 - 16*k2(vb)*
     & q1tq2t*q1tk1t*P1*m1 + 8*k2(vb)*q1tq2t*q1tk1t*P1*m2 - 8*
     & k2(vb)*q1tq2t*q1tk1t*P2*m2 + 8*k2(vb)*q1tq2t*q1tk2t*P1*m1 +
     & 8*k2(vb)*q1tq2t*q1tk2t*P2*m1 + 16*k2(vb)*q1tq2t*q2tk1t*P1*m1
     &  - 8*k2(vb)*q1tq2t*q2tk1t*P1*m2 + 8*k2(vb)*q1tq2t*q2tk1t*P2*
     & m2 - 8*k2(vb)*q1tq2t*q2tk2t*P1*m1 - 8*k2(vb)*q1tq2t*q2tk2t*
     & P2*m1 - 8*k2(vb)*q1tq2t*P1*m1*q1k1 + 4*k2(vb)*q1tq2t*P1*m1*
     & q1k2 + 8*k2(vb)*q1tq2t*P1*m1*q2k1 - 4*k2(vb)*q1tq2t*P1*m1*q2k2
     &  + 4*k2(vb)*q1tq2t*P1*m2*q1k1 - 4*k2(vb)*q1tq2t*P1*m2*q2k1 + 4
     & *k2(vb)*q1tq2t*P2*m1*q1k2 - 4*k2(vb)*q1tq2t*P2*m1*q2k2 - 4*
     & k2(vb)*q1tq2t*P2*m2*q1k1 + 4*k2(vb)*q1tq2t*P2*m2*q2k1 + 8*
     & k2(vb)*q1tk1t*t2*P1*m1 - 4*k2(vb)*q1tk1t*t2*P1*m2 + 4*k2(vb)*
     & q1tk1t*t2*P2*m2 - 4*k2(vb)*q1tk2t*t2*P1*m1 - 4*k2(vb)*q1tk2t*
     & t2*P2*m1 - 8*k2(vb)*q2tk1t*t1*P1*m1 + 4*k2(vb)*q2tk1t*t1*P1*m2
     &  - 4*k2(vb)*q2tk1t*t1*P2*m2 + 4*k2(vb)*q2tk2t*t1*P1*m1 + 4*
     & k2(vb)*q2tk2t*t1*P2*m1
      A58(1,vb) = A58(1,vb) + 8*kw(vb)*q1tq2t*q1tk1t*P1*m2 - 8*kw(vb)
     & *q1tq2t*q1tk1t*P2*m2 + 8*kw(vb)*q1tq2t*q1tk2t*P1*m1 - 8*
     & kw(vb)*q1tq2t*q1tk2t*P2*m1 - 8*kw(vb)*q1tq2t*q2tk1t*P1*m2 +
     & 8*kw(vb)*q1tq2t*q2tk1t*P2*m2 - 8*kw(vb)*q1tq2t*q2tk2t*P1*m1
     &  + 8*kw(vb)*q1tq2t*q2tk2t*P2*m1 + 4*kw(vb)*q1tq2t*P1*m1*q1k2
     &  - 4*kw(vb)*q1tq2t*P1*m1*q2k2 + 4*kw(vb)*q1tq2t*P1*m2*q1k1 - 4
     & *kw(vb)*q1tq2t*P1*m2*q2k1 - 4*kw(vb)*q1tq2t*P2*m1*q1k2 + 4*
     & kw(vb)*q1tq2t*P2*m1*q2k2 - 4*kw(vb)*q1tq2t*P2*m2*q1k1 + 4*
     & kw(vb)*q1tq2t*P2*m2*q2k1 - 4*kw(vb)*q1tk1t*t2*P1*m2 + 4*kw(vb)
     & *q1tk1t*t2*P2*m2 - 4*kw(vb)*q1tk2t*t2*P1*m1 + 4*kw(vb)*q1tk2t
     & *t2*P2*m1 + 4*kw(vb)*q2tk1t*t1*P1*m2 - 4*kw(vb)*q2tk1t*t1*P2*
     & m2 + 4*kw(vb)*q2tk2t*t1*P1*m1 - 4*kw(vb)*q2tk2t*t1*P2*m1 - 4*
     & q1(vb)*q1tq2t*P1*m1*q1k2 - 4*q1(vb)*q1tq2t*P1*m1*q2k2 - 4*
     & q1(vb)*q1tq2t*P1*m2*q1k1 - 4*q1(vb)*q1tq2t*P1*m2*q2k1 + 4*
     & q1(vb)*q1tq2t*P2*m1*q1k2 + 4*q1(vb)*q1tq2t*P2*m1*q2k2 + 4*
     & q1(vb)*q1tq2t*P2*m2*q1k1
      A58(1,vb) = A58(1,vb) + 4*q1(vb)*q1tq2t*P2*m2*q2k1 + 4*q2(vb)*
     & q1tq2t*P1*m1*q1k2 + 4*q2(vb)*q1tq2t*P1*m1*q2k2 + 4*q2(vb)*
     & q1tq2t*P1*m2*q1k1 + 4*q2(vb)*q1tq2t*P1*m2*q2k1 - 4*q2(vb)*
     & q1tq2t*P2*m1*q1k2 - 4*q2(vb)*q1tq2t*P2*m1*q2k2 - 4*q2(vb)*
     & q1tq2t*P2*m2*q1k1 - 4*q2(vb)*q1tq2t*P2*m2*q2k1 - 4*q1t(vb)*t2*
     & P1*m1*q1k2 - 4*q1t(vb)*t2*P1*m1*q2k2 - 4*q1t(vb)*t2*P1*m2*q1k1
     &  - 4*q1t(vb)*t2*P1*m2*q2k1 + 4*q1t(vb)*t2*P2*m1*q1k2 + 4*q1t(vb)
     & *t2*P2*m1*q2k2 + 4*q1t(vb)*t2*P2*m2*q1k1 + 4*q1t(vb)*t2*P2*m2*
     & q2k1 + 8*q1t(vb)*q1tq2t*P1*m1*q1k2 + 8*q1t(vb)*q1tq2t*P1*m1*
     & q2k2 + 8*q1t(vb)*q1tq2t*P1*m2*q1k1 + 8*q1t(vb)*q1tq2t*P1*m2*
     & q2k1 - 8*q1t(vb)*q1tq2t*P2*m1*q1k2 - 8*q1t(vb)*q1tq2t*P2*m1*
     & q2k2 - 8*q1t(vb)*q1tq2t*P2*m2*q1k1 - 8*q1t(vb)*q1tq2t*P2*m2*
     & q2k1 + 4*q2t(vb)*t1*P1*m1*q1k2 + 4*q2t(vb)*t1*P1*m1*q2k2 + 4*
     & q2t(vb)*t1*P1*m2*q1k1 + 4*q2t(vb)*t1*P1*m2*q2k1 - 4*q2t(vb)*t1*
     & P2*m1*q1k2 - 4*q2t(vb)*t1*P2*m1*q2k2 - 4*q2t(vb)*t1*P2*m2*q1k1
     &  - 4*q2t(vb)*t1*P2*m2*q2k1
      A58(1,vb) = A58(1,vb) - 8*q2t(vb)*q1tq2t*P1*m1*q1k2 - 8*q2t(vb)*
     & q1tq2t*P1*m1*q2k2 - 8*q2t(vb)*q1tq2t*P1*m2*q1k1 - 8*q2t(vb)*
     & q1tq2t*P1*m2*q2k1 + 8*q2t(vb)*q1tq2t*P2*m1*q1k2 + 8*q2t(vb)*
     & q1tq2t*P2*m1*q2k2 + 8*q2t(vb)*q1tq2t*P2*m2*q1k1 + 8*q2t(vb)*
     & q1tq2t*P2*m2*q2k1

      A58(2,vb) =
     & 8*k1(vb)*q2U*q1tq2t*P1*m1*m2 + 8*k1(vb)*q2U*q1tq2t*P1*q1k2
     &  + 8*k1(vb)*q2U*q1tq2t*P1*q2k2 - 8*k1(vb)*q2U*q1tq2t*P1*k1k2
     &  - 8*k1(vb)*q2U*q1tq2t*P2*m1*m2 + 8*k1(vb)*q2U*q1tq2t*P2*
     & q1k2 + 8*k1(vb)*q2U*q1tq2t*P2*q2k2 + 8*k1(vb)*q2U*q1tq2t*P2*
     & k1k2 + 8*k1(vb)*q1tq2t*q1tU*P1*m1*m2 + 8*k1(vb)*q1tq2t*q1tU*
     & P1*q1k2 + 8*k1(vb)*q1tq2t*q1tU*P1*q2k2 - 8*k1(vb)*q1tq2t*
     & q1tU*P1*k1k2 - 8*k1(vb)*q1tq2t*q1tU*P2*m1*m2 + 8*k1(vb)*
     & q1tq2t*q1tU*P2*q1k2 + 8*k1(vb)*q1tq2t*q1tU*P2*q2k2 + 8*
     & k1(vb)*q1tq2t*q1tU*P2*k1k2 - 8*k1(vb)*q1tq2t*q2tU*P1*m1*m2
     &  - 8*k1(vb)*q1tq2t*q2tU*P1*q1k2 - 8*k1(vb)*q1tq2t*q2tU*P1*
     & q2k2 + 8*k1(vb)*q1tq2t*q2tU*P1*k1k2 + 8*k1(vb)*q1tq2t*q2tU*
     & P2*m1*m2 - 8*k1(vb)*q1tq2t*q2tU*P2*q1k2 - 8*k1(vb)*q1tq2t*
     & q2tU*P2*q2k2 - 8*k1(vb)*q1tq2t*q2tU*P2*k1k2 - 4*k1(vb)*q1tU*
     & t2*P1*m1*m2 - 4*k1(vb)*q1tU*t2*P1*q1k2 - 4*k1(vb)*q1tU*t2*P1*
     & q2k2
      A58(2,vb) = A58(2,vb) + 4*k1(vb)*q1tU*t2*P1*k1k2 + 4*k1(vb)*
     & q1tU*t2*P2*m1*m2 - 4*k1(vb)*q1tU*t2*P2*q1k2 - 4*k1(vb)*q1tU*
     & t2*P2*q2k2 - 4*k1(vb)*q1tU*t2*P2*k1k2 + 4*k1(vb)*q2tU*t1*P1*m1
     & *m2 + 4*k1(vb)*q2tU*t1*P1*q1k2 + 4*k1(vb)*q2tU*t1*P1*q2k2 - 4*
     & k1(vb)*q2tU*t1*P1*k1k2 - 4*k1(vb)*q2tU*t1*P2*m1*m2 + 4*k1(vb)*
     & q2tU*t1*P2*q1k2 + 4*k1(vb)*q2tU*t1*P2*q2k2 + 4*k1(vb)*q2tU*t1
     & *P2*k1k2 + 8*k2(vb)*q2U*q1tq2t*P1*m1*m2 - 8*k2(vb)*q2U*
     & q1tq2t*P1*q1k1 - 8*k2(vb)*q2U*q1tq2t*P1*q2k1 - 8*k2(vb)*q2U*
     & q1tq2t*P1*k1k2 - 8*k2(vb)*q2U*q1tq2t*P2*m1*m2 - 8*k2(vb)*q2U
     & *q1tq2t*P2*q1k1 - 8*k2(vb)*q2U*q1tq2t*P2*q2k1 + 8*k2(vb)*q2U
     & *q1tq2t*P2*k1k2 + 8*k2(vb)*q1tq2t*q1tU*P1*m1*m2 - 8*k2(vb)*
     & q1tq2t*q1tU*P1*q1k1 - 8*k2(vb)*q1tq2t*q1tU*P1*q2k1 - 8*
     & k2(vb)*q1tq2t*q1tU*P1*k1k2 - 8*k2(vb)*q1tq2t*q1tU*P2*m1*m2
     &  - 8*k2(vb)*q1tq2t*q1tU*P2*q1k1 - 8*k2(vb)*q1tq2t*q1tU*P2*
     & q2k1 + 8*k2(vb)*q1tq2t*q1tU*P2*k1k2 - 8*k2(vb)*q1tq2t*q2tU*
     & P1*m1*m2
      A58(2,vb) = A58(2,vb) + 8*k2(vb)*q1tq2t*q2tU*P1*q1k1 + 8*k2(vb)
     & *q1tq2t*q2tU*P1*q2k1 + 8*k2(vb)*q1tq2t*q2tU*P1*k1k2 + 8*
     & k2(vb)*q1tq2t*q2tU*P2*m1*m2 + 8*k2(vb)*q1tq2t*q2tU*P2*q1k1
     &  + 8*k2(vb)*q1tq2t*q2tU*P2*q2k1 - 8*k2(vb)*q1tq2t*q2tU*P2*
     & k1k2 - 4*k2(vb)*q1tU*t2*P1*m1*m2 + 4*k2(vb)*q1tU*t2*P1*q1k1 +
     & 4*k2(vb)*q1tU*t2*P1*q2k1 + 4*k2(vb)*q1tU*t2*P1*k1k2 + 4*k2(vb)
     & *q1tU*t2*P2*m1*m2 + 4*k2(vb)*q1tU*t2*P2*q1k1 + 4*k2(vb)*q1tU*
     & t2*P2*q2k1 - 4*k2(vb)*q1tU*t2*P2*k1k2 + 4*k2(vb)*q2tU*t1*P1*m1
     & *m2 - 4*k2(vb)*q2tU*t1*P1*q1k1 - 4*k2(vb)*q2tU*t1*P1*q2k1 - 4*
     & k2(vb)*q2tU*t1*P1*k1k2 - 4*k2(vb)*q2tU*t1*P2*m1*m2 - 4*k2(vb)*
     & q2tU*t1*P2*q1k1 - 4*k2(vb)*q2tU*t1*P2*q2k1 + 4*k2(vb)*q2tU*t1
     & *P2*k1k2 + 8*kw(vb)*q2U*q1tq2t*P1*m1*m2 - 8*kw(vb)*q2U*
     & q1tq2t*P1*k1k2 - 8*kw(vb)*q2U*q1tq2t*P2*m1*m2 + 8*kw(vb)*q2U
     & *q1tq2t*P2*k1k2 + 8*kw(vb)*q1tq2t*q1tU*P1*m1*m2 - 8*kw(vb)*
     & q1tq2t*q1tU*P1*k1k2 - 8*kw(vb)*q1tq2t*q1tU*P2*m1*m2 + 8*
     & kw(vb)*q1tq2t*q1tU*P2*k1k2
      A58(2,vb) = A58(2,vb) - 8*kw(vb)*q1tq2t*q2tU*P1*m1*m2 + 8*
     & kw(vb)*q1tq2t*q2tU*P1*k1k2 + 8*kw(vb)*q1tq2t*q2tU*P2*m1*m2
     &  - 8*kw(vb)*q1tq2t*q2tU*P2*k1k2 - 4*kw(vb)*q1tU*t2*P1*m1*m2
     &  + 4*kw(vb)*q1tU*t2*P1*k1k2 + 4*kw(vb)*q1tU*t2*P2*m1*m2 - 4*
     & kw(vb)*q1tU*t2*P2*k1k2 + 4*kw(vb)*q2tU*t1*P1*m1*m2 - 4*kw(vb)*
     & q2tU*t1*P1*k1k2 - 4*kw(vb)*q2tU*t1*P2*m1*m2 + 4*kw(vb)*q2tU*
     & t1*P2*k1k2 + 16*U(vb)*q1tq2t*q1tk1t*P1*m1*m2 + 8*U(vb)*q1tq2t
     & *q1tk1t*P1*q1k2 + 8*U(vb)*q1tq2t*q1tk1t*P1*q2k2 - 16*U(vb)*
     & q1tq2t*q1tk1t*P1*k1k2 + 8*U(vb)*q1tq2t*q1tk1t*P2*q1k2 + 8*
     & U(vb)*q1tq2t*q1tk1t*P2*q2k2 - 8*U(vb)*q1tq2t*q1tk2t*P1*q1k1
     &  - 8*U(vb)*q1tq2t*q1tk2t*P1*q2k1 - 16*U(vb)*q1tq2t*q1tk2t*P2
     & *m1*m2 - 8*U(vb)*q1tq2t*q1tk2t*P2*q1k1 - 8*U(vb)*q1tq2t*
     & q1tk2t*P2*q2k1 + 16*U(vb)*q1tq2t*q1tk2t*P2*k1k2 - 16*U(vb)*
     & q1tq2t*q2tk1t*P1*m1*m2 - 8*U(vb)*q1tq2t*q2tk1t*P1*q1k2 - 8*
     & U(vb)*q1tq2t*q2tk1t*P1*q2k2 + 16*U(vb)*q1tq2t*q2tk1t*P1*k1k2
     &  - 8*U(vb)*q1tq2t*q2tk1t*P2*q1k2
      A58(2,vb) = A58(2,vb) - 8*U(vb)*q1tq2t*q2tk1t*P2*q2k2 + 8*U(vb)
     & *q1tq2t*q2tk2t*P1*q1k1 + 8*U(vb)*q1tq2t*q2tk2t*P1*q2k1 + 16*
     & U(vb)*q1tq2t*q2tk2t*P2*m1*m2 + 8*U(vb)*q1tq2t*q2tk2t*P2*q1k1
     &  + 8*U(vb)*q1tq2t*q2tk2t*P2*q2k1 - 16*U(vb)*q1tq2t*q2tk2t*P2
     & *k1k2 + 8*U(vb)*q1tq2t*P1*m1*m2*q1k1 - 8*U(vb)*q1tq2t*P1*m1*m2
     & *q2k1 + 8*U(vb)*q1tq2t*P1*q1k1*q2k2 - 8*U(vb)*q1tq2t*P1*q1k1*
     & k1k2 - 8*U(vb)*q1tq2t*P1*q1k2*q2k1 + 8*U(vb)*q1tq2t*P1*q2k1*
     & k1k2 - 8*U(vb)*q1tq2t*P2*m1*m2*q1k2 + 8*U(vb)*q1tq2t*P2*m1*m2*
     & q2k2 + 8*U(vb)*q1tq2t*P2*q1k1*q2k2 - 8*U(vb)*q1tq2t*P2*q1k2*
     & q2k1 + 8*U(vb)*q1tq2t*P2*q1k2*k1k2 - 8*U(vb)*q1tq2t*P2*q2k2*
     & k1k2 - 8*U(vb)*q1tk1t*t2*P1*m1*m2 - 4*U(vb)*q1tk1t*t2*P1*q1k2
     &  - 4*U(vb)*q1tk1t*t2*P1*q2k2 + 8*U(vb)*q1tk1t*t2*P1*k1k2 - 4*
     & U(vb)*q1tk1t*t2*P2*q1k2 - 4*U(vb)*q1tk1t*t2*P2*q2k2 + 4*U(vb)*
     & q1tk2t*t2*P1*q1k1 + 4*U(vb)*q1tk2t*t2*P1*q2k1 + 8*U(vb)*
     & q1tk2t*t2*P2*m1*m2 + 4*U(vb)*q1tk2t*t2*P2*q1k1 + 4*U(vb)*
     & q1tk2t*t2*P2*q2k1
      A58(2,vb) = A58(2,vb) - 8*U(vb)*q1tk2t*t2*P2*k1k2 + 8*U(vb)*
     & q2tk1t*t1*P1*m1*m2 + 4*U(vb)*q2tk1t*t1*P1*q1k2 + 4*U(vb)*
     & q2tk1t*t1*P1*q2k2 - 8*U(vb)*q2tk1t*t1*P1*k1k2 + 4*U(vb)*
     & q2tk1t*t1*P2*q1k2 + 4*U(vb)*q2tk1t*t1*P2*q2k2 - 4*U(vb)*
     & q2tk2t*t1*P1*q1k1 - 4*U(vb)*q2tk2t*t1*P1*q2k1 - 8*U(vb)*
     & q2tk2t*t1*P2*m1*m2 - 4*U(vb)*q2tk2t*t1*P2*q1k1 - 4*U(vb)*
     & q2tk2t*t1*P2*q2k1 + 8*U(vb)*q2tk2t*t1*P2*k1k2

      A58(3,vb) =
     & 8*k1(vb)*q1L*q1tq2t*q1tk2t*P1 + 8*k1(vb)*q1L*q1tq2t*q1tk2t
     & *P2 - 8*k1(vb)*q1L*q1tq2t*q2tk2t*P1 - 8*k1(vb)*q1L*q1tq2t*
     & q2tk2t*P2 - 4*k1(vb)*q1L*q1tq2t*P1*m1*m2 - 8*k1(vb)*q1L*
     & q1tq2t*P1*q2k2 + 4*k1(vb)*q1L*q1tq2t*P1*k1k2 + 4*k1(vb)*q1L*
     & q1tq2t*P2*m1*m2 - 8*k1(vb)*q1L*q1tq2t*P2*q2k2 - 4*k1(vb)*q1L
     & *q1tq2t*P2*k1k2 - 4*k1(vb)*q1L*q1tk2t*t2*P1 - 4*k1(vb)*q1L*
     & q1tk2t*t2*P2 + 4*k1(vb)*q1L*q2tk2t*t1*P1 + 4*k1(vb)*q1L*
     & q2tk2t*t1*P2 + 8*k1(vb)*q2L*q1tq2t*q1tk2t*P1 + 8*k1(vb)*q2L
     & *q1tq2t*q1tk2t*P2 - 8*k1(vb)*q2L*q1tq2t*q2tk2t*P1 - 8*
     & k1(vb)*q2L*q1tq2t*q2tk2t*P2 + 4*k1(vb)*q2L*q1tq2t*P1*m1*m2
     &  + 8*k1(vb)*q2L*q1tq2t*P1*q1k2 - 4*k1(vb)*q2L*q1tq2t*P1*k1k2
     &  - 4*k1(vb)*q2L*q1tq2t*P2*m1*m2 + 8*k1(vb)*q2L*q1tq2t*P2*
     & q1k2 + 4*k1(vb)*q2L*q1tq2t*P2*k1k2 - 4*k1(vb)*q2L*q1tk2t*t2*
     & P1 - 4*k1(vb)*q2L*q1tk2t*t2*P2 + 4*k1(vb)*q2L*q2tk2t*t1*P1
     &  + 4*k1(vb)*q2L*q2tk2t*t1*P2
      A58(3,vb) = A58(3,vb) + 8*k1(vb)*q1tq2t*q1tL*P1*m1*m2 + 8*
     & k1(vb)*q1tq2t*q1tL*P1*q1k2 + 8*k1(vb)*q1tq2t*q1tL*P1*q2k2 - 
     & 8*k1(vb)*q1tq2t*q1tL*P1*k1k2 - 8*k1(vb)*q1tq2t*q1tL*P2*m1*m2
     &  + 8*k1(vb)*q1tq2t*q1tL*P2*q1k2 + 8*k1(vb)*q1tq2t*q1tL*P2*
     & q2k2 + 8*k1(vb)*q1tq2t*q1tL*P2*k1k2 - 8*k1(vb)*q1tq2t*q2tL*
     & P1*m1*m2 - 8*k1(vb)*q1tq2t*q2tL*P1*q1k2 - 8*k1(vb)*q1tq2t*
     & q2tL*P1*q2k2 + 8*k1(vb)*q1tq2t*q2tL*P1*k1k2 + 8*k1(vb)*
     & q1tq2t*q2tL*P2*m1*m2 - 8*k1(vb)*q1tq2t*q2tL*P2*q1k2 - 8*
     & k1(vb)*q1tq2t*q2tL*P2*q2k2 - 8*k1(vb)*q1tq2t*q2tL*P2*k1k2 - 
     & 4*k1(vb)*q1tL*t2*P1*m1*m2 - 4*k1(vb)*q1tL*t2*P1*q1k2 - 4*
     & k1(vb)*q1tL*t2*P1*q2k2 + 4*k1(vb)*q1tL*t2*P1*k1k2 + 4*k1(vb)*
     & q1tL*t2*P2*m1*m2 - 4*k1(vb)*q1tL*t2*P2*q1k2 - 4*k1(vb)*q1tL*
     & t2*P2*q2k2 - 4*k1(vb)*q1tL*t2*P2*k1k2 + 4*k1(vb)*q2tL*t1*P1*m1
     & *m2 + 4*k1(vb)*q2tL*t1*P1*q1k2 + 4*k1(vb)*q2tL*t1*P1*q2k2 - 4*
     & k1(vb)*q2tL*t1*P1*k1k2 - 4*k1(vb)*q2tL*t1*P2*m1*m2 + 4*k1(vb)*
     & q2tL*t1*P2*q1k2
      A58(3,vb) = A58(3,vb) + 4*k1(vb)*q2tL*t1*P2*q2k2 + 4*k1(vb)*
     & q2tL*t1*P2*k1k2 - 8*k2(vb)*q1L*q1tq2t*q1tk1t*P1 - 8*k2(vb)*
     & q1L*q1tq2t*q1tk1t*P2 + 8*k2(vb)*q1L*q1tq2t*q2tk1t*P1 + 8*
     & k2(vb)*q1L*q1tq2t*q2tk1t*P2 - 4*k2(vb)*q1L*q1tq2t*P1*m1*m2
     &  + 8*k2(vb)*q1L*q1tq2t*P1*q2k1 + 4*k2(vb)*q1L*q1tq2t*P1*k1k2
     &  + 4*k2(vb)*q1L*q1tq2t*P2*m1*m2 + 8*k2(vb)*q1L*q1tq2t*P2*
     & q2k1 - 4*k2(vb)*q1L*q1tq2t*P2*k1k2 + 4*k2(vb)*q1L*q1tk1t*t2*
     & P1 + 4*k2(vb)*q1L*q1tk1t*t2*P2 - 4*k2(vb)*q1L*q2tk1t*t1*P1
     &  - 4*k2(vb)*q1L*q2tk1t*t1*P2 - 8*k2(vb)*q2L*q1tq2t*q1tk1t*
     & P1 - 8*k2(vb)*q2L*q1tq2t*q1tk1t*P2 + 8*k2(vb)*q2L*q1tq2t*
     & q2tk1t*P1 + 8*k2(vb)*q2L*q1tq2t*q2tk1t*P2 + 4*k2(vb)*q2L*
     & q1tq2t*P1*m1*m2 - 8*k2(vb)*q2L*q1tq2t*P1*q1k1 - 4*k2(vb)*q2L
     & *q1tq2t*P1*k1k2 - 4*k2(vb)*q2L*q1tq2t*P2*m1*m2 - 8*k2(vb)*
     & q2L*q1tq2t*P2*q1k1 + 4*k2(vb)*q2L*q1tq2t*P2*k1k2 + 4*k2(vb)*
     & q2L*q1tk1t*t2*P1 + 4*k2(vb)*q2L*q1tk1t*t2*P2 - 4*k2(vb)*q2L
     & *q2tk1t*t1*P1
      A58(3,vb) = A58(3,vb) - 4*k2(vb)*q2L*q2tk1t*t1*P2 + 8*k2(vb)*
     & q1tq2t*q1tL*P1*m1*m2 - 8*k2(vb)*q1tq2t*q1tL*P1*q1k1 - 8*
     & k2(vb)*q1tq2t*q1tL*P1*q2k1 - 8*k2(vb)*q1tq2t*q1tL*P1*k1k2 - 
     & 8*k2(vb)*q1tq2t*q1tL*P2*m1*m2 - 8*k2(vb)*q1tq2t*q1tL*P2*q1k1
     &  - 8*k2(vb)*q1tq2t*q1tL*P2*q2k1 + 8*k2(vb)*q1tq2t*q1tL*P2*
     & k1k2 - 8*k2(vb)*q1tq2t*q2tL*P1*m1*m2 + 8*k2(vb)*q1tq2t*q2tL*
     & P1*q1k1 + 8*k2(vb)*q1tq2t*q2tL*P1*q2k1 + 8*k2(vb)*q1tq2t*
     & q2tL*P1*k1k2 + 8*k2(vb)*q1tq2t*q2tL*P2*m1*m2 + 8*k2(vb)*
     & q1tq2t*q2tL*P2*q1k1 + 8*k2(vb)*q1tq2t*q2tL*P2*q2k1 - 8*
     & k2(vb)*q1tq2t*q2tL*P2*k1k2 - 4*k2(vb)*q1tL*t2*P1*m1*m2 + 4*
     & k2(vb)*q1tL*t2*P1*q1k1 + 4*k2(vb)*q1tL*t2*P1*q2k1 + 4*k2(vb)*
     & q1tL*t2*P1*k1k2 + 4*k2(vb)*q1tL*t2*P2*m1*m2 + 4*k2(vb)*q1tL*
     & t2*P2*q1k1 + 4*k2(vb)*q1tL*t2*P2*q2k1 - 4*k2(vb)*q1tL*t2*P2*
     & k1k2 + 4*k2(vb)*q2tL*t1*P1*m1*m2 - 4*k2(vb)*q2tL*t1*P1*q1k1 - 
     & 4*k2(vb)*q2tL*t1*P1*q2k1 - 4*k2(vb)*q2tL*t1*P1*k1k2 - 4*k2(vb)
     & *q2tL*t1*P2*m1*m2
      A58(3,vb) = A58(3,vb) - 4*k2(vb)*q2tL*t1*P2*q1k1 - 4*k2(vb)*
     & q2tL*t1*P2*q2k1 + 4*k2(vb)*q2tL*t1*P2*k1k2 - 4*kw(vb)*q1L*
     & q1tq2t*P1*m1*m2 + 4*kw(vb)*q1L*q1tq2t*P1*k1k2 + 4*kw(vb)*q1L
     & *q1tq2t*P2*m1*m2 - 4*kw(vb)*q1L*q1tq2t*P2*k1k2 + 4*kw(vb)*
     & q2L*q1tq2t*P1*m1*m2 - 4*kw(vb)*q2L*q1tq2t*P1*k1k2 - 4*kw(vb)
     & *q2L*q1tq2t*P2*m1*m2 + 4*kw(vb)*q2L*q1tq2t*P2*k1k2 + 8*
     & kw(vb)*q1tq2t*q1tL*P1*m1*m2 - 8*kw(vb)*q1tq2t*q1tL*P1*k1k2
     &  - 8*kw(vb)*q1tq2t*q1tL*P2*m1*m2 + 8*kw(vb)*q1tq2t*q1tL*P2*
     & k1k2 - 8*kw(vb)*q1tq2t*q2tL*P1*m1*m2 + 8*kw(vb)*q1tq2t*q2tL*
     & P1*k1k2 + 8*kw(vb)*q1tq2t*q2tL*P2*m1*m2 - 8*kw(vb)*q1tq2t*
     & q2tL*P2*k1k2 - 4*kw(vb)*q1tL*t2*P1*m1*m2 + 4*kw(vb)*q1tL*t2*
     & P1*k1k2 + 4*kw(vb)*q1tL*t2*P2*m1*m2 - 4*kw(vb)*q1tL*t2*P2*k1k2
     &  + 4*kw(vb)*q2tL*t1*P1*m1*m2 - 4*kw(vb)*q2tL*t1*P1*k1k2 - 4*
     & kw(vb)*q2tL*t1*P2*m1*m2 + 4*kw(vb)*q2tL*t1*P2*k1k2 + 4*q1(vb)*
     & q1L*q1tq2t*P1*m1*m2 - 4*q1(vb)*q1L*q1tq2t*P1*k1k2 - 4*q1(vb)
     & *q1L*q1tq2t*P2*m1*m2
      A58(3,vb) = A58(3,vb) + 4*q1(vb)*q1L*q1tq2t*P2*k1k2 + 4*q1(vb)*
     & q2L*q1tq2t*P1*m1*m2 - 4*q1(vb)*q2L*q1tq2t*P1*k1k2 - 4*q1(vb)
     & *q2L*q1tq2t*P2*m1*m2 + 4*q1(vb)*q2L*q1tq2t*P2*k1k2 - 4*
     & q2(vb)*q1L*q1tq2t*P1*m1*m2 + 4*q2(vb)*q1L*q1tq2t*P1*k1k2 + 4
     & *q2(vb)*q1L*q1tq2t*P2*m1*m2 - 4*q2(vb)*q1L*q1tq2t*P2*k1k2 - 
     & 4*q2(vb)*q2L*q1tq2t*P1*m1*m2 + 4*q2(vb)*q2L*q1tq2t*P1*k1k2
     &  + 4*q2(vb)*q2L*q1tq2t*P2*m1*m2 - 4*q2(vb)*q2L*q1tq2t*P2*
     & k1k2 - 8*q1t(vb)*q1L*q1tq2t*P1*m1*m2 + 8*q1t(vb)*q1L*q1tq2t*
     & P1*k1k2 + 8*q1t(vb)*q1L*q1tq2t*P2*m1*m2 - 8*q1t(vb)*q1L*
     & q1tq2t*P2*k1k2 + 4*q1t(vb)*q1L*t2*P1*m1*m2 - 4*q1t(vb)*q1L*t2
     & *P1*k1k2 - 4*q1t(vb)*q1L*t2*P2*m1*m2 + 4*q1t(vb)*q1L*t2*P2*
     & k1k2 - 8*q1t(vb)*q2L*q1tq2t*P1*m1*m2 + 8*q1t(vb)*q2L*q1tq2t*
     & P1*k1k2 + 8*q1t(vb)*q2L*q1tq2t*P2*m1*m2 - 8*q1t(vb)*q2L*
     & q1tq2t*P2*k1k2 + 4*q1t(vb)*q2L*t2*P1*m1*m2 - 4*q1t(vb)*q2L*t2
     & *P1*k1k2 - 4*q1t(vb)*q2L*t2*P2*m1*m2 + 4*q1t(vb)*q2L*t2*P2*
     & k1k2
      A58(3,vb) = A58(3,vb) + 8*q2t(vb)*q1L*q1tq2t*P1*m1*m2 - 8*
     & q2t(vb)*q1L*q1tq2t*P1*k1k2 - 8*q2t(vb)*q1L*q1tq2t*P2*m1*m2
     &  + 8*q2t(vb)*q1L*q1tq2t*P2*k1k2 - 4*q2t(vb)*q1L*t1*P1*m1*m2
     &  + 4*q2t(vb)*q1L*t1*P1*k1k2 + 4*q2t(vb)*q1L*t1*P2*m1*m2 - 4*
     & q2t(vb)*q1L*t1*P2*k1k2 + 8*q2t(vb)*q2L*q1tq2t*P1*m1*m2 - 8*
     & q2t(vb)*q2L*q1tq2t*P1*k1k2 - 8*q2t(vb)*q2L*q1tq2t*P2*m1*m2
     &  + 8*q2t(vb)*q2L*q1tq2t*P2*k1k2 - 4*q2t(vb)*q2L*t1*P1*m1*m2
     &  + 4*q2t(vb)*q2L*t1*P1*k1k2 + 4*q2t(vb)*q2L*t1*P2*m1*m2 - 4*
     & q2t(vb)*q2L*t1*P2*k1k2 + 16*L(vb)*q1tq2t*q1tk1t*P1*m1*m2 + 8*
     & L(vb)*q1tq2t*q1tk1t*P1*q1k2 + 8*L(vb)*q1tq2t*q1tk1t*P1*q2k2
     &  - 16*L(vb)*q1tq2t*q1tk1t*P1*k1k2 + 8*L(vb)*q1tq2t*q1tk1t*P2
     & *q1k2 + 8*L(vb)*q1tq2t*q1tk1t*P2*q2k2 - 8*L(vb)*q1tq2t*
     & q1tk2t*P1*q1k1 - 8*L(vb)*q1tq2t*q1tk2t*P1*q2k1 - 16*L(vb)*
     & q1tq2t*q1tk2t*P2*m1*m2 - 8*L(vb)*q1tq2t*q1tk2t*P2*q1k1 - 8*
     & L(vb)*q1tq2t*q1tk2t*P2*q2k1 + 16*L(vb)*q1tq2t*q1tk2t*P2*k1k2
     &  - 16*L(vb)*q1tq2t*q2tk1t*P1*m1*m2
      A58(3,vb) = A58(3,vb) - 8*L(vb)*q1tq2t*q2tk1t*P1*q1k2 - 8*L(vb)
     & *q1tq2t*q2tk1t*P1*q2k2 + 16*L(vb)*q1tq2t*q2tk1t*P1*k1k2 - 8*
     & L(vb)*q1tq2t*q2tk1t*P2*q1k2 - 8*L(vb)*q1tq2t*q2tk1t*P2*q2k2
     &  + 8*L(vb)*q1tq2t*q2tk2t*P1*q1k1 + 8*L(vb)*q1tq2t*q2tk2t*P1*
     & q2k1 + 16*L(vb)*q1tq2t*q2tk2t*P2*m1*m2 + 8*L(vb)*q1tq2t*
     & q2tk2t*P2*q1k1 + 8*L(vb)*q1tq2t*q2tk2t*P2*q2k1 - 16*L(vb)*
     & q1tq2t*q2tk2t*P2*k1k2 + 8*L(vb)*q1tq2t*P1*m1*m2*q1k1 - 8*
     & L(vb)*q1tq2t*P1*m1*m2*q2k1 + 8*L(vb)*q1tq2t*P1*q1k1*q2k2 - 8*
     & L(vb)*q1tq2t*P1*q1k1*k1k2 - 8*L(vb)*q1tq2t*P1*q1k2*q2k1 + 8*
     & L(vb)*q1tq2t*P1*q2k1*k1k2 - 8*L(vb)*q1tq2t*P2*m1*m2*q1k2 + 8*
     & L(vb)*q1tq2t*P2*m1*m2*q2k2 + 8*L(vb)*q1tq2t*P2*q1k1*q2k2 - 8*
     & L(vb)*q1tq2t*P2*q1k2*q2k1 + 8*L(vb)*q1tq2t*P2*q1k2*k1k2 - 8*
     & L(vb)*q1tq2t*P2*q2k2*k1k2 - 8*L(vb)*q1tk1t*t2*P1*m1*m2 - 4*
     & L(vb)*q1tk1t*t2*P1*q1k2 - 4*L(vb)*q1tk1t*t2*P1*q2k2 + 8*L(vb)*
     & q1tk1t*t2*P1*k1k2 - 4*L(vb)*q1tk1t*t2*P2*q1k2 - 4*L(vb)*
     & q1tk1t*t2*P2*q2k2
      A58(3,vb) = A58(3,vb) + 4*L(vb)*q1tk2t*t2*P1*q1k1 + 4*L(vb)*
     & q1tk2t*t2*P1*q2k1 + 8*L(vb)*q1tk2t*t2*P2*m1*m2 + 4*L(vb)*
     & q1tk2t*t2*P2*q1k1 + 4*L(vb)*q1tk2t*t2*P2*q2k1 - 8*L(vb)*
     & q1tk2t*t2*P2*k1k2 + 8*L(vb)*q2tk1t*t1*P1*m1*m2 + 4*L(vb)*
     & q2tk1t*t1*P1*q1k2 + 4*L(vb)*q2tk1t*t1*P1*q2k2 - 8*L(vb)*
     & q2tk1t*t1*P1*k1k2 + 4*L(vb)*q2tk1t*t1*P2*q1k2 + 4*L(vb)*
     & q2tk1t*t1*P2*q2k2 - 4*L(vb)*q2tk2t*t1*P1*q1k1 - 4*L(vb)*
     & q2tk2t*t1*P1*q2k1 - 8*L(vb)*q2tk2t*t1*P2*m1*m2 - 4*L(vb)*
     & q2tk2t*t1*P2*q1k1 - 4*L(vb)*q2tk2t*t1*P2*q2k1 + 8*L(vb)*
     & q2tk2t*t1*P2*k1k2

      A58(4,vb) =
     & 8*k1(vb)*q1L*q2U*q1tq2t*P1*m2 + 8*k1(vb)*q1L*q2U*q1tq2t*P2
     & *m2 + 8*k1(vb)*q1L*q1tq2t*q1tU*P1*m2 + 8*k1(vb)*q1L*q1tq2t*
     & q1tU*P2*m2 - 8*k1(vb)*q1L*q1tq2t*q2tU*P1*m2 - 8*k1(vb)*q1L*
     & q1tq2t*q2tU*P2*m2 - 4*k1(vb)*q1L*q1tU*t2*P1*m2 - 4*k1(vb)*
     & q1L*q1tU*t2*P2*m2 + 4*k1(vb)*q1L*q2tU*t1*P1*m2 + 4*k1(vb)*
     & q1L*q2tU*t1*P2*m2 + 8*k1(vb)*q2U*q2L*q1tq2t*P1*m2 + 8*
     & k1(vb)*q2U*q2L*q1tq2t*P2*m2 + 8*k1(vb)*q2L*q1tq2t*q1tU*P1*
     & m2 + 8*k1(vb)*q2L*q1tq2t*q1tU*P2*m2 - 8*k1(vb)*q2L*q1tq2t*
     & q2tU*P1*m2 - 8*k1(vb)*q2L*q1tq2t*q2tU*P2*m2 - 4*k1(vb)*q2L*
     & q1tU*t2*P1*m2 - 4*k1(vb)*q2L*q1tU*t2*P2*m2 + 4*k1(vb)*q2L*
     & q2tU*t1*P1*m2 + 4*k1(vb)*q2L*q2tU*t1*P2*m2 + 8*k2(vb)*q1L*
     & q2U*q1tq2t*P1*m1 + 8*k2(vb)*q1L*q2U*q1tq2t*P2*m1 + 8*k2(vb)
     & *q1L*q1tq2t*q1tU*P1*m1 + 8*k2(vb)*q1L*q1tq2t*q1tU*P2*m1 - 
     & 8*k2(vb)*q1L*q1tq2t*q2tU*P1*m1 - 8*k2(vb)*q1L*q1tq2t*q2tU*
     & P2*m1
      A58(4,vb) = A58(4,vb) - 4*k2(vb)*q1L*q1tU*t2*P1*m1 - 4*k2(vb)*
     & q1L*q1tU*t2*P2*m1 + 4*k2(vb)*q1L*q2tU*t1*P1*m1 + 4*k2(vb)*
     & q1L*q2tU*t1*P2*m1 + 8*k2(vb)*q2U*q2L*q1tq2t*P1*m1 + 8*
     & k2(vb)*q2U*q2L*q1tq2t*P2*m1 + 8*k2(vb)*q2L*q1tq2t*q1tU*P1*
     & m1 + 8*k2(vb)*q2L*q1tq2t*q1tU*P2*m1 - 8*k2(vb)*q2L*q1tq2t*
     & q2tU*P1*m1 - 8*k2(vb)*q2L*q1tq2t*q2tU*P2*m1 - 4*k2(vb)*q2L*
     & q1tU*t2*P1*m1 - 4*k2(vb)*q2L*q1tU*t2*P2*m1 + 4*k2(vb)*q2L*
     & q2tU*t1*P1*m1 + 4*k2(vb)*q2L*q2tU*t1*P2*m1 + 8*U(vb)*q1L*
     & q1tq2t*q1tk1t*P1*m2 + 8*U(vb)*q1L*q1tq2t*q1tk1t*P2*m2 + 8*
     & U(vb)*q1L*q1tq2t*q1tk2t*P1*m1 + 8*U(vb)*q1L*q1tq2t*q1tk2t*
     & P2*m1 - 8*U(vb)*q1L*q1tq2t*q2tk1t*P1*m2 - 8*U(vb)*q1L*
     & q1tq2t*q2tk1t*P2*m2 - 8*U(vb)*q1L*q1tq2t*q2tk2t*P1*m1 - 8*
     & U(vb)*q1L*q1tq2t*q2tk2t*P2*m1 - 8*U(vb)*q1L*q1tq2t*P1*m1*
     & q2k2 - 8*U(vb)*q1L*q1tq2t*P1*m2*q2k1 - 8*U(vb)*q1L*q1tq2t*P2
     & *m1*q2k2 - 8*U(vb)*q1L*q1tq2t*P2*m2*q2k1 - 4*U(vb)*q1L*
     & q1tk1t*t2*P1*m2
      A58(4,vb) = A58(4,vb) - 4*U(vb)*q1L*q1tk1t*t2*P2*m2 - 4*U(vb)*
     & q1L*q1tk2t*t2*P1*m1 - 4*U(vb)*q1L*q1tk2t*t2*P2*m1 + 4*U(vb)*
     & q1L*q2tk1t*t1*P1*m2 + 4*U(vb)*q1L*q2tk1t*t1*P2*m2 + 4*U(vb)*
     & q1L*q2tk2t*t1*P1*m1 + 4*U(vb)*q1L*q2tk2t*t1*P2*m1 + 8*U(vb)*
     & q2L*q1tq2t*q1tk1t*P1*m2 + 8*U(vb)*q2L*q1tq2t*q1tk1t*P2*m2
     &  + 8*U(vb)*q2L*q1tq2t*q1tk2t*P1*m1 + 8*U(vb)*q2L*q1tq2t*
     & q1tk2t*P2*m1 - 8*U(vb)*q2L*q1tq2t*q2tk1t*P1*m2 - 8*U(vb)*
     & q2L*q1tq2t*q2tk1t*P2*m2 - 8*U(vb)*q2L*q1tq2t*q2tk2t*P1*m1
     &  - 8*U(vb)*q2L*q1tq2t*q2tk2t*P2*m1 + 8*U(vb)*q2L*q1tq2t*P1*
     & m1*q1k2 + 8*U(vb)*q2L*q1tq2t*P1*m2*q1k1 + 8*U(vb)*q2L*q1tq2t
     & *P2*m1*q1k2 + 8*U(vb)*q2L*q1tq2t*P2*m2*q1k1 - 4*U(vb)*q2L*
     & q1tk1t*t2*P1*m2 - 4*U(vb)*q2L*q1tk1t*t2*P2*m2 - 4*U(vb)*q2L*
     & q1tk2t*t2*P1*m1 - 4*U(vb)*q2L*q1tk2t*t2*P2*m1 + 4*U(vb)*q2L*
     & q2tk1t*t1*P1*m2 + 4*U(vb)*q2L*q2tk1t*t1*P2*m2 + 4*U(vb)*q2L*
     & q2tk2t*t1*P1*m1 + 4*U(vb)*q2L*q2tk2t*t1*P2*m1 + 8*U(vb)*
     & q1tq2t*q1tL*P1*m1*q1k2
      A58(4,vb) = A58(4,vb) + 8*U(vb)*q1tq2t*q1tL*P1*m1*q2k2 + 8*
     & U(vb)*q1tq2t*q1tL*P1*m2*q1k1 + 8*U(vb)*q1tq2t*q1tL*P1*m2*
     & q2k1 + 8*U(vb)*q1tq2t*q1tL*P2*m1*q1k2 + 8*U(vb)*q1tq2t*q1tL*
     & P2*m1*q2k2 + 8*U(vb)*q1tq2t*q1tL*P2*m2*q1k1 + 8*U(vb)*q1tq2t*
     & q1tL*P2*m2*q2k1 - 8*U(vb)*q1tq2t*q2tL*P1*m1*q1k2 - 8*U(vb)*
     & q1tq2t*q2tL*P1*m1*q2k2 - 8*U(vb)*q1tq2t*q2tL*P1*m2*q1k1 - 8*
     & U(vb)*q1tq2t*q2tL*P1*m2*q2k1 - 8*U(vb)*q1tq2t*q2tL*P2*m1*
     & q1k2 - 8*U(vb)*q1tq2t*q2tL*P2*m1*q2k2 - 8*U(vb)*q1tq2t*q2tL*
     & P2*m2*q1k1 - 8*U(vb)*q1tq2t*q2tL*P2*m2*q2k1 - 4*U(vb)*q1tL*t2
     & *P1*m1*q1k2 - 4*U(vb)*q1tL*t2*P1*m1*q2k2 - 4*U(vb)*q1tL*t2*P1*
     & m2*q1k1 - 4*U(vb)*q1tL*t2*P1*m2*q2k1 - 4*U(vb)*q1tL*t2*P2*m1*
     & q1k2 - 4*U(vb)*q1tL*t2*P2*m1*q2k2 - 4*U(vb)*q1tL*t2*P2*m2*q1k1
     &  - 4*U(vb)*q1tL*t2*P2*m2*q2k1 + 4*U(vb)*q2tL*t1*P1*m1*q1k2 + 4
     & *U(vb)*q2tL*t1*P1*m1*q2k2 + 4*U(vb)*q2tL*t1*P1*m2*q1k1 + 4*
     & U(vb)*q2tL*t1*P1*m2*q2k1 + 4*U(vb)*q2tL*t1*P2*m1*q1k2 + 4*
     & U(vb)*q2tL*t1*P2*m1*q2k2
      A58(4,vb) = A58(4,vb) + 4*U(vb)*q2tL*t1*P2*m2*q1k1 + 4*U(vb)*
     & q2tL*t1*P2*m2*q2k1 - 8*L(vb)*q2U*q1tq2t*P1*m1*q1k2 - 8*L(vb)*
     & q2U*q1tq2t*P1*m1*q2k2 - 8*L(vb)*q2U*q1tq2t*P1*m2*q1k1 - 8*
     & L(vb)*q2U*q1tq2t*P1*m2*q2k1 - 8*L(vb)*q2U*q1tq2t*P2*m1*q1k2
     &  - 8*L(vb)*q2U*q1tq2t*P2*m1*q2k2 - 8*L(vb)*q2U*q1tq2t*P2*m2*
     & q1k1 - 8*L(vb)*q2U*q1tq2t*P2*m2*q2k1 - 8*L(vb)*q1tq2t*q1tU*
     & P1*m1*q1k2 - 8*L(vb)*q1tq2t*q1tU*P1*m1*q2k2 - 8*L(vb)*q1tq2t*
     & q1tU*P1*m2*q1k1 - 8*L(vb)*q1tq2t*q1tU*P1*m2*q2k1 - 8*L(vb)*
     & q1tq2t*q1tU*P2*m1*q1k2 - 8*L(vb)*q1tq2t*q1tU*P2*m1*q2k2 - 8*
     & L(vb)*q1tq2t*q1tU*P2*m2*q1k1 - 8*L(vb)*q1tq2t*q1tU*P2*m2*
     & q2k1 + 8*L(vb)*q1tq2t*q2tU*P1*m1*q1k2 + 8*L(vb)*q1tq2t*q2tU*
     & P1*m1*q2k2 + 8*L(vb)*q1tq2t*q2tU*P1*m2*q1k1 + 8*L(vb)*q1tq2t*
     & q2tU*P1*m2*q2k1 + 8*L(vb)*q1tq2t*q2tU*P2*m1*q1k2 + 8*L(vb)*
     & q1tq2t*q2tU*P2*m1*q2k2 + 8*L(vb)*q1tq2t*q2tU*P2*m2*q1k1 + 8*
     & L(vb)*q1tq2t*q2tU*P2*m2*q2k1 + 4*L(vb)*q1tU*t2*P1*m1*q1k2 + 4
     & *L(vb)*q1tU*t2*P1*m1*q2k2
      A58(4,vb) = A58(4,vb) + 4*L(vb)*q1tU*t2*P1*m2*q1k1 + 4*L(vb)*
     & q1tU*t2*P1*m2*q2k1 + 4*L(vb)*q1tU*t2*P2*m1*q1k2 + 4*L(vb)*
     & q1tU*t2*P2*m1*q2k2 + 4*L(vb)*q1tU*t2*P2*m2*q1k1 + 4*L(vb)*
     & q1tU*t2*P2*m2*q2k1 - 4*L(vb)*q2tU*t1*P1*m1*q1k2 - 4*L(vb)*
     & q2tU*t1*P1*m1*q2k2 - 4*L(vb)*q2tU*t1*P1*m2*q1k1 - 4*L(vb)*
     & q2tU*t1*P1*m2*q2k1 - 4*L(vb)*q2tU*t1*P2*m1*q1k2 - 4*L(vb)*
     & q2tU*t1*P2*m1*q2k2 - 4*L(vb)*q2tU*t1*P2*m2*q1k1 - 4*L(vb)*
     & q2tU*t1*P2*m2*q2k1

      A51(1,vb) =
     & 4*k1(vb)*t1*t2*P3*m2 + 4*k1(vb)*t1*t2*P4*m2 - 4*k1(vb)*t1*t2*P7*
     & m2 - 12*k1(vb)*q1tq2t*q1tk1t*P3*m2 + 4*k1(vb)*q1tq2t*q1tk2t*
     & P3*m1 + 4*k1(vb)*q1tq2t*q2tk1t*P4*m2 + 4*k1(vb)*q1tq2t*
     & q2tk2t*P4*m1 - 16*k1(vb)*q1tq2t*q2tk2t*P4*m2 - 4*k1(vb)*
     & q1tq2t*t1*P3*m2 - 4*k1(vb)*q1tq2t*t2*P4*m2 + 4*k1(vb)*q1tq2t*
     & P3*m1*q1k2 + 4*k1(vb)*q1tq2t*P3*m2*q1q2 - 4*k1(vb)*q1tq2t*P3*
     & m2*q1k1 + 4*k1(vb)*q1tq2t*P4*m1*q2k2 + 4*k1(vb)*q1tq2t*P4*m2*
     & q1q2 + 4*k1(vb)*q1tq2t*P4*m2*q2k1 - 8*k1(vb)*q1tq2t*P4*m2*q2k2
     &  + 4*k1(vb)*q1tq2t*P7*m2*q1q2 + 4*k1(vb)*q1tq2t**2*P3*m2 + 4*
     & k1(vb)*q1tq2t**2*P4*m2 + 4*k1(vb)*q1tq2t**2*P7*m2 + 8*k1(vb)*
     & q1tk1t*q2tk1t*P3*m2 - 8*k1(vb)*q1tk1t*q2tk2t*P3*m1 - 8*
     & k1(vb)*q1tk1t*q2tk2t*P4*m2 - 16*k1(vb)*q1tk1t*q2tk2t*P7*m2
     &  + 8*k1(vb)*q1tk1t*t2*P3*m2 - 4*k1(vb)*q1tk1t*t2*P4*m2 - 8*
     & k1(vb)*q1tk1t*t2*P7*m2 - 8*k1(vb)*q1tk2t*q2tk2t*P4*m1 + 16*
     & k1(vb)*q1tk2t*q2tk2t*P4*m2
      A51(1,vb) = A51(1,vb) - 4*k1(vb)*q1tk2t*t2*P4*m1 + 8*k1(vb)*
     & q1tk2t*t2*P4*m2 + 4*k1(vb)*q2tk1t*t1*P3*m2 - 4*k1(vb)*q2tk2t*
     & t1*P3*m1 + 8*k1(vb)*q2tk2t*t1*P4*m2 - 8*k1(vb)*q2tk2t*t1*P7*m2
     &  + 4*k2(vb)*t1*t2*P3*m1 + 4*k2(vb)*t1*t2*P4*m1 - 4*k2(vb)*t1*t2*
     & P7*m1 - 16*k2(vb)*q1tq2t*q1tk1t*P3*m1 + 4*k2(vb)*q1tq2t*
     & q1tk1t*P3*m2 + 4*k2(vb)*q1tq2t*q1tk2t*P3*m1 + 4*k2(vb)*
     & q1tq2t*q2tk1t*P4*m2 - 12*k2(vb)*q1tq2t*q2tk2t*P4*m1 - 4*
     & k2(vb)*q1tq2t*t1*P3*m1 - 4*k2(vb)*q1tq2t*t2*P4*m1 + 4*k2(vb)*
     & q1tq2t*P3*m1*q1q2 - 8*k2(vb)*q1tq2t*P3*m1*q1k1 + 4*k2(vb)*
     & q1tq2t*P3*m1*q1k2 + 4*k2(vb)*q1tq2t*P3*m2*q1k1 + 4*k2(vb)*
     & q1tq2t*P4*m1*q1q2 - 4*k2(vb)*q1tq2t*P4*m1*q2k2 + 4*k2(vb)*
     & q1tq2t*P4*m2*q2k1 + 4*k2(vb)*q1tq2t*P7*m1*q1q2 + 4*k2(vb)*
     & q1tq2t**2*P3*m1 + 4*k2(vb)*q1tq2t**2*P4*m1 + 4*k2(vb)*
     & q1tq2t**2*P7*m1 + 16*k2(vb)*q1tk1t*q2tk1t*P3*m1 - 8*k2(vb)*
     & q1tk1t*q2tk1t*P3*m2 - 8*k2(vb)*q1tk1t*q2tk2t*P3*m1 - 8*
     & k2(vb)*q1tk1t*q2tk2t*P4*m2
      A51(1,vb) = A51(1,vb) - 16*k2(vb)*q1tk1t*q2tk2t*P7*m1 + 8*
     & k2(vb)*q1tk1t*t2*P3*m1 - 4*k2(vb)*q1tk1t*t2*P4*m2 - 8*k2(vb)*
     & q1tk1t*t2*P7*m1 + 8*k2(vb)*q1tk2t*q2tk2t*P4*m1 + 4*k2(vb)*
     & q1tk2t*t2*P4*m1 + 8*k2(vb)*q2tk1t*t1*P3*m1 - 4*k2(vb)*q2tk1t*
     & t1*P3*m2 - 4*k2(vb)*q2tk2t*t1*P3*m1 + 8*k2(vb)*q2tk2t*t1*P4*m1
     &  - 8*k2(vb)*q2tk2t*t1*P7*m1 + 4*kw(vb)*q1tq2t*q1tk1t*P3*m2 + 
     & 4*kw(vb)*q1tq2t*q1tk2t*P3*m1 + 4*kw(vb)*q1tq2t*q2tk1t*P4*m2
     &  + 4*kw(vb)*q1tq2t*q2tk2t*P4*m1 + 4*kw(vb)*q1tq2t*P3*m1*q1k2
     &  + 4*kw(vb)*q1tq2t*P3*m2*q1k1 + 4*kw(vb)*q1tq2t*P4*m1*q2k2 + 4
     & *kw(vb)*q1tq2t*P4*m2*q2k1 - 8*kw(vb)*q1tk1t*q2tk1t*P3*m2 - 8*
     & kw(vb)*q1tk1t*q2tk2t*P3*m1 - 8*kw(vb)*q1tk1t*q2tk2t*P4*m2 - 
     & 4*kw(vb)*q1tk1t*t2*P4*m2 - 8*kw(vb)*q1tk2t*q2tk2t*P4*m1 - 4*
     & kw(vb)*q1tk2t*t2*P4*m1 - 4*kw(vb)*q2tk1t*t1*P3*m2 - 4*kw(vb)*
     & q2tk2t*t1*P3*m1 + 4*q1(vb)*q1tq2t*q1tk1t*P3*m2 + 4*q1(vb)*
     & q1tq2t*q1tk2t*P3*m1 - 4*q1(vb)*q1tq2t*q2tk1t*P3*m2 - 4*
     & q1(vb)*q1tq2t*q2tk1t*P7*m2
      A51(1,vb) = A51(1,vb) - 4*q1(vb)*q1tq2t*q2tk2t*P3*m1 - 4*
     & q1(vb)*q1tq2t*q2tk2t*P7*m1 - 4*q1(vb)*q1tq2t*P3*m1*q1k2 - 4*
     & q1(vb)*q1tq2t*P3*m1*q2k2 - 4*q1(vb)*q1tq2t*P3*m2*q1k1 - 4*
     & q1(vb)*q1tq2t*P3*m2*q2k1 - 4*q1(vb)*q1tq2t*P7*m1*q2k2 - 4*
     & q1(vb)*q1tq2t*P7*m2*q2k1 + 8*q1(vb)*q1tk1t*q2tk2t*P3*m1 + 8*
     & q1(vb)*q1tk1t*q2tk2t*P7*m2 - 4*q1(vb)*q1tk1t*t2*P3*m2 + 4*
     & q1(vb)*q1tk1t*t2*P7*m2 - 8*q1(vb)*q1tk2t*q2tk1t*P3*m1 + 8*
     & q1(vb)*q1tk2t*q2tk2t*P7*m1 - 4*q1(vb)*q1tk2t*t2*P3*m1 + 4*
     & q1(vb)*q1tk2t*t2*P7*m1 + 4*q1(vb)*q2tk1t*t1*P3*m2 + 4*q1(vb)*
     & q2tk2t*t1*P3*m1 - 4*q2(vb)*q1tq2t*q1tk1t*P4*m2 - 4*q2(vb)*
     & q1tq2t*q1tk1t*P7*m2 - 4*q2(vb)*q1tq2t*q1tk2t*P4*m1 - 4*
     & q2(vb)*q1tq2t*q1tk2t*P7*m1 + 4*q2(vb)*q1tq2t*q2tk1t*P4*m2 + 
     & 4*q2(vb)*q1tq2t*q2tk2t*P4*m1 - 4*q2(vb)*q1tq2t*P4*m1*q1k2 - 4
     & *q2(vb)*q1tq2t*P4*m1*q2k2 - 4*q2(vb)*q1tq2t*P4*m2*q1k1 - 4*
     & q2(vb)*q1tq2t*P4*m2*q2k1 - 4*q2(vb)*q1tq2t*P7*m1*q1k2 - 4*
     & q2(vb)*q1tq2t*P7*m2*q1k1
      A51(1,vb) = A51(1,vb) + 8*q2(vb)*q1tk1t*q2tk1t*P7*m2 + 8*
     & q2(vb)*q1tk1t*q2tk2t*P4*m2 + 8*q2(vb)*q1tk1t*q2tk2t*P7*m1 + 
     & 4*q2(vb)*q1tk1t*t2*P4*m2 - 8*q2(vb)*q1tk2t*q2tk1t*P4*m2 + 4*
     & q2(vb)*q1tk2t*t2*P4*m1 - 4*q2(vb)*q2tk1t*t1*P4*m2 + 4*q2(vb)*
     & q2tk1t*t1*P7*m2 - 4*q2(vb)*q2tk2t*t1*P4*m1 + 4*q2(vb)*q2tk2t*
     & t1*P7*m1 - 4*q1t(vb)*t2*P3*m1*q1k2 - 4*q1t(vb)*t2*P3*m2*q1k1 - 4
     & *q1t(vb)*t2*P4*m1*q1k2 - 4*q1t(vb)*t2*P4*m2*q1k1 + 4*q1t(vb)*t2*
     & P7*m1*q1k2 + 4*q1t(vb)*t2*P7*m2*q1k1 + 8*q1t(vb)*q1tq2t*P3*m1*
     & q1k2 + 4*q1t(vb)*q1tq2t*P3*m1*q2k2 + 8*q1t(vb)*q1tq2t*P3*m2*
     & q1k1 + 4*q1t(vb)*q1tq2t*P3*m2*q2k1 - 4*q1t(vb)*q1tq2t*P4*m1*
     & q2k2 - 4*q1t(vb)*q1tq2t*P4*m2*q2k1 + 4*q1t(vb)*q1tq2t*P7*m1*
     & q2k2 + 4*q1t(vb)*q1tq2t*P7*m2*q2k1 - 4*q1t(vb)*q2tk1t*t1*P3*m2
     &  + 4*q1t(vb)*q2tk1t*t2*P4*m2 - 8*q1t(vb)*q2tk1t*P3*m1*q1k2 - 4
     & *q1t(vb)*q2tk1t*P3*m2*q1q2 + 4*q1t(vb)*q2tk1t*P4*m2*q1q2 - 8*
     & q1t(vb)*q2tk1t*P4*m2*q2k2 - 4*q1t(vb)*q2tk1t*P7*m2*q1q2 - 4*
     & q1t(vb)*q2tk2t*t1*P3*m1
      A51(1,vb) = A51(1,vb) + 4*q1t(vb)*q2tk2t*t2*P4*m1 - 4*q1t(vb)*
     & q2tk2t*P3*m1*q1q2 + 8*q1t(vb)*q2tk2t*P3*m1*q1k1 + 4*q1t(vb)*
     & q2tk2t*P4*m1*q1q2 - 8*q1t(vb)*q2tk2t*P4*m1*q1k2 - 8*q1t(vb)*
     & q2tk2t*P4*m1*q2k2 - 8*q1t(vb)*q2tk2t*P4*m2*q1k1 - 4*q1t(vb)*
     & q2tk2t*P7*m1*q1q2 + 8*q1t(vb)*q2tk2t*P7*m1*q1k2 + 8*q1t(vb)*
     & q2tk2t*P7*m2*q1k1 - 4*q2t(vb)*t1*P3*m1*q2k2 - 4*q2t(vb)*t1*P3*
     & m2*q2k1 - 4*q2t(vb)*t1*P4*m1*q2k2 - 4*q2t(vb)*t1*P4*m2*q2k1 + 4*
     & q2t(vb)*t1*P7*m1*q2k2 + 4*q2t(vb)*t1*P7*m2*q2k1 - 4*q2t(vb)*
     & q1tq2t*P3*m1*q1k2 - 4*q2t(vb)*q1tq2t*P3*m2*q1k1 + 4*q2t(vb)*
     & q1tq2t*P4*m1*q1k2 + 8*q2t(vb)*q1tq2t*P4*m1*q2k2 + 4*q2t(vb)*
     & q1tq2t*P4*m2*q1k1 + 8*q2t(vb)*q1tq2t*P4*m2*q2k1 + 4*q2t(vb)*
     & q1tq2t*P7*m1*q1k2 + 4*q2t(vb)*q1tq2t*P7*m2*q1k1 + 4*q2t(vb)*
     & q1tk1t*t1*P3*m2 - 4*q2t(vb)*q1tk1t*t2*P4*m2 - 8*q2t(vb)*
     & q1tk1t*P3*m1*q2k2 + 4*q2t(vb)*q1tk1t*P3*m2*q1q2 - 8*q2t(vb)*
     & q1tk1t*P3*m2*q1k1 - 8*q2t(vb)*q1tk1t*P3*m2*q2k1 - 4*q2t(vb)*
     & q1tk1t*P4*m2*q1q2
      A51(1,vb) = A51(1,vb) + 8*q2t(vb)*q1tk1t*P4*m2*q2k2 + 8*
     & q2t(vb)*q1tk1t*P7*m1*q2k2 - 4*q2t(vb)*q1tk1t*P7*m2*q1q2 + 8*
     & q2t(vb)*q1tk1t*P7*m2*q2k1 + 4*q2t(vb)*q1tk2t*t1*P3*m1 - 4*
     & q2t(vb)*q1tk2t*t2*P4*m1 + 4*q2t(vb)*q1tk2t*P3*m1*q1q2 - 8*
     & q2t(vb)*q1tk2t*P3*m1*q1k1 - 4*q2t(vb)*q1tk2t*P4*m1*q1q2 - 8*
     & q2t(vb)*q1tk2t*P4*m2*q2k1 - 4*q2t(vb)*q1tk2t*P7*m1*q1q2

      A51(2,vb) =
     &  - 4*k1(vb)*q2U*q1tq2t*q1tk2t*P3 + 4*k1(vb)*q2U*q1tq2t*
     & q1tk2t*P4 + 4*k1(vb)*q2U*q1tq2t*q1tk2t*P7 + 4*k1(vb)*q2U*
     & q1tq2t*q2tk2t*P3 - 4*k1(vb)*q2U*q1tq2t*q2tk2t*P4 + 4*k1(vb)
     & *q2U*q1tq2t*q2tk2t*P7 + 4*k1(vb)*q2U*q1tq2t*P3*m1*m2 + 4*
     & k1(vb)*q2U*q1tq2t*P3*q1k2 + 4*k1(vb)*q2U*q1tq2t*P3*q2k2 - 4*
     & k1(vb)*q2U*q1tq2t*P3*k1k2 - 4*k1(vb)*q2U*q1tq2t*P4*m1*m2 + 4
     & *k1(vb)*q2U*q1tq2t*P4*q1k2 + 4*k1(vb)*q2U*q1tq2t*P4*q2k2 + 4
     & *k1(vb)*q2U*q1tq2t*P4*k1k2 + 4*k1(vb)*q2U*q1tq2t*P7*q1k2 + 4
     & *k1(vb)*q2U*q1tq2t*P7*q2k2 - 4*k1(vb)*q2U*q1tk1t*q2tk2t*P3
     &  - 4*k1(vb)*q2U*q1tk1t*q2tk2t*P4 - 8*k1(vb)*q2U*q1tk1t*
     & q2tk2t*P7 + 4*k1(vb)*q2U*q1tk2t*q2tk1t*P3 + 4*k1(vb)*q2U*
     & q1tk2t*q2tk1t*P4 - 8*k1(vb)*q2U*q1tk2t*q2tk2t*P7 + 4*k1(vb)
     & *q2U*q1tk2t*t2*P3 - 4*k1(vb)*q2U*q1tk2t*t2*P4 - 4*k1(vb)*
     & q2U*q1tk2t*t2*P7 - 4*k1(vb)*q2U*q2tk2t*t1*P3 + 4*k1(vb)*q2U
     & *q2tk2t*t1*P4
      A51(2,vb) = A51(2,vb) - 4*k1(vb)*q2U*q2tk2t*t1*P7 + 4*k1(vb)*
     & q1tq2t*q1tU*P3*m1*m2 + 8*k1(vb)*q1tq2t*q1tU*P3*q1k2 + 4*
     & k1(vb)*q1tq2t*q1tU*P3*q2k2 - 4*k1(vb)*q1tq2t*q1tU*P3*k1k2 + 
     & 4*k1(vb)*q1tq2t*q1tU*P4*q2k2 + 4*k1(vb)*q1tq2t*q1tU*P7*q2k2
     &  - 4*k1(vb)*q1tq2t*q2tU*P3*q1k2 + 4*k1(vb)*q1tq2t*q2tU*P4*m1
     & *m2 - 4*k1(vb)*q1tq2t*q2tU*P4*q1k2 - 8*k1(vb)*q1tq2t*q2tU*P4
     & *q2k2 - 4*k1(vb)*q1tq2t*q2tU*P4*k1k2 - 4*k1(vb)*q1tq2t*q2tU*
     & P7*q1k2 - 8*k1(vb)*q1tk1t*q2tU*P3*m1*m2 - 4*k1(vb)*q1tk1t*
     & q2tU*P3*q1k2 - 8*k1(vb)*q1tk1t*q2tU*P3*q2k2 + 8*k1(vb)*
     & q1tk1t*q2tU*P3*k1k2 - 4*k1(vb)*q1tk1t*q2tU*P4*q2k2 - 8*
     & k1(vb)*q1tk1t*q2tU*P7*q2k2 + 4*k1(vb)*q1tk2t*q2tU*t1*P3 + 4*
     & k1(vb)*q1tk2t*q2tU*t2*P4 + 4*k1(vb)*q1tk2t*q2tU*P3*q1q2 - 4*
     & k1(vb)*q1tk2t*q2tU*P3*q1k1 + 4*k1(vb)*q1tk2t*q2tU*P4*q1q2 +
     & 4*k1(vb)*q1tk2t*q2tU*P4*q2k1 + 4*k1(vb)*q1tk2t*q2tU*P7*q1q2
     &  - 4*k1(vb)*q1tU*q2tk1t*P3*q1k2 + 4*k1(vb)*q1tU*q2tk1t*P4*
     & q2k2
      A51(2,vb) = A51(2,vb) - 4*k1(vb)*q1tU*q2tk2t*t1*P3 - 4*k1(vb)
     & *q1tU*q2tk2t*t2*P4 - 4*k1(vb)*q1tU*q2tk2t*P3*q1q2 + 4*k1(vb)
     & *q1tU*q2tk2t*P3*q1k1 - 8*k1(vb)*q1tU*q2tk2t*P4*m1*m2 - 4*
     & k1(vb)*q1tU*q2tk2t*P4*q1q2 + 8*k1(vb)*q1tU*q2tk2t*P4*q1k2 - 
     & 4*k1(vb)*q1tU*q2tk2t*P4*q2k1 + 8*k1(vb)*q1tU*q2tk2t*P4*q2k2
     &  + 8*k1(vb)*q1tU*q2tk2t*P4*k1k2 - 4*k1(vb)*q1tU*q2tk2t*P7*
     & q1q2 + 8*k1(vb)*q1tU*q2tk2t*P7*q1k2 - 4*k1(vb)*q1tU*t2*P3*
     & q1k2 - 4*k1(vb)*q1tU*t2*P4*m1*m2 + 4*k1(vb)*q1tU*t2*P4*q1k2 + 
     & 4*k1(vb)*q1tU*t2*P4*k1k2 + 4*k1(vb)*q1tU*t2*P7*q1k2 - 4*k1(vb)
     & *q2tU*t1*P3*m1*m2 - 4*k1(vb)*q2tU*t1*P3*q2k2 + 4*k1(vb)*q2tU*
     & t1*P3*k1k2 + 4*k1(vb)*q2tU*t1*P4*q2k2 - 4*k1(vb)*q2tU*t1*P7*
     & q2k2 + 4*k2(vb)*q2U*q1tq2t*q1tk1t*P3 - 4*k2(vb)*q2U*q1tq2t*
     & q1tk1t*P4 - 4*k2(vb)*q2U*q1tq2t*q1tk1t*P7 - 4*k2(vb)*q2U*
     & q1tq2t*q2tk1t*P3 + 4*k2(vb)*q2U*q1tq2t*q2tk1t*P4 - 4*k2(vb)
     & *q2U*q1tq2t*q2tk1t*P7 + 4*k2(vb)*q2U*q1tq2t*P3*m1*m2 - 4*
     & k2(vb)*q2U*q1tq2t*P3*q1k1
      A51(2,vb) = A51(2,vb) - 4*k2(vb)*q2U*q1tq2t*P3*q2k1 - 4*
     & k2(vb)*q2U*q1tq2t*P3*k1k2 - 4*k2(vb)*q2U*q1tq2t*P4*m1*m2 - 4
     & *k2(vb)*q2U*q1tq2t*P4*q1k1 - 4*k2(vb)*q2U*q1tq2t*P4*q2k1 + 4
     & *k2(vb)*q2U*q1tq2t*P4*k1k2 - 4*k2(vb)*q2U*q1tq2t*P7*q1k1 - 4
     & *k2(vb)*q2U*q1tq2t*P7*q2k1 + 8*k2(vb)*q2U*q1tk1t*q2tk1t*P7
     &  + 4*k2(vb)*q2U*q1tk1t*q2tk2t*P3 + 4*k2(vb)*q2U*q1tk1t*
     & q2tk2t*P4 + 8*k2(vb)*q2U*q1tk1t*q2tk2t*P7 - 4*k2(vb)*q2U*
     & q1tk1t*t2*P3 + 4*k2(vb)*q2U*q1tk1t*t2*P4 + 4*k2(vb)*q2U*
     & q1tk1t*t2*P7 - 4*k2(vb)*q2U*q1tk2t*q2tk1t*P3 - 4*k2(vb)*q2U
     & *q1tk2t*q2tk1t*P4 + 4*k2(vb)*q2U*q2tk1t*t1*P3 - 4*k2(vb)*
     & q2U*q2tk1t*t1*P4 + 4*k2(vb)*q2U*q2tk1t*t1*P7 + 4*k2(vb)*
     & q1tq2t*q1tU*P3*m1*m2 - 8*k2(vb)*q1tq2t*q1tU*P3*q1k1 - 4*
     & k2(vb)*q1tq2t*q1tU*P3*q2k1 - 4*k2(vb)*q1tq2t*q1tU*P3*k1k2 - 
     & 4*k2(vb)*q1tq2t*q1tU*P4*q2k1 - 4*k2(vb)*q1tq2t*q1tU*P7*q2k1
     &  + 4*k2(vb)*q1tq2t*q2tU*P3*q1k1 + 4*k2(vb)*q1tq2t*q2tU*P4*m1
     & *m2
      A51(2,vb) = A51(2,vb) + 4*k2(vb)*q1tq2t*q2tU*P4*q1k1 + 8*
     & k2(vb)*q1tq2t*q2tU*P4*q2k1 - 4*k2(vb)*q1tq2t*q2tU*P4*k1k2 + 
     & 4*k2(vb)*q1tq2t*q2tU*P7*q1k1 - 4*k2(vb)*q1tk1t*q2tU*t1*P3 - 
     & 4*k2(vb)*q1tk1t*q2tU*t2*P4 - 8*k2(vb)*q1tk1t*q2tU*P3*m1*m2
     &  - 4*k2(vb)*q1tk1t*q2tU*P3*q1q2 + 8*k2(vb)*q1tk1t*q2tU*P3*
     & q1k1 - 4*k2(vb)*q1tk1t*q2tU*P3*q1k2 + 8*k2(vb)*q1tk1t*q2tU*
     & P3*q2k1 + 8*k2(vb)*q1tk1t*q2tU*P3*k1k2 - 4*k2(vb)*q1tk1t*
     & q2tU*P4*q1q2 + 4*k2(vb)*q1tk1t*q2tU*P4*q2k2 - 4*k2(vb)*
     & q1tk1t*q2tU*P7*q1q2 + 8*k2(vb)*q1tk1t*q2tU*P7*q2k1 + 4*
     & k2(vb)*q1tk2t*q2tU*P3*q1k1 - 4*k2(vb)*q1tk2t*q2tU*P4*q2k1 + 
     & 4*k2(vb)*q1tU*q2tk1t*t1*P3 + 4*k2(vb)*q1tU*q2tk1t*t2*P4 + 4*
     & k2(vb)*q1tU*q2tk1t*P3*q1q2 + 4*k2(vb)*q1tU*q2tk1t*P3*q1k2 + 
     & 4*k2(vb)*q1tU*q2tk1t*P4*q1q2 - 4*k2(vb)*q1tU*q2tk1t*P4*q2k2
     &  + 4*k2(vb)*q1tU*q2tk1t*P7*q1q2 - 4*k2(vb)*q1tU*q2tk2t*P3*
     & q1k1 - 8*k2(vb)*q1tU*q2tk2t*P4*m1*m2 - 8*k2(vb)*q1tU*q2tk2t*
     & P4*q1k1
      A51(2,vb) = A51(2,vb) - 4*k2(vb)*q1tU*q2tk2t*P4*q2k1 + 8*
     & k2(vb)*q1tU*q2tk2t*P4*k1k2 - 8*k2(vb)*q1tU*q2tk2t*P7*q1k1 + 
     & 4*k2(vb)*q1tU*t2*P3*q1k1 - 4*k2(vb)*q1tU*t2*P4*m1*m2 - 4*
     & k2(vb)*q1tU*t2*P4*q1k1 + 4*k2(vb)*q1tU*t2*P4*k1k2 - 4*k2(vb)*
     & q1tU*t2*P7*q1k1 - 4*k2(vb)*q2tU*t1*P3*m1*m2 + 4*k2(vb)*q2tU*
     & t1*P3*q2k1 + 4*k2(vb)*q2tU*t1*P3*k1k2 - 4*k2(vb)*q2tU*t1*P4*
     & q2k1 + 4*k2(vb)*q2tU*t1*P7*q2k1 + 4*kw(vb)*q2U*q1tq2t*P3*m1*
     & m2 - 4*kw(vb)*q2U*q1tq2t*P3*k1k2 - 4*kw(vb)*q2U*q1tq2t*P4*m1
     & *m2 + 4*kw(vb)*q2U*q1tq2t*P4*k1k2 + 4*kw(vb)*q2U*q1tk1t*
     & q2tk2t*P3 - 4*kw(vb)*q2U*q1tk1t*q2tk2t*P4 - 4*kw(vb)*q2U*
     & q1tk2t*q2tk1t*P3 + 4*kw(vb)*q2U*q1tk2t*q2tk1t*P4 + 4*kw(vb)
     & *q1tq2t*q1tU*P3*m1*m2 - 4*kw(vb)*q1tq2t*q1tU*P3*k1k2 + 4*
     & kw(vb)*q1tq2t*q2tU*P4*m1*m2 - 4*kw(vb)*q1tq2t*q2tU*P4*k1k2
     &  - 8*kw(vb)*q1tk1t*q2tU*P3*m1*m2 - 4*kw(vb)*q1tk1t*q2tU*P3*
     & q1k2 + 8*kw(vb)*q1tk1t*q2tU*P3*k1k2 - 4*kw(vb)*q1tk1t*q2tU*
     & P4*q2k2
      A51(2,vb) = A51(2,vb) + 4*kw(vb)*q1tk2t*q2tU*P3*q1k1 + 4*
     & kw(vb)*q1tk2t*q2tU*P4*q2k1 + 4*kw(vb)*q1tU*q2tk1t*P3*q1k2 + 
     & 4*kw(vb)*q1tU*q2tk1t*P4*q2k2 - 4*kw(vb)*q1tU*q2tk2t*P3*q1k1
     &  - 8*kw(vb)*q1tU*q2tk2t*P4*m1*m2 - 4*kw(vb)*q1tU*q2tk2t*P4*
     & q2k1 + 8*kw(vb)*q1tU*q2tk2t*P4*k1k2 - 4*kw(vb)*q1tU*t2*P4*m1*
     & m2 + 4*kw(vb)*q1tU*t2*P4*k1k2 - 4*kw(vb)*q2tU*t1*P3*m1*m2 + 4*
     & kw(vb)*q2tU*t1*P3*k1k2 + 4*q1(vb)*q2U*q1tq2t*P7*m1*m2 - 4*
     & q1(vb)*q2U*q1tq2t*P7*k1k2 + 4*q1(vb)*q2U*q1tk1t*q2tk2t*P7
     &  - 4*q1(vb)*q2U*q1tk2t*q2tk1t*P7 + 4*q1(vb)*q1tq2t*q1tU*P3*
     & m1*m2 - 4*q1(vb)*q1tq2t*q1tU*P3*k1k2 - 4*q1(vb)*q1tq2t*q2tU*
     & P3*m1*m2 + 4*q1(vb)*q1tq2t*q2tU*P3*k1k2 - 4*q1(vb)*q1tq2t*
     & q2tU*P7*m1*m2 + 4*q1(vb)*q1tq2t*q2tU*P7*k1k2 + 8*q1(vb)*
     & q1tk1t*q2tU*P3*m1*m2 + 4*q1(vb)*q1tk1t*q2tU*P3*q1k2 + 4*
     & q1(vb)*q1tk1t*q2tU*P3*q2k2 - 8*q1(vb)*q1tk1t*q2tU*P3*k1k2 + 
     & 4*q1(vb)*q1tk1t*q2tU*P7*q2k2 - 4*q1(vb)*q1tk2t*q2tU*P3*q1k1
     &  - 4*q1(vb)*q1tk2t*q2tU*P3*q2k1
      A51(2,vb) = A51(2,vb) - 4*q1(vb)*q1tk2t*q2tU*P7*q2k1 - 8*
     & q1(vb)*q1tU*q2tk1t*P3*m1*m2 - 4*q1(vb)*q1tU*q2tk1t*P3*q1k2
     &  - 4*q1(vb)*q1tU*q2tk1t*P3*q2k2 + 8*q1(vb)*q1tU*q2tk1t*P3*
     & k1k2 - 4*q1(vb)*q1tU*q2tk1t*P7*q2k2 + 4*q1(vb)*q1tU*q2tk2t*
     & P3*q1k1 + 4*q1(vb)*q1tU*q2tk2t*P3*q2k1 + 8*q1(vb)*q1tU*
     & q2tk2t*P7*m1*m2 + 4*q1(vb)*q1tU*q2tk2t*P7*q2k1 - 8*q1(vb)*
     & q1tU*q2tk2t*P7*k1k2 - 4*q1(vb)*q1tU*t2*P3*m1*m2 + 4*q1(vb)*
     & q1tU*t2*P3*k1k2 + 4*q1(vb)*q1tU*t2*P7*m1*m2 - 4*q1(vb)*q1tU*
     & t2*P7*k1k2 + 4*q1(vb)*q2tU*t1*P3*m1*m2 - 4*q1(vb)*q2tU*t1*P3*
     & k1k2 - 4*q2(vb)*q2U*q1tq2t*P7*m1*m2 + 4*q2(vb)*q2U*q1tq2t*P7
     & *k1k2 - 4*q2(vb)*q2U*q1tk1t*q2tk2t*P7 + 4*q2(vb)*q2U*q1tk2t
     & *q2tk1t*P7 - 4*q2(vb)*q1tq2t*q1tU*P4*m1*m2 + 4*q2(vb)*q1tq2t
     & *q1tU*P4*k1k2 - 4*q2(vb)*q1tq2t*q1tU*P7*m1*m2 + 4*q2(vb)*
     & q1tq2t*q1tU*P7*k1k2 + 4*q2(vb)*q1tq2t*q2tU*P4*m1*m2 - 4*
     & q2(vb)*q1tq2t*q2tU*P4*k1k2 + 4*q2(vb)*q1tk1t*q2tU*P4*q1k2 + 
     & 4*q2(vb)*q1tk1t*q2tU*P4*q2k2
      A51(2,vb) = A51(2,vb) + 8*q2(vb)*q1tk1t*q2tU*P7*m1*m2 + 4*
     & q2(vb)*q1tk1t*q2tU*P7*q1k2 - 8*q2(vb)*q1tk1t*q2tU*P7*k1k2 - 
     & 8*q2(vb)*q1tk2t*q2tU*P4*m1*m2 - 4*q2(vb)*q1tk2t*q2tU*P4*q1k1
     &  - 4*q2(vb)*q1tk2t*q2tU*P4*q2k1 + 8*q2(vb)*q1tk2t*q2tU*P4*
     & k1k2 - 4*q2(vb)*q1tk2t*q2tU*P7*q1k1 - 4*q2(vb)*q1tU*q2tk1t*
     & P4*q1k2 - 4*q2(vb)*q1tU*q2tk1t*P4*q2k2 - 4*q2(vb)*q1tU*
     & q2tk1t*P7*q1k2 + 8*q2(vb)*q1tU*q2tk2t*P4*m1*m2 + 4*q2(vb)*
     & q1tU*q2tk2t*P4*q1k1 + 4*q2(vb)*q1tU*q2tk2t*P4*q2k1 - 8*
     & q2(vb)*q1tU*q2tk2t*P4*k1k2 + 4*q2(vb)*q1tU*q2tk2t*P7*q1k1 + 
     & 4*q2(vb)*q1tU*t2*P4*m1*m2 - 4*q2(vb)*q1tU*t2*P4*k1k2 - 4*
     & q2(vb)*q2tU*t1*P4*m1*m2 + 4*q2(vb)*q2tU*t1*P4*k1k2 + 4*q2(vb)*
     & q2tU*t1*P7*m1*m2 - 4*q2(vb)*q2tU*t1*P7*k1k2 + 4*q1t(vb)*q2U*
     & q1tq2t*P3*m1*m2 - 4*q1t(vb)*q2U*q1tq2t*P3*k1k2 + 4*q1t(vb)*
     & q2U*q1tq2t*P4*m1*m2 - 4*q1t(vb)*q2U*q1tq2t*P4*k1k2 - 4*
     & q1t(vb)*q2U*q1tq2t*P7*m1*m2 + 4*q1t(vb)*q2U*q1tq2t*P7*k1k2
     &  - 8*q1t(vb)*q2U*q2tk1t*P3*m1*m2
      A51(2,vb) = A51(2,vb) - 4*q1t(vb)*q2U*q2tk1t*P3*q1k2 - 4*
     & q1t(vb)*q2U*q2tk1t*P3*q2k2 + 8*q1t(vb)*q2U*q2tk1t*P3*k1k2 + 
     & 4*q1t(vb)*q2U*q2tk1t*P4*q1k2 + 4*q1t(vb)*q2U*q2tk1t*P4*q2k2
     &  - 4*q1t(vb)*q2U*q2tk1t*P7*q1k2 - 4*q1t(vb)*q2U*q2tk1t*P7*
     & q2k2 + 4*q1t(vb)*q2U*q2tk2t*P3*q1k1 + 4*q1t(vb)*q2U*q2tk2t*
     & P3*q2k1 - 8*q1t(vb)*q2U*q2tk2t*P4*m1*m2 - 4*q1t(vb)*q2U*
     & q2tk2t*P4*q1k1 - 4*q1t(vb)*q2U*q2tk2t*P4*q2k1 + 8*q1t(vb)*
     & q2U*q2tk2t*P4*k1k2 + 8*q1t(vb)*q2U*q2tk2t*P7*m1*m2 + 4*
     & q1t(vb)*q2U*q2tk2t*P7*q1k1 + 4*q1t(vb)*q2U*q2tk2t*P7*q2k1 -
     & 8*q1t(vb)*q2U*q2tk2t*P7*k1k2 - 4*q1t(vb)*q2U*t2*P3*m1*m2 + 4*
     & q1t(vb)*q2U*t2*P3*k1k2 - 4*q1t(vb)*q2U*t2*P4*m1*m2 + 4*q1t(vb)
     & *q2U*t2*P4*k1k2 + 4*q1t(vb)*q2U*t2*P7*m1*m2 - 4*q1t(vb)*q2U*
     & t2*P7*k1k2 - 4*q1t(vb)*q2tU*t1*P3*m1*m2 + 4*q1t(vb)*q2tU*t1*P3
     & *k1k2 + 4*q1t(vb)*q2tU*t2*P4*m1*m2 - 4*q1t(vb)*q2tU*t2*P4*k1k2
     &  - 4*q1t(vb)*q2tU*P3*m1*m2*q1q2 + 8*q1t(vb)*q2tU*P3*m1*m2*q1k1
     &  + 4*q1t(vb)*q2tU*P3*q1q2*k1k2
      A51(2,vb) = A51(2,vb) + 4*q1t(vb)*q2tU*P3*q1k1*q2k2 - 8*
     & q1t(vb)*q2tU*P3*q1k1*k1k2 - 4*q1t(vb)*q2tU*P3*q1k2*q2k1 + 4*
     & q1t(vb)*q2tU*P4*m1*m2*q1q2 - 8*q1t(vb)*q2tU*P4*m1*m2*q2k2 - 4*
     & q1t(vb)*q2tU*P4*q1q2*k1k2 - 4*q1t(vb)*q2tU*P4*q1k1*q2k2 + 4*
     & q1t(vb)*q2tU*P4*q1k2*q2k1 + 8*q1t(vb)*q2tU*P4*q2k2*k1k2 - 4*
     & q1t(vb)*q2tU*P7*m1*m2*q1q2 + 4*q1t(vb)*q2tU*P7*q1q2*k1k2 + 4*
     & q1t(vb)*q2tU*P7*q1k1*q2k2 - 4*q1t(vb)*q2tU*P7*q1k2*q2k1 - 4*
     & q2t(vb)*q2U*q1tq2t*P3*m1*m2 + 4*q2t(vb)*q2U*q1tq2t*P3*k1k2
     &  - 4*q2t(vb)*q2U*q1tq2t*P4*m1*m2 + 4*q2t(vb)*q2U*q1tq2t*P4*
     & k1k2 + 4*q2t(vb)*q2U*q1tq2t*P7*m1*m2 - 4*q2t(vb)*q2U*q1tq2t*
     & P7*k1k2 + 8*q2t(vb)*q2U*q1tk1t*P3*m1*m2 + 4*q2t(vb)*q2U*
     & q1tk1t*P3*q1k2 + 4*q2t(vb)*q2U*q1tk1t*P3*q2k2 - 8*q2t(vb)*
     & q2U*q1tk1t*P3*k1k2 - 4*q2t(vb)*q2U*q1tk1t*P4*q1k2 - 4*
     & q2t(vb)*q2U*q1tk1t*P4*q2k2 - 8*q2t(vb)*q2U*q1tk1t*P7*m1*m2
     &  - 4*q2t(vb)*q2U*q1tk1t*P7*q1k2 - 4*q2t(vb)*q2U*q1tk1t*P7*
     & q2k2
      A51(2,vb) = A51(2,vb) + 8*q2t(vb)*q2U*q1tk1t*P7*k1k2 - 4*
     & q2t(vb)*q2U*q1tk2t*P3*q1k1 - 4*q2t(vb)*q2U*q1tk2t*P3*q2k1 + 
     & 8*q2t(vb)*q2U*q1tk2t*P4*m1*m2 + 4*q2t(vb)*q2U*q1tk2t*P4*q1k1
     &  + 4*q2t(vb)*q2U*q1tk2t*P4*q2k1 - 8*q2t(vb)*q2U*q1tk2t*P4*
     & k1k2 + 4*q2t(vb)*q2U*q1tk2t*P7*q1k1 + 4*q2t(vb)*q2U*q1tk2t*
     & P7*q2k1 + 4*q2t(vb)*q2U*t1*P3*m1*m2 - 4*q2t(vb)*q2U*t1*P3*k1k2
     &  + 4*q2t(vb)*q2U*t1*P4*m1*m2 - 4*q2t(vb)*q2U*t1*P4*k1k2 - 4*
     & q2t(vb)*q2U*t1*P7*m1*m2 + 4*q2t(vb)*q2U*t1*P7*k1k2 + 4*q2t(vb)
     & *q1tU*t1*P3*m1*m2 - 4*q2t(vb)*q1tU*t1*P3*k1k2 - 4*q2t(vb)*
     & q1tU*t2*P4*m1*m2 + 4*q2t(vb)*q1tU*t2*P4*k1k2 + 4*q2t(vb)*q1tU
     & *P3*m1*m2*q1q2 - 8*q2t(vb)*q1tU*P3*m1*m2*q1k1 - 4*q2t(vb)*q1tU
     & *P3*q1q2*k1k2 - 4*q2t(vb)*q1tU*P3*q1k1*q2k2 + 8*q2t(vb)*q1tU*
     & P3*q1k1*k1k2 + 4*q2t(vb)*q1tU*P3*q1k2*q2k1 - 4*q2t(vb)*q1tU*P4
     & *m1*m2*q1q2 + 8*q2t(vb)*q1tU*P4*m1*m2*q2k2 + 4*q2t(vb)*q1tU*P4
     & *q1q2*k1k2 + 4*q2t(vb)*q1tU*P4*q1k1*q2k2 - 4*q2t(vb)*q1tU*P4*
     & q1k2*q2k1
      A51(2,vb) = A51(2,vb) - 8*q2t(vb)*q1tU*P4*q2k2*k1k2 - 4*
     & q2t(vb)*q1tU*P7*m1*m2*q1q2 + 4*q2t(vb)*q1tU*P7*q1q2*k1k2 + 4*
     & q2t(vb)*q1tU*P7*q1k1*q2k2 - 4*q2t(vb)*q1tU*P7*q1k2*q2k1 - 4*
     & U(vb)*t1*t2*P3*m1*m2 + 4*U(vb)*t1*t2*P3*k1k2 - 4*U(vb)*t1*t2*P4*
     & m1*m2 + 4*U(vb)*t1*t2*P4*k1k2 + 4*U(vb)*t1*t2*P7*m1*m2 - 4*U(vb)
     & *t1*t2*P7*k1k2 + 16*U(vb)*q1tq2t*q1tk1t*P3*m1*m2 + 8*U(vb)*
     & q1tq2t*q1tk1t*P3*q1k2 + 4*U(vb)*q1tq2t*q1tk1t*P3*q2k2 - 16*
     & U(vb)*q1tq2t*q1tk1t*P3*k1k2 + 4*U(vb)*q1tq2t*q1tk1t*P4*q2k2
     &  + 4*U(vb)*q1tq2t*q1tk1t*P7*q2k2 - 8*U(vb)*q1tq2t*q1tk2t*P3*
     & q1k1 - 4*U(vb)*q1tq2t*q1tk2t*P3*q2k1 - 4*U(vb)*q1tq2t*q1tk2t
     & *P4*q2k1 - 4*U(vb)*q1tq2t*q1tk2t*P7*q2k1 - 4*U(vb)*q1tq2t*
     & q2tk1t*P3*q1k2 - 4*U(vb)*q1tq2t*q2tk1t*P4*q1k2 - 8*U(vb)*
     & q1tq2t*q2tk1t*P4*q2k2 - 4*U(vb)*q1tq2t*q2tk1t*P7*q1k2 + 4*
     & U(vb)*q1tq2t*q2tk2t*P3*q1k1 + 16*U(vb)*q1tq2t*q2tk2t*P4*m1*
     & m2 + 4*U(vb)*q1tq2t*q2tk2t*P4*q1k1 + 8*U(vb)*q1tq2t*q2tk2t*
     & P4*q2k1
      A51(2,vb) = A51(2,vb) - 16*U(vb)*q1tq2t*q2tk2t*P4*k1k2 + 4*
     & U(vb)*q1tq2t*q2tk2t*P7*q1k1 + 4*U(vb)*q1tq2t*t1*P3*m1*m2 - 4*
     & U(vb)*q1tq2t*t1*P3*k1k2 + 4*U(vb)*q1tq2t*t2*P4*m1*m2 - 4*U(vb)
     & *q1tq2t*t2*P4*k1k2 - 4*U(vb)*q1tq2t*P3*m1*m2*q1q2 + 8*U(vb)*
     & q1tq2t*P3*m1*m2*q1k1 + 4*U(vb)*q1tq2t*P3*q1q2*k1k2 + 4*U(vb)*
     & q1tq2t*P3*q1k1*q2k2 - 8*U(vb)*q1tq2t*P3*q1k1*k1k2 - 4*U(vb)*
     & q1tq2t*P3*q1k2*q2k1 - 4*U(vb)*q1tq2t*P4*m1*m2*q1q2 + 8*U(vb)*
     & q1tq2t*P4*m1*m2*q2k2 + 4*U(vb)*q1tq2t*P4*q1q2*k1k2 + 4*U(vb)*
     & q1tq2t*P4*q1k1*q2k2 - 4*U(vb)*q1tq2t*P4*q1k2*q2k1 - 8*U(vb)*
     & q1tq2t*P4*q2k2*k1k2 - 4*U(vb)*q1tq2t*P7*m1*m2*q1q2 + 4*U(vb)*
     & q1tq2t*P7*q1q2*k1k2 + 4*U(vb)*q1tq2t*P7*q1k1*q2k2 - 4*U(vb)*
     & q1tq2t*P7*q1k2*q2k1 - 4*U(vb)*q1tq2t**2*P3*m1*m2 + 4*U(vb)*
     & q1tq2t**2*P3*k1k2 - 4*U(vb)*q1tq2t**2*P4*m1*m2 + 4*U(vb)*
     & q1tq2t**2*P4*k1k2 - 4*U(vb)*q1tq2t**2*P7*m1*m2 + 4*U(vb)*
     & q1tq2t**2*P7*k1k2 - 16*U(vb)*q1tk1t*q2tk1t*P3*m1*m2 - 8*U(vb)
     & *q1tk1t*q2tk1t*P3*q1k2
      A51(2,vb) = A51(2,vb) - 8*U(vb)*q1tk1t*q2tk1t*P3*q2k2 + 16*
     & U(vb)*q1tk1t*q2tk1t*P3*k1k2 - 8*U(vb)*q1tk1t*q2tk1t*P7*q2k2
     &  - 4*U(vb)*q1tk1t*q2tk2t*t1*P3 - 4*U(vb)*q1tk1t*q2tk2t*t2*P4
     &  - 4*U(vb)*q1tk1t*q2tk2t*P3*q1q2 + 8*U(vb)*q1tk1t*q2tk2t*P3*
     & q1k1 + 8*U(vb)*q1tk1t*q2tk2t*P3*q2k1 - 4*U(vb)*q1tk1t*q2tk2t
     & *P4*q1q2 + 8*U(vb)*q1tk1t*q2tk2t*P4*q1k2 + 8*U(vb)*q1tk1t*
     & q2tk2t*P4*q2k2 + 16*U(vb)*q1tk1t*q2tk2t*P7*m1*m2 - 4*U(vb)*
     & q1tk1t*q2tk2t*P7*q1q2 + 8*U(vb)*q1tk1t*q2tk2t*P7*q1k2 + 8*
     & U(vb)*q1tk1t*q2tk2t*P7*q2k1 - 16*U(vb)*q1tk1t*q2tk2t*P7*k1k2
     &  - 8*U(vb)*q1tk1t*t2*P3*m1*m2 - 4*U(vb)*q1tk1t*t2*P3*q1k2 + 8*
     & U(vb)*q1tk1t*t2*P3*k1k2 + 4*U(vb)*q1tk1t*t2*P4*q1k2 + 8*U(vb)*
     & q1tk1t*t2*P7*m1*m2 + 4*U(vb)*q1tk1t*t2*P7*q1k2 - 8*U(vb)*
     & q1tk1t*t2*P7*k1k2 + 4*U(vb)*q1tk2t*q2tk1t*t1*P3 + 4*U(vb)*
     & q1tk2t*q2tk1t*t2*P4 + 4*U(vb)*q1tk2t*q2tk1t*P3*q1q2 + 4*
     & U(vb)*q1tk2t*q2tk1t*P4*q1q2 + 4*U(vb)*q1tk2t*q2tk1t*P7*q1q2
     &  - 16*U(vb)*q1tk2t*q2tk2t*P4*m1*m2
      A51(2,vb) = A51(2,vb) - 8*U(vb)*q1tk2t*q2tk2t*P4*q1k1 - 8*
     & U(vb)*q1tk2t*q2tk2t*P4*q2k1 + 16*U(vb)*q1tk2t*q2tk2t*P4*k1k2
     &  - 8*U(vb)*q1tk2t*q2tk2t*P7*q1k1 + 4*U(vb)*q1tk2t*t2*P3*q1k1
     &  - 8*U(vb)*q1tk2t*t2*P4*m1*m2 - 4*U(vb)*q1tk2t*t2*P4*q1k1 + 8*
     & U(vb)*q1tk2t*t2*P4*k1k2 - 4*U(vb)*q1tk2t*t2*P7*q1k1 - 8*U(vb)*
     & q2tk1t*t1*P3*m1*m2 - 4*U(vb)*q2tk1t*t1*P3*q2k2 + 8*U(vb)*
     & q2tk1t*t1*P3*k1k2 + 4*U(vb)*q2tk1t*t1*P4*q2k2 - 4*U(vb)*
     & q2tk1t*t1*P7*q2k2 + 4*U(vb)*q2tk2t*t1*P3*q2k1 - 8*U(vb)*
     & q2tk2t*t1*P4*m1*m2 - 4*U(vb)*q2tk2t*t1*P4*q2k1 + 8*U(vb)*
     & q2tk2t*t1*P4*k1k2 + 8*U(vb)*q2tk2t*t1*P7*m1*m2 + 4*U(vb)*
     & q2tk2t*t1*P7*q2k1 - 8*U(vb)*q2tk2t*t1*P7*k1k2

      A51(3,vb) =
     & 8*k1(vb)*q1L*q1tq2t*q1tk2t*P3 - 4*k1(vb)*q1L*q1tq2t*q2tk2t
     & *P3 - 4*k1(vb)*q1L*q1tq2t*q2tk2t*P4 - 4*k1(vb)*q1L*q1tq2t*
     & q2tk2t*P7 - 4*k1(vb)*q1L*q1tq2t*P3*m1*m2 - 4*k1(vb)*q1L*
     & q1tq2t*P3*q2k2 + 4*k1(vb)*q1L*q1tq2t*P3*k1k2 - 4*k1(vb)*q1L*
     & q1tq2t*P4*q2k2 - 4*k1(vb)*q1L*q1tq2t*P7*q2k2 - 4*k1(vb)*q1L*
     & q1tk1t*q2tk2t*P3 - 4*k1(vb)*q1L*q1tk2t*q2tk1t*P3 + 8*k1(vb)
     & *q1L*q1tk2t*q2tk2t*P4 + 8*k1(vb)*q1L*q1tk2t*q2tk2t*P7 - 4*
     & k1(vb)*q1L*q1tk2t*t2*P3 + 4*k1(vb)*q1L*q1tk2t*t2*P4 + 4*
     & k1(vb)*q1L*q1tk2t*t2*P7 + 4*k1(vb)*q2L*q1tq2t*q1tk2t*P3 + 4
     & *k1(vb)*q2L*q1tq2t*q1tk2t*P4 + 4*k1(vb)*q2L*q1tq2t*q1tk2t*
     & P7 - 8*k1(vb)*q2L*q1tq2t*q2tk2t*P4 + 4*k1(vb)*q2L*q1tq2t*P3
     & *q1k2 - 4*k1(vb)*q2L*q1tq2t*P4*m1*m2 + 4*k1(vb)*q2L*q1tq2t*
     & P4*q1k2 + 4*k1(vb)*q2L*q1tq2t*P4*k1k2 + 4*k1(vb)*q2L*q1tq2t*
     & P7*q1k2 - 8*k1(vb)*q2L*q1tk1t*q2tk2t*P3 - 4*k1(vb)*q2L*
     & q1tk1t*q2tk2t*P4
      A51(3,vb) = A51(3,vb) - 8*k1(vb)*q2L*q1tk1t*q2tk2t*P7 + 4*
     & k1(vb)*q2L*q1tk2t*q2tk1t*P4 + 8*k1(vb)*q2L*q1tk2t*q2tk2t*
     & P4 - 4*k1(vb)*q2L*q2tk2t*t1*P3 + 4*k1(vb)*q2L*q2tk2t*t1*P4
     &  - 4*k1(vb)*q2L*q2tk2t*t1*P7 + 4*k1(vb)*q1tq2t*q1tL*P3*m1*m2
     &  + 8*k1(vb)*q1tq2t*q1tL*P3*q1k2 + 4*k1(vb)*q1tq2t*q1tL*P3*
     & q2k2 - 4*k1(vb)*q1tq2t*q1tL*P3*k1k2 + 4*k1(vb)*q1tq2t*q1tL*
     & P4*q2k2 + 4*k1(vb)*q1tq2t*q1tL*P7*q2k2 - 4*k1(vb)*q1tq2t*
     & q2tL*P3*q1k2 + 4*k1(vb)*q1tq2t*q2tL*P4*m1*m2 - 4*k1(vb)*
     & q1tq2t*q2tL*P4*q1k2 - 8*k1(vb)*q1tq2t*q2tL*P4*q2k2 - 4*
     & k1(vb)*q1tq2t*q2tL*P4*k1k2 - 4*k1(vb)*q1tq2t*q2tL*P7*q1k2 - 
     & 8*k1(vb)*q1tk1t*q2tL*P3*m1*m2 - 4*k1(vb)*q1tk1t*q2tL*P3*q1k2
     &  - 8*k1(vb)*q1tk1t*q2tL*P3*q2k2 + 8*k1(vb)*q1tk1t*q2tL*P3*
     & k1k2 - 4*k1(vb)*q1tk1t*q2tL*P4*q2k2 - 8*k1(vb)*q1tk1t*q2tL*
     & P7*q2k2 + 4*k1(vb)*q1tk2t*q2tL*t1*P3 + 4*k1(vb)*q1tk2t*q2tL*
     & t2*P4 + 4*k1(vb)*q1tk2t*q2tL*P3*q1q2 - 4*k1(vb)*q1tk2t*q2tL*
     & P3*q1k1
      A51(3,vb) = A51(3,vb) + 4*k1(vb)*q1tk2t*q2tL*P4*q1q2 + 4*
     & k1(vb)*q1tk2t*q2tL*P4*q2k1 + 4*k1(vb)*q1tk2t*q2tL*P7*q1q2 - 
     & 4*k1(vb)*q1tL*q2tk1t*P3*q1k2 + 4*k1(vb)*q1tL*q2tk1t*P4*q2k2
     &  - 4*k1(vb)*q1tL*q2tk2t*t1*P3 - 4*k1(vb)*q1tL*q2tk2t*t2*P4
     &  - 4*k1(vb)*q1tL*q2tk2t*P3*q1q2 + 4*k1(vb)*q1tL*q2tk2t*P3*
     & q1k1 - 8*k1(vb)*q1tL*q2tk2t*P4*m1*m2 - 4*k1(vb)*q1tL*q2tk2t*
     & P4*q1q2 + 8*k1(vb)*q1tL*q2tk2t*P4*q1k2 - 4*k1(vb)*q1tL*
     & q2tk2t*P4*q2k1 + 8*k1(vb)*q1tL*q2tk2t*P4*q2k2 + 8*k1(vb)*
     & q1tL*q2tk2t*P4*k1k2 - 4*k1(vb)*q1tL*q2tk2t*P7*q1q2 + 8*
     & k1(vb)*q1tL*q2tk2t*P7*q1k2 - 4*k1(vb)*q1tL*t2*P3*q1k2 - 4*
     & k1(vb)*q1tL*t2*P4*m1*m2 + 4*k1(vb)*q1tL*t2*P4*q1k2 + 4*k1(vb)*
     & q1tL*t2*P4*k1k2 + 4*k1(vb)*q1tL*t2*P7*q1k2 - 4*k1(vb)*q2tL*t1
     & *P3*m1*m2 - 4*k1(vb)*q2tL*t1*P3*q2k2 + 4*k1(vb)*q2tL*t1*P3*
     & k1k2 + 4*k1(vb)*q2tL*t1*P4*q2k2 - 4*k1(vb)*q2tL*t1*P7*q2k2 - 8
     & *k2(vb)*q1L*q1tq2t*q1tk1t*P3 + 4*k2(vb)*q1L*q1tq2t*q2tk1t*
     & P3
      A51(3,vb) = A51(3,vb) + 4*k2(vb)*q1L*q1tq2t*q2tk1t*P4 + 4*
     & k2(vb)*q1L*q1tq2t*q2tk1t*P7 - 4*k2(vb)*q1L*q1tq2t*P3*m1*m2
     &  + 4*k2(vb)*q1L*q1tq2t*P3*q2k1 + 4*k2(vb)*q1L*q1tq2t*P3*k1k2
     &  + 4*k2(vb)*q1L*q1tq2t*P4*q2k1 + 4*k2(vb)*q1L*q1tq2t*P7*q2k1
     &  + 8*k2(vb)*q1L*q1tk1t*q2tk1t*P3 - 4*k2(vb)*q1L*q1tk1t*
     & q2tk2t*P3 - 8*k2(vb)*q1L*q1tk1t*q2tk2t*P4 - 8*k2(vb)*q1L*
     & q1tk1t*q2tk2t*P7 + 4*k2(vb)*q1L*q1tk1t*t2*P3 - 4*k2(vb)*q1L
     & *q1tk1t*t2*P4 - 4*k2(vb)*q1L*q1tk1t*t2*P7 + 4*k2(vb)*q1L*
     & q1tk2t*q2tk1t*P3 - 4*k2(vb)*q2L*q1tq2t*q1tk1t*P3 - 4*k2(vb)
     & *q2L*q1tq2t*q1tk1t*P4 - 4*k2(vb)*q2L*q1tq2t*q1tk1t*P7 + 8*
     & k2(vb)*q2L*q1tq2t*q2tk1t*P4 - 4*k2(vb)*q2L*q1tq2t*P3*q1k1
     &  - 4*k2(vb)*q2L*q1tq2t*P4*m1*m2 - 4*k2(vb)*q2L*q1tq2t*P4*
     & q1k1 + 4*k2(vb)*q2L*q1tq2t*P4*k1k2 - 4*k2(vb)*q2L*q1tq2t*P7*
     & q1k1 + 8*k2(vb)*q2L*q1tk1t*q2tk1t*P3 + 8*k2(vb)*q2L*q1tk1t*
     & q2tk1t*P7 - 4*k2(vb)*q2L*q1tk1t*q2tk2t*P4 - 4*k2(vb)*q2L*
     & q1tk2t*q2tk1t*P4
      A51(3,vb) = A51(3,vb) + 4*k2(vb)*q2L*q2tk1t*t1*P3 - 4*k2(vb)*
     & q2L*q2tk1t*t1*P4 + 4*k2(vb)*q2L*q2tk1t*t1*P7 + 4*k2(vb)*
     & q1tq2t*q1tL*P3*m1*m2 - 8*k2(vb)*q1tq2t*q1tL*P3*q1k1 - 4*
     & k2(vb)*q1tq2t*q1tL*P3*q2k1 - 4*k2(vb)*q1tq2t*q1tL*P3*k1k2 - 
     & 4*k2(vb)*q1tq2t*q1tL*P4*q2k1 - 4*k2(vb)*q1tq2t*q1tL*P7*q2k1
     &  + 4*k2(vb)*q1tq2t*q2tL*P3*q1k1 + 4*k2(vb)*q1tq2t*q2tL*P4*m1
     & *m2 + 4*k2(vb)*q1tq2t*q2tL*P4*q1k1 + 8*k2(vb)*q1tq2t*q2tL*P4
     & *q2k1 - 4*k2(vb)*q1tq2t*q2tL*P4*k1k2 + 4*k2(vb)*q1tq2t*q2tL*
     & P7*q1k1 - 4*k2(vb)*q1tk1t*q2tL*t1*P3 - 4*k2(vb)*q1tk1t*q2tL*
     & t2*P4 - 8*k2(vb)*q1tk1t*q2tL*P3*m1*m2 - 4*k2(vb)*q1tk1t*q2tL
     & *P3*q1q2 + 8*k2(vb)*q1tk1t*q2tL*P3*q1k1 - 4*k2(vb)*q1tk1t*
     & q2tL*P3*q1k2 + 8*k2(vb)*q1tk1t*q2tL*P3*q2k1 + 8*k2(vb)*
     & q1tk1t*q2tL*P3*k1k2 - 4*k2(vb)*q1tk1t*q2tL*P4*q1q2 + 4*
     & k2(vb)*q1tk1t*q2tL*P4*q2k2 - 4*k2(vb)*q1tk1t*q2tL*P7*q1q2 + 
     & 8*k2(vb)*q1tk1t*q2tL*P7*q2k1 + 4*k2(vb)*q1tk2t*q2tL*P3*q1k1
     &  - 4*k2(vb)*q1tk2t*q2tL*P4*q2k1
      A51(3,vb) = A51(3,vb) + 4*k2(vb)*q1tL*q2tk1t*t1*P3 + 4*k2(vb)
     & *q1tL*q2tk1t*t2*P4 + 4*k2(vb)*q1tL*q2tk1t*P3*q1q2 + 4*k2(vb)
     & *q1tL*q2tk1t*P3*q1k2 + 4*k2(vb)*q1tL*q2tk1t*P4*q1q2 - 4*
     & k2(vb)*q1tL*q2tk1t*P4*q2k2 + 4*k2(vb)*q1tL*q2tk1t*P7*q1q2 - 
     & 4*k2(vb)*q1tL*q2tk2t*P3*q1k1 - 8*k2(vb)*q1tL*q2tk2t*P4*m1*m2
     &  - 8*k2(vb)*q1tL*q2tk2t*P4*q1k1 - 4*k2(vb)*q1tL*q2tk2t*P4*
     & q2k1 + 8*k2(vb)*q1tL*q2tk2t*P4*k1k2 - 8*k2(vb)*q1tL*q2tk2t*
     & P7*q1k1 + 4*k2(vb)*q1tL*t2*P3*q1k1 - 4*k2(vb)*q1tL*t2*P4*m1*m2
     &  - 4*k2(vb)*q1tL*t2*P4*q1k1 + 4*k2(vb)*q1tL*t2*P4*k1k2 - 4*
     & k2(vb)*q1tL*t2*P7*q1k1 - 4*k2(vb)*q2tL*t1*P3*m1*m2 + 4*k2(vb)*
     & q2tL*t1*P3*q2k1 + 4*k2(vb)*q2tL*t1*P3*k1k2 - 4*k2(vb)*q2tL*t1
     & *P4*q2k1 + 4*k2(vb)*q2tL*t1*P7*q2k1 - 4*kw(vb)*q1L*q1tq2t*P3*
     & m1*m2 + 4*kw(vb)*q1L*q1tq2t*P3*k1k2 - 4*kw(vb)*q1L*q1tk1t*
     & q2tk2t*P3 + 4*kw(vb)*q1L*q1tk2t*q2tk1t*P3 - 4*kw(vb)*q2L*
     & q1tq2t*P4*m1*m2 + 4*kw(vb)*q2L*q1tq2t*P4*k1k2 - 4*kw(vb)*q2L
     & *q1tk1t*q2tk2t*P4
      A51(3,vb) = A51(3,vb) + 4*kw(vb)*q2L*q1tk2t*q2tk1t*P4 + 4*
     & kw(vb)*q1tq2t*q1tL*P3*m1*m2 - 4*kw(vb)*q1tq2t*q1tL*P3*k1k2
     &  + 4*kw(vb)*q1tq2t*q2tL*P4*m1*m2 - 4*kw(vb)*q1tq2t*q2tL*P4*
     & k1k2 - 8*kw(vb)*q1tk1t*q2tL*P3*m1*m2 - 4*kw(vb)*q1tk1t*q2tL*
     & P3*q1k2 + 8*kw(vb)*q1tk1t*q2tL*P3*k1k2 - 4*kw(vb)*q1tk1t*
     & q2tL*P4*q2k2 + 4*kw(vb)*q1tk2t*q2tL*P3*q1k1 + 4*kw(vb)*
     & q1tk2t*q2tL*P4*q2k1 + 4*kw(vb)*q1tL*q2tk1t*P3*q1k2 + 4*
     & kw(vb)*q1tL*q2tk1t*P4*q2k2 - 4*kw(vb)*q1tL*q2tk2t*P3*q1k1 - 
     & 8*kw(vb)*q1tL*q2tk2t*P4*m1*m2 - 4*kw(vb)*q1tL*q2tk2t*P4*q2k1
     &  + 8*kw(vb)*q1tL*q2tk2t*P4*k1k2 - 4*kw(vb)*q1tL*t2*P4*m1*m2
     &  + 4*kw(vb)*q1tL*t2*P4*k1k2 - 4*kw(vb)*q2tL*t1*P3*m1*m2 + 4*
     & kw(vb)*q2tL*t1*P3*k1k2 + 4*q1(vb)*q1L*q1tq2t*P3*m1*m2 - 4*
     & q1(vb)*q1L*q1tq2t*P3*k1k2 + 4*q1(vb)*q1L*q1tk1t*q2tk2t*P3
     &  - 4*q1(vb)*q1L*q1tk2t*q2tk1t*P3 + 4*q1(vb)*q2L*q1tq2t*P3*
     & m1*m2 - 4*q1(vb)*q2L*q1tq2t*P3*k1k2 + 4*q1(vb)*q2L*q1tq2t*P7
     & *m1*m2
      A51(3,vb) = A51(3,vb) - 4*q1(vb)*q2L*q1tq2t*P7*k1k2 + 4*
     & q1(vb)*q2L*q1tk1t*q2tk2t*P3 + 4*q1(vb)*q2L*q1tk1t*q2tk2t*
     & P7 - 4*q1(vb)*q2L*q1tk2t*q2tk1t*P3 - 4*q1(vb)*q2L*q1tk2t*
     & q2tk1t*P7 + 4*q1(vb)*q1tq2t*q1tL*P3*m1*m2 - 4*q1(vb)*q1tq2t*
     & q1tL*P3*k1k2 - 4*q1(vb)*q1tq2t*q2tL*P3*m1*m2 + 4*q1(vb)*
     & q1tq2t*q2tL*P3*k1k2 - 4*q1(vb)*q1tq2t*q2tL*P7*m1*m2 + 4*
     & q1(vb)*q1tq2t*q2tL*P7*k1k2 + 8*q1(vb)*q1tk1t*q2tL*P3*m1*m2
     &  + 4*q1(vb)*q1tk1t*q2tL*P3*q1k2 + 4*q1(vb)*q1tk1t*q2tL*P3*
     & q2k2 - 8*q1(vb)*q1tk1t*q2tL*P3*k1k2 + 4*q1(vb)*q1tk1t*q2tL*
     & P7*q2k2 - 4*q1(vb)*q1tk2t*q2tL*P3*q1k1 - 4*q1(vb)*q1tk2t*
     & q2tL*P3*q2k1 - 4*q1(vb)*q1tk2t*q2tL*P7*q2k1 - 8*q1(vb)*q1tL*
     & q2tk1t*P3*m1*m2 - 4*q1(vb)*q1tL*q2tk1t*P3*q1k2 - 4*q1(vb)*
     & q1tL*q2tk1t*P3*q2k2 + 8*q1(vb)*q1tL*q2tk1t*P3*k1k2 - 4*
     & q1(vb)*q1tL*q2tk1t*P7*q2k2 + 4*q1(vb)*q1tL*q2tk2t*P3*q1k1 + 
     & 4*q1(vb)*q1tL*q2tk2t*P3*q2k1 + 8*q1(vb)*q1tL*q2tk2t*P7*m1*m2
     &  + 4*q1(vb)*q1tL*q2tk2t*P7*q2k1
      A51(3,vb) = A51(3,vb) - 8*q1(vb)*q1tL*q2tk2t*P7*k1k2 - 4*
     & q1(vb)*q1tL*t2*P3*m1*m2 + 4*q1(vb)*q1tL*t2*P3*k1k2 + 4*q1(vb)*
     & q1tL*t2*P7*m1*m2 - 4*q1(vb)*q1tL*t2*P7*k1k2 + 4*q1(vb)*q2tL*
     & t1*P3*m1*m2 - 4*q1(vb)*q2tL*t1*P3*k1k2 + 4*q2(vb)*q1L*q1tq2t*
     & P4*m1*m2 - 4*q2(vb)*q1L*q1tq2t*P4*k1k2 + 4*q2(vb)*q1L*q1tq2t
     & *P7*m1*m2 - 4*q2(vb)*q1L*q1tq2t*P7*k1k2 + 4*q2(vb)*q1L*
     & q1tk1t*q2tk2t*P4 + 4*q2(vb)*q1L*q1tk1t*q2tk2t*P7 - 4*q2(vb)
     & *q1L*q1tk2t*q2tk1t*P4 - 4*q2(vb)*q1L*q1tk2t*q2tk1t*P7 + 4*
     & q2(vb)*q2L*q1tq2t*P4*m1*m2 - 4*q2(vb)*q2L*q1tq2t*P4*k1k2 + 4
     & *q2(vb)*q2L*q1tk1t*q2tk2t*P4 - 4*q2(vb)*q2L*q1tk2t*q2tk1t*
     & P4 - 4*q2(vb)*q1tq2t*q1tL*P4*m1*m2 + 4*q2(vb)*q1tq2t*q1tL*P4
     & *k1k2 - 4*q2(vb)*q1tq2t*q1tL*P7*m1*m2 + 4*q2(vb)*q1tq2t*q1tL
     & *P7*k1k2 + 4*q2(vb)*q1tq2t*q2tL*P4*m1*m2 - 4*q2(vb)*q1tq2t*
     & q2tL*P4*k1k2 + 4*q2(vb)*q1tk1t*q2tL*P4*q1k2 + 4*q2(vb)*
     & q1tk1t*q2tL*P4*q2k2 + 8*q2(vb)*q1tk1t*q2tL*P7*m1*m2 + 4*
     & q2(vb)*q1tk1t*q2tL*P7*q1k2
      A51(3,vb) = A51(3,vb) - 8*q2(vb)*q1tk1t*q2tL*P7*k1k2 - 8*
     & q2(vb)*q1tk2t*q2tL*P4*m1*m2 - 4*q2(vb)*q1tk2t*q2tL*P4*q1k1
     &  - 4*q2(vb)*q1tk2t*q2tL*P4*q2k1 + 8*q2(vb)*q1tk2t*q2tL*P4*
     & k1k2 - 4*q2(vb)*q1tk2t*q2tL*P7*q1k1 - 4*q2(vb)*q1tL*q2tk1t*
     & P4*q1k2 - 4*q2(vb)*q1tL*q2tk1t*P4*q2k2 - 4*q2(vb)*q1tL*
     & q2tk1t*P7*q1k2 + 8*q2(vb)*q1tL*q2tk2t*P4*m1*m2 + 4*q2(vb)*
     & q1tL*q2tk2t*P4*q1k1 + 4*q2(vb)*q1tL*q2tk2t*P4*q2k1 - 8*
     & q2(vb)*q1tL*q2tk2t*P4*k1k2 + 4*q2(vb)*q1tL*q2tk2t*P7*q1k1 + 
     & 4*q2(vb)*q1tL*t2*P4*m1*m2 - 4*q2(vb)*q1tL*t2*P4*k1k2 - 4*
     & q2(vb)*q2tL*t1*P4*m1*m2 + 4*q2(vb)*q2tL*t1*P4*k1k2 + 4*q2(vb)*
     & q2tL*t1*P7*m1*m2 - 4*q2(vb)*q2tL*t1*P7*k1k2 - 8*q1t(vb)*q1L*
     & q1tq2t*P3*m1*m2 + 8*q1t(vb)*q1L*q1tq2t*P3*k1k2 + 8*q1t(vb)*
     & q1L*q2tk1t*P3*m1*m2 + 4*q1t(vb)*q1L*q2tk1t*P3*q2k2 - 8*
     & q1t(vb)*q1L*q2tk1t*P3*k1k2 - 4*q1t(vb)*q1L*q2tk1t*P4*q2k2 + 
     & 4*q1t(vb)*q1L*q2tk1t*P7*q2k2 - 4*q1t(vb)*q1L*q2tk2t*P3*q2k1
     &  + 8*q1t(vb)*q1L*q2tk2t*P4*m1*m2
      A51(3,vb) = A51(3,vb) + 4*q1t(vb)*q1L*q2tk2t*P4*q2k1 - 8*
     & q1t(vb)*q1L*q2tk2t*P4*k1k2 - 8*q1t(vb)*q1L*q2tk2t*P7*m1*m2
     &  - 4*q1t(vb)*q1L*q2tk2t*P7*q2k1 + 8*q1t(vb)*q1L*q2tk2t*P7*
     & k1k2 + 4*q1t(vb)*q1L*t2*P3*m1*m2 - 4*q1t(vb)*q1L*t2*P3*k1k2 + 
     & 4*q1t(vb)*q1L*t2*P4*m1*m2 - 4*q1t(vb)*q1L*t2*P4*k1k2 - 4*
     & q1t(vb)*q1L*t2*P7*m1*m2 + 4*q1t(vb)*q1L*t2*P7*k1k2 - 4*q1t(vb)
     & *q2L*q1tq2t*P3*m1*m2 + 4*q1t(vb)*q2L*q1tq2t*P3*k1k2 + 4*
     & q1t(vb)*q2L*q1tq2t*P4*m1*m2 - 4*q1t(vb)*q2L*q1tq2t*P4*k1k2
     &  - 4*q1t(vb)*q2L*q1tq2t*P7*m1*m2 + 4*q1t(vb)*q2L*q1tq2t*P7*
     & k1k2 - 4*q1t(vb)*q2L*q2tk1t*P3*q1k2 + 4*q1t(vb)*q2L*q2tk1t*
     & P4*q1k2 - 4*q1t(vb)*q2L*q2tk1t*P7*q1k2 + 4*q1t(vb)*q2L*
     & q2tk2t*P3*q1k1 - 4*q1t(vb)*q2L*q2tk2t*P4*q1k1 + 4*q1t(vb)*
     & q2L*q2tk2t*P7*q1k1 - 4*q1t(vb)*q2tL*t1*P3*m1*m2 + 4*q1t(vb)*
     & q2tL*t1*P3*k1k2 + 4*q1t(vb)*q2tL*t2*P4*m1*m2 - 4*q1t(vb)*q2tL
     & *t2*P4*k1k2 - 4*q1t(vb)*q2tL*P3*m1*m2*q1q2 + 8*q1t(vb)*q2tL*P3
     & *m1*m2*q1k1
      A51(3,vb) = A51(3,vb) + 4*q1t(vb)*q2tL*P3*q1q2*k1k2 + 4*
     & q1t(vb)*q2tL*P3*q1k1*q2k2 - 8*q1t(vb)*q2tL*P3*q1k1*k1k2 - 4*
     & q1t(vb)*q2tL*P3*q1k2*q2k1 + 4*q1t(vb)*q2tL*P4*m1*m2*q1q2 - 8*
     & q1t(vb)*q2tL*P4*m1*m2*q2k2 - 4*q1t(vb)*q2tL*P4*q1q2*k1k2 - 4*
     & q1t(vb)*q2tL*P4*q1k1*q2k2 + 4*q1t(vb)*q2tL*P4*q1k2*q2k1 + 8*
     & q1t(vb)*q2tL*P4*q2k2*k1k2 - 4*q1t(vb)*q2tL*P7*m1*m2*q1q2 + 4*
     & q1t(vb)*q2tL*P7*q1q2*k1k2 + 4*q1t(vb)*q2tL*P7*q1k1*q2k2 - 4*
     & q1t(vb)*q2tL*P7*q1k2*q2k1 + 4*q2t(vb)*q1L*q1tq2t*P3*m1*m2 - 4
     & *q2t(vb)*q1L*q1tq2t*P3*k1k2 - 4*q2t(vb)*q1L*q1tq2t*P4*m1*m2
     &  + 4*q2t(vb)*q1L*q1tq2t*P4*k1k2 - 4*q2t(vb)*q1L*q1tq2t*P7*m1
     & *m2 + 4*q2t(vb)*q1L*q1tq2t*P7*k1k2 - 4*q2t(vb)*q1L*q1tk1t*P3
     & *q2k2 + 4*q2t(vb)*q1L*q1tk1t*P4*q2k2 + 4*q2t(vb)*q1L*q1tk1t*
     & P7*q2k2 + 4*q2t(vb)*q1L*q1tk2t*P3*q2k1 - 4*q2t(vb)*q1L*
     & q1tk2t*P4*q2k1 - 4*q2t(vb)*q1L*q1tk2t*P7*q2k1 - 8*q2t(vb)*
     & q2L*q1tq2t*P4*m1*m2 + 8*q2t(vb)*q2L*q1tq2t*P4*k1k2 + 8*
     & q2t(vb)*q2L*q1tk1t*P3*m1*m2
      A51(3,vb) = A51(3,vb) + 4*q2t(vb)*q2L*q1tk1t*P3*q1k2 - 8*
     & q2t(vb)*q2L*q1tk1t*P3*k1k2 - 4*q2t(vb)*q2L*q1tk1t*P4*q1k2 - 
     & 8*q2t(vb)*q2L*q1tk1t*P7*m1*m2 - 4*q2t(vb)*q2L*q1tk1t*P7*q1k2
     &  + 8*q2t(vb)*q2L*q1tk1t*P7*k1k2 - 4*q2t(vb)*q2L*q1tk2t*P3*
     & q1k1 + 8*q2t(vb)*q2L*q1tk2t*P4*m1*m2 + 4*q2t(vb)*q2L*q1tk2t*
     & P4*q1k1 - 8*q2t(vb)*q2L*q1tk2t*P4*k1k2 + 4*q2t(vb)*q2L*
     & q1tk2t*P7*q1k1 + 4*q2t(vb)*q2L*t1*P3*m1*m2 - 4*q2t(vb)*q2L*t1
     & *P3*k1k2 + 4*q2t(vb)*q2L*t1*P4*m1*m2 - 4*q2t(vb)*q2L*t1*P4*
     & k1k2 - 4*q2t(vb)*q2L*t1*P7*m1*m2 + 4*q2t(vb)*q2L*t1*P7*k1k2 + 
     & 4*q2t(vb)*q1tL*t1*P3*m1*m2 - 4*q2t(vb)*q1tL*t1*P3*k1k2 - 4*
     & q2t(vb)*q1tL*t2*P4*m1*m2 + 4*q2t(vb)*q1tL*t2*P4*k1k2 + 4*
     & q2t(vb)*q1tL*P3*m1*m2*q1q2 - 8*q2t(vb)*q1tL*P3*m1*m2*q1k1 - 4*
     & q2t(vb)*q1tL*P3*q1q2*k1k2 - 4*q2t(vb)*q1tL*P3*q1k1*q2k2 + 8*
     & q2t(vb)*q1tL*P3*q1k1*k1k2 + 4*q2t(vb)*q1tL*P3*q1k2*q2k1 - 4*
     & q2t(vb)*q1tL*P4*m1*m2*q1q2 + 8*q2t(vb)*q1tL*P4*m1*m2*q2k2 + 4*
     & q2t(vb)*q1tL*P4*q1q2*k1k2
      A51(3,vb) = A51(3,vb) + 4*q2t(vb)*q1tL*P4*q1k1*q2k2 - 4*
     & q2t(vb)*q1tL*P4*q1k2*q2k1 - 8*q2t(vb)*q1tL*P4*q2k2*k1k2 - 4*
     & q2t(vb)*q1tL*P7*m1*m2*q1q2 + 4*q2t(vb)*q1tL*P7*q1q2*k1k2 + 4*
     & q2t(vb)*q1tL*P7*q1k1*q2k2 - 4*q2t(vb)*q1tL*P7*q1k2*q2k1 - 4*
     & L(vb)*t1*t2*P3*m1*m2 + 4*L(vb)*t1*t2*P3*k1k2 - 4*L(vb)*t1*t2*P4*
     & m1*m2 + 4*L(vb)*t1*t2*P4*k1k2 + 4*L(vb)*t1*t2*P7*m1*m2 - 4*L(vb)
     & *t1*t2*P7*k1k2 + 16*L(vb)*q1tq2t*q1tk1t*P3*m1*m2 + 8*L(vb)*
     & q1tq2t*q1tk1t*P3*q1k2 + 4*L(vb)*q1tq2t*q1tk1t*P3*q2k2 - 16*
     & L(vb)*q1tq2t*q1tk1t*P3*k1k2 + 4*L(vb)*q1tq2t*q1tk1t*P4*q2k2
     &  + 4*L(vb)*q1tq2t*q1tk1t*P7*q2k2 - 8*L(vb)*q1tq2t*q1tk2t*P3*
     & q1k1 - 4*L(vb)*q1tq2t*q1tk2t*P3*q2k1 - 4*L(vb)*q1tq2t*q1tk2t
     & *P4*q2k1 - 4*L(vb)*q1tq2t*q1tk2t*P7*q2k1 - 4*L(vb)*q1tq2t*
     & q2tk1t*P3*q1k2 - 4*L(vb)*q1tq2t*q2tk1t*P4*q1k2 - 8*L(vb)*
     & q1tq2t*q2tk1t*P4*q2k2 - 4*L(vb)*q1tq2t*q2tk1t*P7*q1k2 + 4*
     & L(vb)*q1tq2t*q2tk2t*P3*q1k1 + 16*L(vb)*q1tq2t*q2tk2t*P4*m1*
     & m2
      A51(3,vb) = A51(3,vb) + 4*L(vb)*q1tq2t*q2tk2t*P4*q1k1 + 8*
     & L(vb)*q1tq2t*q2tk2t*P4*q2k1 - 16*L(vb)*q1tq2t*q2tk2t*P4*k1k2
     &  + 4*L(vb)*q1tq2t*q2tk2t*P7*q1k1 + 4*L(vb)*q1tq2t*t1*P3*m1*m2
     &  - 4*L(vb)*q1tq2t*t1*P3*k1k2 + 4*L(vb)*q1tq2t*t2*P4*m1*m2 - 4*
     & L(vb)*q1tq2t*t2*P4*k1k2 - 4*L(vb)*q1tq2t*P3*m1*m2*q1q2 + 8*
     & L(vb)*q1tq2t*P3*m1*m2*q1k1 + 4*L(vb)*q1tq2t*P3*q1q2*k1k2 + 4*
     & L(vb)*q1tq2t*P3*q1k1*q2k2 - 8*L(vb)*q1tq2t*P3*q1k1*k1k2 - 4*
     & L(vb)*q1tq2t*P3*q1k2*q2k1 - 4*L(vb)*q1tq2t*P4*m1*m2*q1q2 + 8*
     & L(vb)*q1tq2t*P4*m1*m2*q2k2 + 4*L(vb)*q1tq2t*P4*q1q2*k1k2 + 4*
     & L(vb)*q1tq2t*P4*q1k1*q2k2 - 4*L(vb)*q1tq2t*P4*q1k2*q2k1 - 8*
     & L(vb)*q1tq2t*P4*q2k2*k1k2 - 4*L(vb)*q1tq2t*P7*m1*m2*q1q2 + 4*
     & L(vb)*q1tq2t*P7*q1q2*k1k2 + 4*L(vb)*q1tq2t*P7*q1k1*q2k2 - 4*
     & L(vb)*q1tq2t*P7*q1k2*q2k1 - 4*L(vb)*q1tq2t**2*P3*m1*m2 + 4*
     & L(vb)*q1tq2t**2*P3*k1k2 - 4*L(vb)*q1tq2t**2*P4*m1*m2 + 4*L(vb)
     & *q1tq2t**2*P4*k1k2 - 4*L(vb)*q1tq2t**2*P7*m1*m2 + 4*L(vb)*
     & q1tq2t**2*P7*k1k2
      A51(3,vb) = A51(3,vb) - 16*L(vb)*q1tk1t*q2tk1t*P3*m1*m2 - 8*
     & L(vb)*q1tk1t*q2tk1t*P3*q1k2 - 8*L(vb)*q1tk1t*q2tk1t*P3*q2k2
     &  + 16*L(vb)*q1tk1t*q2tk1t*P3*k1k2 - 8*L(vb)*q1tk1t*q2tk1t*P7
     & *q2k2 - 4*L(vb)*q1tk1t*q2tk2t*t1*P3 - 4*L(vb)*q1tk1t*q2tk2t*
     & t2*P4 - 4*L(vb)*q1tk1t*q2tk2t*P3*q1q2 + 8*L(vb)*q1tk1t*
     & q2tk2t*P3*q1k1 + 8*L(vb)*q1tk1t*q2tk2t*P3*q2k1 - 4*L(vb)*
     & q1tk1t*q2tk2t*P4*q1q2 + 8*L(vb)*q1tk1t*q2tk2t*P4*q1k2 + 8*
     & L(vb)*q1tk1t*q2tk2t*P4*q2k2 + 16*L(vb)*q1tk1t*q2tk2t*P7*m1*
     & m2 - 4*L(vb)*q1tk1t*q2tk2t*P7*q1q2 + 8*L(vb)*q1tk1t*q2tk2t*
     & P7*q1k2 + 8*L(vb)*q1tk1t*q2tk2t*P7*q2k1 - 16*L(vb)*q1tk1t*
     & q2tk2t*P7*k1k2 - 8*L(vb)*q1tk1t*t2*P3*m1*m2 - 4*L(vb)*q1tk1t*
     & t2*P3*q1k2 + 8*L(vb)*q1tk1t*t2*P3*k1k2 + 4*L(vb)*q1tk1t*t2*P4*
     & q1k2 + 8*L(vb)*q1tk1t*t2*P7*m1*m2 + 4*L(vb)*q1tk1t*t2*P7*q1k2
     &  - 8*L(vb)*q1tk1t*t2*P7*k1k2 + 4*L(vb)*q1tk2t*q2tk1t*t1*P3 + 
     & 4*L(vb)*q1tk2t*q2tk1t*t2*P4 + 4*L(vb)*q1tk2t*q2tk1t*P3*q1q2
     &  + 4*L(vb)*q1tk2t*q2tk1t*P4*q1q2
      A51(3,vb) = A51(3,vb) + 4*L(vb)*q1tk2t*q2tk1t*P7*q1q2 - 16*
     & L(vb)*q1tk2t*q2tk2t*P4*m1*m2 - 8*L(vb)*q1tk2t*q2tk2t*P4*q1k1
     &  - 8*L(vb)*q1tk2t*q2tk2t*P4*q2k1 + 16*L(vb)*q1tk2t*q2tk2t*P4
     & *k1k2 - 8*L(vb)*q1tk2t*q2tk2t*P7*q1k1 + 4*L(vb)*q1tk2t*t2*P3*
     & q1k1 - 8*L(vb)*q1tk2t*t2*P4*m1*m2 - 4*L(vb)*q1tk2t*t2*P4*q1k1
     &  + 8*L(vb)*q1tk2t*t2*P4*k1k2 - 4*L(vb)*q1tk2t*t2*P7*q1k1 - 8*
     & L(vb)*q2tk1t*t1*P3*m1*m2 - 4*L(vb)*q2tk1t*t1*P3*q2k2 + 8*L(vb)
     & *q2tk1t*t1*P3*k1k2 + 4*L(vb)*q2tk1t*t1*P4*q2k2 - 4*L(vb)*
     & q2tk1t*t1*P7*q2k2 + 4*L(vb)*q2tk2t*t1*P3*q2k1 - 8*L(vb)*
     & q2tk2t*t1*P4*m1*m2 - 4*L(vb)*q2tk2t*t1*P4*q2k1 + 8*L(vb)*
     & q2tk2t*t1*P4*k1k2 + 8*L(vb)*q2tk2t*t1*P7*m1*m2 + 4*L(vb)*
     & q2tk2t*t1*P7*q2k1 - 8*L(vb)*q2tk2t*t1*P7*k1k2

      A51(4,vb) =
     & 4*k1(vb)*q1L*q2U*q1tq2t*P3*m2 + 4*k1(vb)*q1L*q2U*q1tq2t*P4
     & *m2 + 4*k1(vb)*q1L*q2U*q1tq2t*P7*m2 + 8*k1(vb)*q1L*q1tq2t*
     & q1tU*P3*m2 - 4*k1(vb)*q1L*q1tq2t*q2tU*P3*m2 - 4*k1(vb)*q1L*
     & q1tq2t*q2tU*P4*m2 - 4*k1(vb)*q1L*q1tq2t*q2tU*P7*m2 - 4*
     & k1(vb)*q1L*q1tk1t*q2tU*P3*m2 - 4*k1(vb)*q1L*q1tk2t*q2tU*P3
     & *m1 - 4*k1(vb)*q1L*q1tU*q2tk1t*P3*m2 + 4*k1(vb)*q1L*q1tU*
     & q2tk2t*P3*m1 + 8*k1(vb)*q1L*q1tU*q2tk2t*P4*m2 + 8*k1(vb)*
     & q1L*q1tU*q2tk2t*P7*m2 - 4*k1(vb)*q1L*q1tU*t2*P3*m2 + 4*
     & k1(vb)*q1L*q1tU*t2*P4*m2 + 4*k1(vb)*q1L*q1tU*t2*P7*m2 + 4*
     & k1(vb)*q2U*q2L*q1tq2t*P3*m2 + 4*k1(vb)*q2U*q2L*q1tq2t*P4*
     & m2 + 4*k1(vb)*q2U*q2L*q1tq2t*P7*m2 + 4*k1(vb)*q2U*q1tq2t*
     & q1tL*P3*m2 - 4*k1(vb)*q2U*q1tq2t*q1tL*P4*m2 - 4*k1(vb)*q2U*
     & q1tq2t*q1tL*P7*m2 - 4*k1(vb)*q2U*q1tq2t*q2tL*P3*m2 + 4*
     & k1(vb)*q2U*q1tq2t*q2tL*P4*m2 - 4*k1(vb)*q2U*q1tq2t*q2tL*P7
     & *m2
      A51(4,vb) = A51(4,vb) + 4*k1(vb)*q2U*q1tk1t*q2tL*P3*m2 + 4*
     & k1(vb)*q2U*q1tk1t*q2tL*P4*m2 + 8*k1(vb)*q2U*q1tk1t*q2tL*P7
     & *m2 - 4*k1(vb)*q2U*q1tk2t*q2tL*P3*m1 + 4*k1(vb)*q2U*q1tk2t*
     & q2tL*P4*m1 - 8*k1(vb)*q2U*q1tk2t*q2tL*P4*m2 - 4*k1(vb)*q2U*
     & q1tL*q2tk1t*P3*m2 - 4*k1(vb)*q2U*q1tL*q2tk1t*P4*m2 + 4*
     & k1(vb)*q2U*q1tL*q2tk2t*P3*m1 - 4*k1(vb)*q2U*q1tL*q2tk2t*P4
     & *m1 + 8*k1(vb)*q2U*q1tL*q2tk2t*P4*m2 + 8*k1(vb)*q2U*q1tL*
     & q2tk2t*P7*m2 - 4*k1(vb)*q2U*q1tL*t2*P3*m2 + 4*k1(vb)*q2U*
     & q1tL*t2*P4*m2 + 4*k1(vb)*q2U*q1tL*t2*P7*m2 + 4*k1(vb)*q2U*
     & q2tL*t1*P3*m2 - 4*k1(vb)*q2U*q2tL*t1*P4*m2 + 4*k1(vb)*q2U*
     & q2tL*t1*P7*m2 + 4*k1(vb)*q2L*q1tq2t*q1tU*P3*m2 + 4*k1(vb)*
     & q2L*q1tq2t*q1tU*P4*m2 + 4*k1(vb)*q2L*q1tq2t*q1tU*P7*m2 - 8
     & *k1(vb)*q2L*q1tq2t*q2tU*P4*m2 - 8*k1(vb)*q2L*q1tk1t*q2tU*
     & P3*m2 - 4*k1(vb)*q2L*q1tk1t*q2tU*P4*m2 - 8*k1(vb)*q2L*
     & q1tk1t*q2tU*P7*m2 - 4*k1(vb)*q2L*q1tk2t*q2tU*P4*m1 + 8*
     & k1(vb)*q2L*q1tk2t*q2tU*P4*m2
      A51(4,vb) = A51(4,vb) + 4*k1(vb)*q2L*q1tU*q2tk1t*P4*m2 + 4*
     & k1(vb)*q2L*q1tU*q2tk2t*P4*m1 - 4*k1(vb)*q2L*q2tU*t1*P3*m2
     &  + 4*k1(vb)*q2L*q2tU*t1*P4*m2 - 4*k1(vb)*q2L*q2tU*t1*P7*m2
     &  + 4*k1(vb)*q1tU*q2tL*t1*P3*m2 + 4*k1(vb)*q1tU*q2tL*t2*P4*m2
     &  + 4*k1(vb)*q1tU*q2tL*P3*m1*q1k2 + 4*k1(vb)*q1tU*q2tL*P3*m2*
     & q1q2 - 4*k1(vb)*q1tU*q2tL*P3*m2*q1k1 + 4*k1(vb)*q1tU*q2tL*P4
     & *m1*q2k2 + 4*k1(vb)*q1tU*q2tL*P4*m2*q1q2 + 4*k1(vb)*q1tU*
     & q2tL*P4*m2*q2k1 - 8*k1(vb)*q1tU*q2tL*P4*m2*q2k2 + 4*k1(vb)*
     & q1tU*q2tL*P7*m2*q1q2 - 4*k1(vb)*q1tL*q2tU*t1*P3*m2 - 4*
     & k1(vb)*q1tL*q2tU*t2*P4*m2 - 4*k1(vb)*q1tL*q2tU*P3*m1*q1k2 - 
     & 4*k1(vb)*q1tL*q2tU*P3*m2*q1q2 + 4*k1(vb)*q1tL*q2tU*P3*m2*
     & q1k1 - 4*k1(vb)*q1tL*q2tU*P4*m1*q2k2 - 4*k1(vb)*q1tL*q2tU*P4
     & *m2*q1q2 - 4*k1(vb)*q1tL*q2tU*P4*m2*q2k1 + 8*k1(vb)*q1tL*
     & q2tU*P4*m2*q2k2 - 4*k1(vb)*q1tL*q2tU*P7*m2*q1q2 + 4*k2(vb)*
     & q1L*q2U*q1tq2t*P3*m1 + 4*k2(vb)*q1L*q2U*q1tq2t*P4*m1 + 4*
     & k2(vb)*q1L*q2U*q1tq2t*P7*m1
      A51(4,vb) = A51(4,vb) + 8*k2(vb)*q1L*q1tq2t*q1tU*P3*m1 - 4*
     & k2(vb)*q1L*q1tq2t*q2tU*P3*m1 - 4*k2(vb)*q1L*q1tq2t*q2tU*P4
     & *m1 - 4*k2(vb)*q1L*q1tq2t*q2tU*P7*m1 - 4*k2(vb)*q1L*q1tk1t*
     & q2tU*P3*m2 - 4*k2(vb)*q1L*q1tk2t*q2tU*P3*m1 - 8*k2(vb)*q1L*
     & q1tU*q2tk1t*P3*m1 + 4*k2(vb)*q1L*q1tU*q2tk1t*P3*m2 + 4*
     & k2(vb)*q1L*q1tU*q2tk2t*P3*m1 + 8*k2(vb)*q1L*q1tU*q2tk2t*P4
     & *m1 + 8*k2(vb)*q1L*q1tU*q2tk2t*P7*m1 - 4*k2(vb)*q1L*q1tU*t2
     & *P3*m1 + 4*k2(vb)*q1L*q1tU*t2*P4*m1 + 4*k2(vb)*q1L*q1tU*t2*
     & P7*m1 + 4*k2(vb)*q2U*q2L*q1tq2t*P3*m1 + 4*k2(vb)*q2U*q2L*
     & q1tq2t*P4*m1 + 4*k2(vb)*q2U*q2L*q1tq2t*P7*m1 + 4*k2(vb)*q2U
     & *q1tq2t*q1tL*P3*m1 - 4*k2(vb)*q2U*q1tq2t*q1tL*P4*m1 - 4*
     & k2(vb)*q2U*q1tq2t*q1tL*P7*m1 - 4*k2(vb)*q2U*q1tq2t*q2tL*P3
     & *m1 + 4*k2(vb)*q2U*q1tq2t*q2tL*P4*m1 - 4*k2(vb)*q2U*q1tq2t*
     & q2tL*P7*m1 + 8*k2(vb)*q2U*q1tk1t*q2tL*P3*m1 - 4*k2(vb)*q2U*
     & q1tk1t*q2tL*P3*m2 + 4*k2(vb)*q2U*q1tk1t*q2tL*P4*m2 + 8*
     & k2(vb)*q2U*q1tk1t*q2tL*P7*m1
      A51(4,vb) = A51(4,vb) - 4*k2(vb)*q2U*q1tk2t*q2tL*P3*m1 - 4*
     & k2(vb)*q2U*q1tk2t*q2tL*P4*m1 - 8*k2(vb)*q2U*q1tL*q2tk1t*P3
     & *m1 + 4*k2(vb)*q2U*q1tL*q2tk1t*P3*m2 - 4*k2(vb)*q2U*q1tL*
     & q2tk1t*P4*m2 + 4*k2(vb)*q2U*q1tL*q2tk2t*P3*m1 + 4*k2(vb)*
     & q2U*q1tL*q2tk2t*P4*m1 + 8*k2(vb)*q2U*q1tL*q2tk2t*P7*m1 - 4
     & *k2(vb)*q2U*q1tL*t2*P3*m1 + 4*k2(vb)*q2U*q1tL*t2*P4*m1 + 4*
     & k2(vb)*q2U*q1tL*t2*P7*m1 + 4*k2(vb)*q2U*q2tL*t1*P3*m1 - 4*
     & k2(vb)*q2U*q2tL*t1*P4*m1 + 4*k2(vb)*q2U*q2tL*t1*P7*m1 + 4*
     & k2(vb)*q2L*q1tq2t*q1tU*P3*m1 + 4*k2(vb)*q2L*q1tq2t*q1tU*P4
     & *m1 + 4*k2(vb)*q2L*q1tq2t*q1tU*P7*m1 - 8*k2(vb)*q2L*q1tq2t*
     & q2tU*P4*m1 - 8*k2(vb)*q2L*q1tk1t*q2tU*P3*m1 - 4*k2(vb)*q2L*
     & q1tk1t*q2tU*P4*m2 - 8*k2(vb)*q2L*q1tk1t*q2tU*P7*m1 + 4*
     & k2(vb)*q2L*q1tk2t*q2tU*P4*m1 + 4*k2(vb)*q2L*q1tU*q2tk1t*P4
     & *m2 + 4*k2(vb)*q2L*q1tU*q2tk2t*P4*m1 - 4*k2(vb)*q2L*q2tU*t1
     & *P3*m1 + 4*k2(vb)*q2L*q2tU*t1*P4*m1 - 4*k2(vb)*q2L*q2tU*t1*
     & P7*m1
      A51(4,vb) = A51(4,vb) + 4*k2(vb)*q1tU*q2tL*t1*P3*m1 + 4*
     & k2(vb)*q1tU*q2tL*t2*P4*m1 + 4*k2(vb)*q1tU*q2tL*P3*m1*q1q2 - 
     & 8*k2(vb)*q1tU*q2tL*P3*m1*q1k1 + 4*k2(vb)*q1tU*q2tL*P3*m1*
     & q1k2 + 4*k2(vb)*q1tU*q2tL*P3*m2*q1k1 + 4*k2(vb)*q1tU*q2tL*P4
     & *m1*q1q2 - 4*k2(vb)*q1tU*q2tL*P4*m1*q2k2 + 4*k2(vb)*q1tU*
     & q2tL*P4*m2*q2k1 + 4*k2(vb)*q1tU*q2tL*P7*m1*q1q2 - 4*k2(vb)*
     & q1tL*q2tU*t1*P3*m1 - 4*k2(vb)*q1tL*q2tU*t2*P4*m1 - 4*k2(vb)*
     & q1tL*q2tU*P3*m1*q1q2 + 8*k2(vb)*q1tL*q2tU*P3*m1*q1k1 - 4*
     & k2(vb)*q1tL*q2tU*P3*m1*q1k2 - 4*k2(vb)*q1tL*q2tU*P3*m2*q1k1
     &  - 4*k2(vb)*q1tL*q2tU*P4*m1*q1q2 + 4*k2(vb)*q1tL*q2tU*P4*m1*
     & q2k2 - 4*k2(vb)*q1tL*q2tU*P4*m2*q2k1 - 4*k2(vb)*q1tL*q2tU*P7
     & *m1*q1q2 - 4*kw(vb)*q1L*q1tk1t*q2tU*P3*m2 - 4*kw(vb)*q1L*
     & q1tk2t*q2tU*P3*m1 + 4*kw(vb)*q1L*q1tU*q2tk1t*P3*m2 + 4*
     & kw(vb)*q1L*q1tU*q2tk2t*P3*m1 - 4*kw(vb)*q2U*q1tk1t*q2tL*P3
     & *m2 + 4*kw(vb)*q2U*q1tk1t*q2tL*P4*m2 - 4*kw(vb)*q2U*q1tk2t*
     & q2tL*P3*m1
      A51(4,vb) = A51(4,vb) + 4*kw(vb)*q2U*q1tk2t*q2tL*P4*m1 + 4*
     & kw(vb)*q2U*q1tL*q2tk1t*P3*m2 - 4*kw(vb)*q2U*q1tL*q2tk1t*P4
     & *m2 + 4*kw(vb)*q2U*q1tL*q2tk2t*P3*m1 - 4*kw(vb)*q2U*q1tL*
     & q2tk2t*P4*m1 - 4*kw(vb)*q2L*q1tk1t*q2tU*P4*m2 - 4*kw(vb)*
     & q2L*q1tk2t*q2tU*P4*m1 + 4*kw(vb)*q2L*q1tU*q2tk1t*P4*m2 + 4
     & *kw(vb)*q2L*q1tU*q2tk2t*P4*m1 + 4*kw(vb)*q1tU*q2tL*P3*m1*
     & q1k2 + 4*kw(vb)*q1tU*q2tL*P3*m2*q1k1 + 4*kw(vb)*q1tU*q2tL*P4
     & *m1*q2k2 + 4*kw(vb)*q1tU*q2tL*P4*m2*q2k1 - 4*kw(vb)*q1tL*
     & q2tU*P3*m1*q1k2 - 4*kw(vb)*q1tL*q2tU*P3*m2*q1k1 - 4*kw(vb)*
     & q1tL*q2tU*P4*m1*q2k2 - 4*kw(vb)*q1tL*q2tU*P4*m2*q2k1 + 4*
     & q1(vb)*q1L*q1tk1t*q2tU*P3*m2 + 4*q1(vb)*q1L*q1tk2t*q2tU*P3
     & *m1 - 4*q1(vb)*q1L*q1tU*q2tk1t*P3*m2 - 4*q1(vb)*q1L*q1tU*
     & q2tk2t*P3*m1 - 4*q1(vb)*q2U*q1tk1t*q2tL*P7*m2 - 4*q1(vb)*
     & q2U*q1tk2t*q2tL*P7*m1 + 4*q1(vb)*q2U*q1tL*q2tk1t*P7*m2 + 4
     & *q1(vb)*q2U*q1tL*q2tk2t*P7*m1 + 4*q1(vb)*q2L*q1tk1t*q2tU*
     & P3*m2
      A51(4,vb) = A51(4,vb) + 4*q1(vb)*q2L*q1tk1t*q2tU*P7*m2 + 4*
     & q1(vb)*q2L*q1tk2t*q2tU*P3*m1 + 4*q1(vb)*q2L*q1tk2t*q2tU*P7
     & *m1 - 4*q1(vb)*q2L*q1tU*q2tk1t*P3*m2 - 4*q1(vb)*q2L*q1tU*
     & q2tk1t*P7*m2 - 4*q1(vb)*q2L*q1tU*q2tk2t*P3*m1 - 4*q1(vb)*
     & q2L*q1tU*q2tk2t*P7*m1 - 4*q1(vb)*q1tU*q2tL*P3*m1*q1k2 - 4*
     & q1(vb)*q1tU*q2tL*P3*m1*q2k2 - 4*q1(vb)*q1tU*q2tL*P3*m2*q1k1
     &  - 4*q1(vb)*q1tU*q2tL*P3*m2*q2k1 - 4*q1(vb)*q1tU*q2tL*P7*m1*
     & q2k2 - 4*q1(vb)*q1tU*q2tL*P7*m2*q2k1 + 4*q1(vb)*q1tL*q2tU*P3
     & *m1*q1k2 + 4*q1(vb)*q1tL*q2tU*P3*m1*q2k2 + 4*q1(vb)*q1tL*
     & q2tU*P3*m2*q1k1 + 4*q1(vb)*q1tL*q2tU*P3*m2*q2k1 + 4*q1(vb)*
     & q1tL*q2tU*P7*m1*q2k2 + 4*q1(vb)*q1tL*q2tU*P7*m2*q2k1 + 4*
     & q2(vb)*q1L*q1tk1t*q2tU*P4*m2 + 4*q2(vb)*q1L*q1tk1t*q2tU*P7
     & *m2 + 4*q2(vb)*q1L*q1tk2t*q2tU*P4*m1 + 4*q2(vb)*q1L*q1tk2t*
     & q2tU*P7*m1 - 4*q2(vb)*q1L*q1tU*q2tk1t*P4*m2 - 4*q2(vb)*q1L*
     & q1tU*q2tk1t*P7*m2 - 4*q2(vb)*q1L*q1tU*q2tk2t*P4*m1 - 4*
     & q2(vb)*q1L*q1tU*q2tk2t*P7*m1
      A51(4,vb) = A51(4,vb) + 4*q2(vb)*q2U*q1tk1t*q2tL*P7*m2 + 4*
     & q2(vb)*q2U*q1tk2t*q2tL*P7*m1 - 4*q2(vb)*q2U*q1tL*q2tk1t*P7
     & *m2 - 4*q2(vb)*q2U*q1tL*q2tk2t*P7*m1 + 4*q2(vb)*q2L*q1tk1t*
     & q2tU*P4*m2 + 4*q2(vb)*q2L*q1tk2t*q2tU*P4*m1 - 4*q2(vb)*q2L*
     & q1tU*q2tk1t*P4*m2 - 4*q2(vb)*q2L*q1tU*q2tk2t*P4*m1 - 4*
     & q2(vb)*q1tU*q2tL*P4*m1*q1k2 - 4*q2(vb)*q1tU*q2tL*P4*m1*q2k2
     &  - 4*q2(vb)*q1tU*q2tL*P4*m2*q1k1 - 4*q2(vb)*q1tU*q2tL*P4*m2*
     & q2k1 - 4*q2(vb)*q1tU*q2tL*P7*m1*q1k2 - 4*q2(vb)*q1tU*q2tL*P7
     & *m2*q1k1 + 4*q2(vb)*q1tL*q2tU*P4*m1*q1k2 + 4*q2(vb)*q1tL*
     & q2tU*P4*m1*q2k2 + 4*q2(vb)*q1tL*q2tU*P4*m2*q1k1 + 4*q2(vb)*
     & q1tL*q2tU*P4*m2*q2k1 + 4*q2(vb)*q1tL*q2tU*P7*m1*q1k2 + 4*
     & q2(vb)*q1tL*q2tU*P7*m2*q1k1 - 4*q1t(vb)*q1L*q2U*q2tk1t*P3*
     & m2 + 4*q1t(vb)*q1L*q2U*q2tk1t*P4*m2 - 4*q1t(vb)*q1L*q2U*
     & q2tk1t*P7*m2 - 4*q1t(vb)*q1L*q2U*q2tk2t*P3*m1 + 4*q1t(vb)*
     & q1L*q2U*q2tk2t*P4*m1 - 4*q1t(vb)*q1L*q2U*q2tk2t*P7*m1 - 4*
     & q1t(vb)*q1L*q2tU*P3*m1*q2k2
      A51(4,vb) = A51(4,vb) - 4*q1t(vb)*q1L*q2tU*P3*m2*q2k1 + 4*
     & q1t(vb)*q1L*q2tU*P4*m1*q2k2 + 4*q1t(vb)*q1L*q2tU*P4*m2*q2k1
     &  - 4*q1t(vb)*q1L*q2tU*P7*m1*q2k2 - 4*q1t(vb)*q1L*q2tU*P7*m2*
     & q2k1 - 4*q1t(vb)*q2U*q2L*q2tk1t*P3*m2 + 4*q1t(vb)*q2U*q2L*
     & q2tk1t*P4*m2 - 4*q1t(vb)*q2U*q2L*q2tk1t*P7*m2 - 4*q1t(vb)*
     & q2U*q2L*q2tk2t*P3*m1 + 4*q1t(vb)*q2U*q2L*q2tk2t*P4*m1 - 4*
     & q1t(vb)*q2U*q2L*q2tk2t*P7*m1 - 4*q1t(vb)*q2U*q2tL*P3*m1*
     & q1k2 - 4*q1t(vb)*q2U*q2tL*P3*m1*q2k2 - 4*q1t(vb)*q2U*q2tL*P3
     & *m2*q1k1 - 4*q1t(vb)*q2U*q2tL*P3*m2*q2k1 + 4*q1t(vb)*q2U*
     & q2tL*P4*m1*q1k2 + 4*q1t(vb)*q2U*q2tL*P4*m1*q2k2 + 4*q1t(vb)*
     & q2U*q2tL*P4*m2*q1k1 + 4*q1t(vb)*q2U*q2tL*P4*m2*q2k1 - 4*
     & q1t(vb)*q2U*q2tL*P7*m1*q1k2 - 4*q1t(vb)*q2U*q2tL*P7*m1*q2k2
     &  - 4*q1t(vb)*q2U*q2tL*P7*m2*q1k1 - 4*q1t(vb)*q2U*q2tL*P7*m2*
     & q2k1 + 4*q1t(vb)*q2L*q2tU*P3*m1*q1k2 + 4*q1t(vb)*q2L*q2tU*P3
     & *m2*q1k1 - 4*q1t(vb)*q2L*q2tU*P4*m1*q1k2 - 4*q1t(vb)*q2L*
     & q2tU*P4*m2*q1k1
      A51(4,vb) = A51(4,vb) + 4*q1t(vb)*q2L*q2tU*P7*m1*q1k2 + 4*
     & q1t(vb)*q2L*q2tU*P7*m2*q1k1 + 4*q2t(vb)*q1L*q2U*q1tk1t*P3*
     & m2 - 4*q2t(vb)*q1L*q2U*q1tk1t*P4*m2 - 4*q2t(vb)*q1L*q2U*
     & q1tk1t*P7*m2 + 4*q2t(vb)*q1L*q2U*q1tk2t*P3*m1 - 4*q2t(vb)*
     & q1L*q2U*q1tk2t*P4*m1 - 4*q2t(vb)*q1L*q2U*q1tk2t*P7*m1 + 4*
     & q2t(vb)*q1L*q1tU*P3*m1*q2k2 + 4*q2t(vb)*q1L*q1tU*P3*m2*q2k1
     &  - 4*q2t(vb)*q1L*q1tU*P4*m1*q2k2 - 4*q2t(vb)*q1L*q1tU*P4*m2*
     & q2k1 - 4*q2t(vb)*q1L*q1tU*P7*m1*q2k2 - 4*q2t(vb)*q1L*q1tU*P7
     & *m2*q2k1 + 4*q2t(vb)*q2U*q2L*q1tk1t*P3*m2 - 4*q2t(vb)*q2U*
     & q2L*q1tk1t*P4*m2 - 4*q2t(vb)*q2U*q2L*q1tk1t*P7*m2 + 4*
     & q2t(vb)*q2U*q2L*q1tk2t*P3*m1 - 4*q2t(vb)*q2U*q2L*q1tk2t*P4
     & *m1 - 4*q2t(vb)*q2U*q2L*q1tk2t*P7*m1 + 4*q2t(vb)*q2U*q1tL*
     & P3*m1*q1k2 + 4*q2t(vb)*q2U*q1tL*P3*m1*q2k2 + 4*q2t(vb)*q2U*
     & q1tL*P3*m2*q1k1 + 4*q2t(vb)*q2U*q1tL*P3*m2*q2k1 - 4*q2t(vb)*
     & q2U*q1tL*P4*m1*q1k2 - 4*q2t(vb)*q2U*q1tL*P4*m1*q2k2 - 4*
     & q2t(vb)*q2U*q1tL*P4*m2*q1k1
      A51(4,vb) = A51(4,vb) - 4*q2t(vb)*q2U*q1tL*P4*m2*q2k1 - 4*
     & q2t(vb)*q2U*q1tL*P7*m1*q1k2 - 4*q2t(vb)*q2U*q1tL*P7*m1*q2k2
     &  - 4*q2t(vb)*q2U*q1tL*P7*m2*q1k1 - 4*q2t(vb)*q2U*q1tL*P7*m2*
     & q2k1 - 4*q2t(vb)*q2L*q1tU*P3*m1*q1k2 - 4*q2t(vb)*q2L*q1tU*P3
     & *m2*q1k1 + 4*q2t(vb)*q2L*q1tU*P4*m1*q1k2 + 4*q2t(vb)*q2L*
     & q1tU*P4*m2*q1k1 + 4*q2t(vb)*q2L*q1tU*P7*m1*q1k2 + 4*q2t(vb)*
     & q2L*q1tU*P7*m2*q1k1 + 8*U(vb)*q1L*q1tq2t*q1tk1t*P3*m2 + 8*
     & U(vb)*q1L*q1tq2t*q1tk2t*P3*m1 - 4*U(vb)*q1L*q1tq2t*q2tk1t*
     & P3*m2 - 4*U(vb)*q1L*q1tq2t*q2tk1t*P4*m2 - 4*U(vb)*q1L*
     & q1tq2t*q2tk1t*P7*m2 - 4*U(vb)*q1L*q1tq2t*q2tk2t*P3*m1 - 4*
     & U(vb)*q1L*q1tq2t*q2tk2t*P4*m1 - 4*U(vb)*q1L*q1tq2t*q2tk2t*
     & P7*m1 - 4*U(vb)*q1L*q1tq2t*P3*m1*q2k2 - 4*U(vb)*q1L*q1tq2t*
     & P3*m2*q2k1 - 4*U(vb)*q1L*q1tq2t*P4*m1*q2k2 - 4*U(vb)*q1L*
     & q1tq2t*P4*m2*q2k1 - 4*U(vb)*q1L*q1tq2t*P7*m1*q2k2 - 4*U(vb)*
     & q1L*q1tq2t*P7*m2*q2k1 - 8*U(vb)*q1L*q1tk1t*q2tk1t*P3*m2 + 8
     & *U(vb)*q1L*q1tk1t*q2tk2t*P4*m2
      A51(4,vb) = A51(4,vb) + 8*U(vb)*q1L*q1tk1t*q2tk2t*P7*m2 - 4*
     & U(vb)*q1L*q1tk1t*t2*P3*m2 + 4*U(vb)*q1L*q1tk1t*t2*P4*m2 + 4*
     & U(vb)*q1L*q1tk1t*t2*P7*m2 - 8*U(vb)*q1L*q1tk2t*q2tk1t*P3*m1
     &  + 8*U(vb)*q1L*q1tk2t*q2tk2t*P4*m1 + 8*U(vb)*q1L*q1tk2t*
     & q2tk2t*P7*m1 - 4*U(vb)*q1L*q1tk2t*t2*P3*m1 + 4*U(vb)*q1L*
     & q1tk2t*t2*P4*m1 + 4*U(vb)*q1L*q1tk2t*t2*P7*m1 + 4*U(vb)*q2L*
     & q1tq2t*q1tk1t*P3*m2 + 4*U(vb)*q2L*q1tq2t*q1tk1t*P4*m2 + 4*
     & U(vb)*q2L*q1tq2t*q1tk1t*P7*m2 + 4*U(vb)*q2L*q1tq2t*q1tk2t*
     & P3*m1 + 4*U(vb)*q2L*q1tq2t*q1tk2t*P4*m1 + 4*U(vb)*q2L*
     & q1tq2t*q1tk2t*P7*m1 - 8*U(vb)*q2L*q1tq2t*q2tk1t*P4*m2 - 8*
     & U(vb)*q2L*q1tq2t*q2tk2t*P4*m1 + 4*U(vb)*q2L*q1tq2t*P3*m1*
     & q1k2 + 4*U(vb)*q2L*q1tq2t*P3*m2*q1k1 + 4*U(vb)*q2L*q1tq2t*P4
     & *m1*q1k2 + 4*U(vb)*q2L*q1tq2t*P4*m2*q1k1 + 4*U(vb)*q2L*
     & q1tq2t*P7*m1*q1k2 + 4*U(vb)*q2L*q1tq2t*P7*m2*q1k1 - 8*U(vb)*
     & q2L*q1tk1t*q2tk1t*P3*m2 - 8*U(vb)*q2L*q1tk1t*q2tk1t*P7*m2
     &  - 8*U(vb)*q2L*q1tk1t*q2tk2t*P3*m1
      A51(4,vb) = A51(4,vb) - 8*U(vb)*q2L*q1tk1t*q2tk2t*P7*m1 + 8*
     & U(vb)*q2L*q1tk2t*q2tk1t*P4*m2 + 8*U(vb)*q2L*q1tk2t*q2tk2t*
     & P4*m1 - 4*U(vb)*q2L*q2tk1t*t1*P3*m2 + 4*U(vb)*q2L*q2tk1t*t1*
     & P4*m2 - 4*U(vb)*q2L*q2tk1t*t1*P7*m2 - 4*U(vb)*q2L*q2tk2t*t1*
     & P3*m1 + 4*U(vb)*q2L*q2tk2t*t1*P4*m1 - 4*U(vb)*q2L*q2tk2t*t1*
     & P7*m1 + 8*U(vb)*q1tq2t*q1tL*P3*m1*q1k2 + 4*U(vb)*q1tq2t*q1tL
     & *P3*m1*q2k2 + 8*U(vb)*q1tq2t*q1tL*P3*m2*q1k1 + 4*U(vb)*q1tq2t
     & *q1tL*P3*m2*q2k1 + 4*U(vb)*q1tq2t*q1tL*P4*m1*q2k2 + 4*U(vb)*
     & q1tq2t*q1tL*P4*m2*q2k1 + 4*U(vb)*q1tq2t*q1tL*P7*m1*q2k2 + 4*
     & U(vb)*q1tq2t*q1tL*P7*m2*q2k1 - 4*U(vb)*q1tq2t*q2tL*P3*m1*
     & q1k2 - 4*U(vb)*q1tq2t*q2tL*P3*m2*q1k1 - 4*U(vb)*q1tq2t*q2tL*
     & P4*m1*q1k2 - 8*U(vb)*q1tq2t*q2tL*P4*m1*q2k2 - 4*U(vb)*q1tq2t*
     & q2tL*P4*m2*q1k1 - 8*U(vb)*q1tq2t*q2tL*P4*m2*q2k1 - 4*U(vb)*
     & q1tq2t*q2tL*P7*m1*q1k2 - 4*U(vb)*q1tq2t*q2tL*P7*m2*q1k1 + 4*
     & U(vb)*q1tk1t*q2tL*t1*P3*m2 + 4*U(vb)*q1tk1t*q2tL*t2*P4*m2 - 
     & 8*U(vb)*q1tk1t*q2tL*P3*m1*q2k2
      A51(4,vb) = A51(4,vb) + 4*U(vb)*q1tk1t*q2tL*P3*m2*q1q2 - 8*
     & U(vb)*q1tk1t*q2tL*P3*m2*q1k1 - 8*U(vb)*q1tk1t*q2tL*P3*m2*
     & q2k1 + 4*U(vb)*q1tk1t*q2tL*P4*m2*q1q2 - 8*U(vb)*q1tk1t*q2tL*
     & P4*m2*q2k2 - 8*U(vb)*q1tk1t*q2tL*P7*m1*q2k2 + 4*U(vb)*q1tk1t*
     & q2tL*P7*m2*q1q2 - 8*U(vb)*q1tk1t*q2tL*P7*m2*q2k1 + 4*U(vb)*
     & q1tk2t*q2tL*t1*P3*m1 + 4*U(vb)*q1tk2t*q2tL*t2*P4*m1 + 4*
     & U(vb)*q1tk2t*q2tL*P3*m1*q1q2 - 8*U(vb)*q1tk2t*q2tL*P3*m1*
     & q1k1 + 4*U(vb)*q1tk2t*q2tL*P4*m1*q1q2 + 8*U(vb)*q1tk2t*q2tL*
     & P4*m2*q2k1 + 4*U(vb)*q1tk2t*q2tL*P7*m1*q1q2 - 4*U(vb)*q1tL*
     & q2tk1t*t1*P3*m2 - 4*U(vb)*q1tL*q2tk1t*t2*P4*m2 - 8*U(vb)*
     & q1tL*q2tk1t*P3*m1*q1k2 - 4*U(vb)*q1tL*q2tk1t*P3*m2*q1q2 - 4*
     & U(vb)*q1tL*q2tk1t*P4*m2*q1q2 + 8*U(vb)*q1tL*q2tk1t*P4*m2*
     & q2k2 - 4*U(vb)*q1tL*q2tk1t*P7*m2*q1q2 - 4*U(vb)*q1tL*q2tk2t*
     & t1*P3*m1 - 4*U(vb)*q1tL*q2tk2t*t2*P4*m1 - 4*U(vb)*q1tL*
     & q2tk2t*P3*m1*q1q2 + 8*U(vb)*q1tL*q2tk2t*P3*m1*q1k1 - 4*U(vb)*
     & q1tL*q2tk2t*P4*m1*q1q2
      A51(4,vb) = A51(4,vb) + 8*U(vb)*q1tL*q2tk2t*P4*m1*q1k2 + 8*
     & U(vb)*q1tL*q2tk2t*P4*m1*q2k2 + 8*U(vb)*q1tL*q2tk2t*P4*m2*
     & q1k1 - 4*U(vb)*q1tL*q2tk2t*P7*m1*q1q2 + 8*U(vb)*q1tL*q2tk2t*
     & P7*m1*q1k2 + 8*U(vb)*q1tL*q2tk2t*P7*m2*q1k1 - 4*U(vb)*q1tL*t2
     & *P3*m1*q1k2 - 4*U(vb)*q1tL*t2*P3*m2*q1k1 + 4*U(vb)*q1tL*t2*P4*
     & m1*q1k2 + 4*U(vb)*q1tL*t2*P4*m2*q1k1 + 4*U(vb)*q1tL*t2*P7*m1*
     & q1k2 + 4*U(vb)*q1tL*t2*P7*m2*q1k1 - 4*U(vb)*q2tL*t1*P3*m1*q2k2
     &  - 4*U(vb)*q2tL*t1*P3*m2*q2k1 + 4*U(vb)*q2tL*t1*P4*m1*q2k2 + 4
     & *U(vb)*q2tL*t1*P4*m2*q2k1 - 4*U(vb)*q2tL*t1*P7*m1*q2k2 - 4*
     & U(vb)*q2tL*t1*P7*m2*q2k1 + 4*L(vb)*q2U*q1tq2t*q1tk1t*P3*m2
     &  - 4*L(vb)*q2U*q1tq2t*q1tk1t*P4*m2 - 4*L(vb)*q2U*q1tq2t*
     & q1tk1t*P7*m2 + 4*L(vb)*q2U*q1tq2t*q1tk2t*P3*m1 - 4*L(vb)*
     & q2U*q1tq2t*q1tk2t*P4*m1 - 4*L(vb)*q2U*q1tq2t*q1tk2t*P7*m1
     &  - 4*L(vb)*q2U*q1tq2t*q2tk1t*P3*m2 + 4*L(vb)*q2U*q1tq2t*
     & q2tk1t*P4*m2 - 4*L(vb)*q2U*q1tq2t*q2tk1t*P7*m2 - 4*L(vb)*
     & q2U*q1tq2t*q2tk2t*P3*m1
      A51(4,vb) = A51(4,vb) + 4*L(vb)*q2U*q1tq2t*q2tk2t*P4*m1 - 4*
     & L(vb)*q2U*q1tq2t*q2tk2t*P7*m1 - 4*L(vb)*q2U*q1tq2t*P3*m1*
     & q1k2 - 4*L(vb)*q2U*q1tq2t*P3*m1*q2k2 - 4*L(vb)*q2U*q1tq2t*P3
     & *m2*q1k1 - 4*L(vb)*q2U*q1tq2t*P3*m2*q2k1 - 4*L(vb)*q2U*
     & q1tq2t*P4*m1*q1k2 - 4*L(vb)*q2U*q1tq2t*P4*m1*q2k2 - 4*L(vb)*
     & q2U*q1tq2t*P4*m2*q1k1 - 4*L(vb)*q2U*q1tq2t*P4*m2*q2k1 - 4*
     & L(vb)*q2U*q1tq2t*P7*m1*q1k2 - 4*L(vb)*q2U*q1tq2t*P7*m1*q2k2
     &  - 4*L(vb)*q2U*q1tq2t*P7*m2*q1k1 - 4*L(vb)*q2U*q1tq2t*P7*m2*
     & q2k1 + 8*L(vb)*q2U*q1tk1t*q2tk1t*P7*m2 + 8*L(vb)*q2U*q1tk1t
     & *q2tk2t*P3*m1 + 8*L(vb)*q2U*q1tk1t*q2tk2t*P4*m2 + 8*L(vb)*
     & q2U*q1tk1t*q2tk2t*P7*m1 + 8*L(vb)*q2U*q1tk1t*q2tk2t*P7*m2
     &  - 4*L(vb)*q2U*q1tk1t*t2*P3*m2 + 4*L(vb)*q2U*q1tk1t*t2*P4*m2
     &  + 4*L(vb)*q2U*q1tk1t*t2*P7*m2 - 8*L(vb)*q2U*q1tk2t*q2tk1t*
     & P3*m1 - 8*L(vb)*q2U*q1tk2t*q2tk1t*P4*m2 + 8*L(vb)*q2U*
     & q1tk2t*q2tk2t*P7*m1 - 4*L(vb)*q2U*q1tk2t*t2*P3*m1 + 4*L(vb)*
     & q2U*q1tk2t*t2*P4*m1
      A51(4,vb) = A51(4,vb) + 4*L(vb)*q2U*q1tk2t*t2*P7*m1 + 4*L(vb)
     & *q2U*q2tk1t*t1*P3*m2 - 4*L(vb)*q2U*q2tk1t*t1*P4*m2 + 4*L(vb)
     & *q2U*q2tk1t*t1*P7*m2 + 4*L(vb)*q2U*q2tk2t*t1*P3*m1 - 4*L(vb)
     & *q2U*q2tk2t*t1*P4*m1 + 4*L(vb)*q2U*q2tk2t*t1*P7*m1 - 8*L(vb)
     & *q1tq2t*q1tU*P3*m1*q1k2 - 4*L(vb)*q1tq2t*q1tU*P3*m1*q2k2 - 8
     & *L(vb)*q1tq2t*q1tU*P3*m2*q1k1 - 4*L(vb)*q1tq2t*q1tU*P3*m2*
     & q2k1 - 4*L(vb)*q1tq2t*q1tU*P4*m1*q2k2 - 4*L(vb)*q1tq2t*q1tU*
     & P4*m2*q2k1 - 4*L(vb)*q1tq2t*q1tU*P7*m1*q2k2 - 4*L(vb)*q1tq2t*
     & q1tU*P7*m2*q2k1 + 4*L(vb)*q1tq2t*q2tU*P3*m1*q1k2 + 4*L(vb)*
     & q1tq2t*q2tU*P3*m2*q1k1 + 4*L(vb)*q1tq2t*q2tU*P4*m1*q1k2 + 8*
     & L(vb)*q1tq2t*q2tU*P4*m1*q2k2 + 4*L(vb)*q1tq2t*q2tU*P4*m2*
     & q1k1 + 8*L(vb)*q1tq2t*q2tU*P4*m2*q2k1 + 4*L(vb)*q1tq2t*q2tU*
     & P7*m1*q1k2 + 4*L(vb)*q1tq2t*q2tU*P7*m2*q1k1 - 4*L(vb)*q1tk1t*
     & q2tU*t1*P3*m2 - 4*L(vb)*q1tk1t*q2tU*t2*P4*m2 + 8*L(vb)*
     & q1tk1t*q2tU*P3*m1*q2k2 - 4*L(vb)*q1tk1t*q2tU*P3*m2*q1q2 + 8*
     & L(vb)*q1tk1t*q2tU*P3*m2*q1k1
      A51(4,vb) = A51(4,vb) + 8*L(vb)*q1tk1t*q2tU*P3*m2*q2k1 - 4*
     & L(vb)*q1tk1t*q2tU*P4*m2*q1q2 + 8*L(vb)*q1tk1t*q2tU*P4*m2*
     & q2k2 + 8*L(vb)*q1tk1t*q2tU*P7*m1*q2k2 - 4*L(vb)*q1tk1t*q2tU*
     & P7*m2*q1q2 + 8*L(vb)*q1tk1t*q2tU*P7*m2*q2k1 - 4*L(vb)*q1tk2t*
     & q2tU*t1*P3*m1 - 4*L(vb)*q1tk2t*q2tU*t2*P4*m1 - 4*L(vb)*
     & q1tk2t*q2tU*P3*m1*q1q2 + 8*L(vb)*q1tk2t*q2tU*P3*m1*q1k1 - 4*
     & L(vb)*q1tk2t*q2tU*P4*m1*q1q2 - 8*L(vb)*q1tk2t*q2tU*P4*m2*
     & q2k1 - 4*L(vb)*q1tk2t*q2tU*P7*m1*q1q2 + 4*L(vb)*q1tU*q2tk1t*
     & t1*P3*m2 + 4*L(vb)*q1tU*q2tk1t*t2*P4*m2 + 8*L(vb)*q1tU*
     & q2tk1t*P3*m1*q1k2 + 4*L(vb)*q1tU*q2tk1t*P3*m2*q1q2 + 4*L(vb)*
     & q1tU*q2tk1t*P4*m2*q1q2 - 8*L(vb)*q1tU*q2tk1t*P4*m2*q2k2 + 4*
     & L(vb)*q1tU*q2tk1t*P7*m2*q1q2 + 4*L(vb)*q1tU*q2tk2t*t1*P3*m1
     &  + 4*L(vb)*q1tU*q2tk2t*t2*P4*m1 + 4*L(vb)*q1tU*q2tk2t*P3*m1*
     & q1q2 - 8*L(vb)*q1tU*q2tk2t*P3*m1*q1k1 + 4*L(vb)*q1tU*q2tk2t*
     & P4*m1*q1q2 - 8*L(vb)*q1tU*q2tk2t*P4*m1*q1k2 - 8*L(vb)*q1tU*
     & q2tk2t*P4*m1*q2k2
      A51(4,vb) = A51(4,vb) - 8*L(vb)*q1tU*q2tk2t*P4*m2*q1k1 + 4*
     & L(vb)*q1tU*q2tk2t*P7*m1*q1q2 - 8*L(vb)*q1tU*q2tk2t*P7*m1*
     & q1k2 - 8*L(vb)*q1tU*q2tk2t*P7*m2*q1k1 + 4*L(vb)*q1tU*t2*P3*m1
     & *q1k2 + 4*L(vb)*q1tU*t2*P3*m2*q1k1 - 4*L(vb)*q1tU*t2*P4*m1*
     & q1k2 - 4*L(vb)*q1tU*t2*P4*m2*q1k1 - 4*L(vb)*q1tU*t2*P7*m1*q1k2
     &  - 4*L(vb)*q1tU*t2*P7*m2*q1k1 + 4*L(vb)*q2tU*t1*P3*m1*q2k2 + 4
     & *L(vb)*q2tU*t1*P3*m2*q2k1 - 4*L(vb)*q2tU*t1*P4*m1*q2k2 - 4*
     & L(vb)*q2tU*t1*P4*m2*q2k1 + 4*L(vb)*q2tU*t1*P7*m1*q2k2 + 4*
     & L(vb)*q2tU*t1*P7*m2*q2k1

      A52(1,vb) =
     & 4*k1(vb)*t1*t2*P5*m2 + 4*k1(vb)*t1*t2*P6*m2 - 4*k1(vb)*t1*t2*P8*
     & m2 + 4*k1(vb)*q1tq2t*q1tk1t*P6*m2 + 4*k1(vb)*q1tq2t*q1tk2t*
     & P6*m1 - 16*k1(vb)*q1tq2t*q1tk2t*P6*m2 - 12*k1(vb)*q1tq2t*
     & q2tk1t*P5*m2 + 4*k1(vb)*q1tq2t*q2tk2t*P5*m1 - 4*k1(vb)*
     & q1tq2t*t1*P6*m2 - 4*k1(vb)*q1tq2t*t2*P5*m2 + 4*k1(vb)*q1tq2t*
     & P5*m1*q2k2 + 4*k1(vb)*q1tq2t*P5*m2*q1q2 - 4*k1(vb)*q1tq2t*P5*
     & m2*q2k1 + 4*k1(vb)*q1tq2t*P6*m1*q1k2 + 4*k1(vb)*q1tq2t*P6*m2*
     & q1q2 + 4*k1(vb)*q1tq2t*P6*m2*q1k1 - 8*k1(vb)*q1tq2t*P6*m2*q1k2
     &  + 4*k1(vb)*q1tq2t*P8*m2*q1q2 + 4*k1(vb)*q1tq2t**2*P5*m2 + 4*
     & k1(vb)*q1tq2t**2*P6*m2 + 4*k1(vb)*q1tq2t**2*P8*m2 + 8*k1(vb)*
     & q1tk1t*q2tk1t*P5*m2 + 4*k1(vb)*q1tk1t*t2*P5*m2 - 8*k1(vb)*
     & q1tk2t*q2tk1t*P5*m1 - 8*k1(vb)*q1tk2t*q2tk1t*P6*m2 - 16*
     & k1(vb)*q1tk2t*q2tk1t*P8*m2 - 8*k1(vb)*q1tk2t*q2tk2t*P6*m1 + 
     & 16*k1(vb)*q1tk2t*q2tk2t*P6*m2 - 4*k1(vb)*q1tk2t*t2*P5*m1 + 8*
     & k1(vb)*q1tk2t*t2*P6*m2
      A52(1,vb) = A52(1,vb) - 8*k1(vb)*q1tk2t*t2*P8*m2 + 8*k1(vb)*
     & q2tk1t*t1*P5*m2 - 4*k1(vb)*q2tk1t*t1*P6*m2 - 8*k1(vb)*q2tk1t*
     & t1*P8*m2 - 4*k1(vb)*q2tk2t*t1*P6*m1 + 8*k1(vb)*q2tk2t*t1*P6*m2
     &  + 4*k2(vb)*t1*t2*P5*m1 + 4*k2(vb)*t1*t2*P6*m1 - 4*k2(vb)*t1*t2*
     & P8*m1 + 4*k2(vb)*q1tq2t*q1tk1t*P6*m2 - 12*k2(vb)*q1tq2t*
     & q1tk2t*P6*m1 - 16*k2(vb)*q1tq2t*q2tk1t*P5*m1 + 4*k2(vb)*
     & q1tq2t*q2tk1t*P5*m2 + 4*k2(vb)*q1tq2t*q2tk2t*P5*m1 - 4*
     & k2(vb)*q1tq2t*t1*P6*m1 - 4*k2(vb)*q1tq2t*t2*P5*m1 + 4*k2(vb)*
     & q1tq2t*P5*m1*q1q2 - 8*k2(vb)*q1tq2t*P5*m1*q2k1 + 4*k2(vb)*
     & q1tq2t*P5*m1*q2k2 + 4*k2(vb)*q1tq2t*P5*m2*q2k1 + 4*k2(vb)*
     & q1tq2t*P6*m1*q1q2 - 4*k2(vb)*q1tq2t*P6*m1*q1k2 + 4*k2(vb)*
     & q1tq2t*P6*m2*q1k1 + 4*k2(vb)*q1tq2t*P8*m1*q1q2 + 4*k2(vb)*
     & q1tq2t**2*P5*m1 + 4*k2(vb)*q1tq2t**2*P6*m1 + 4*k2(vb)*
     & q1tq2t**2*P8*m1 + 16*k2(vb)*q1tk1t*q2tk1t*P5*m1 - 8*k2(vb)*
     & q1tk1t*q2tk1t*P5*m2 + 8*k2(vb)*q1tk1t*t2*P5*m1 - 4*k2(vb)*
     & q1tk1t*t2*P5*m2
      A52(1,vb) = A52(1,vb) - 8*k2(vb)*q1tk2t*q2tk1t*P5*m1 - 8*
     & k2(vb)*q1tk2t*q2tk1t*P6*m2 - 16*k2(vb)*q1tk2t*q2tk1t*P8*m1
     &  + 8*k2(vb)*q1tk2t*q2tk2t*P6*m1 - 4*k2(vb)*q1tk2t*t2*P5*m1 + 
     & 8*k2(vb)*q1tk2t*t2*P6*m1 - 8*k2(vb)*q1tk2t*t2*P8*m1 + 8*k2(vb)
     & *q2tk1t*t1*P5*m1 - 4*k2(vb)*q2tk1t*t1*P6*m2 - 8*k2(vb)*q2tk1t
     & *t1*P8*m1 + 4*k2(vb)*q2tk2t*t1*P6*m1 + 4*kw(vb)*q1tq2t*q1tk1t
     & *P6*m2 + 4*kw(vb)*q1tq2t*q1tk2t*P6*m1 + 4*kw(vb)*q1tq2t*
     & q2tk1t*P5*m2 + 4*kw(vb)*q1tq2t*q2tk2t*P5*m1 + 4*kw(vb)*
     & q1tq2t*P5*m1*q2k2 + 4*kw(vb)*q1tq2t*P5*m2*q2k1 + 4*kw(vb)*
     & q1tq2t*P6*m1*q1k2 + 4*kw(vb)*q1tq2t*P6*m2*q1k1 - 8*kw(vb)*
     & q1tk1t*q2tk1t*P5*m2 - 4*kw(vb)*q1tk1t*t2*P5*m2 - 8*kw(vb)*
     & q1tk2t*q2tk1t*P5*m1 - 8*kw(vb)*q1tk2t*q2tk1t*P6*m2 - 8*
     & kw(vb)*q1tk2t*q2tk2t*P6*m1 - 4*kw(vb)*q1tk2t*t2*P5*m1 - 4*
     & kw(vb)*q2tk1t*t1*P6*m2 - 4*kw(vb)*q2tk2t*t1*P6*m1 + 4*q1(vb)*
     & q1tq2t*q1tk1t*P6*m2 + 4*q1(vb)*q1tq2t*q1tk2t*P6*m1 - 4*
     & q1(vb)*q1tq2t*q2tk1t*P6*m2
      A52(1,vb) = A52(1,vb) - 4*q1(vb)*q1tq2t*q2tk1t*P8*m2 - 4*
     & q1(vb)*q1tq2t*q2tk2t*P6*m1 - 4*q1(vb)*q1tq2t*q2tk2t*P8*m1 - 
     & 4*q1(vb)*q1tq2t*P6*m1*q1k2 - 4*q1(vb)*q1tq2t*P6*m1*q2k2 - 4*
     & q1(vb)*q1tq2t*P6*m2*q1k1 - 4*q1(vb)*q1tq2t*P6*m2*q2k1 - 4*
     & q1(vb)*q1tq2t*P8*m1*q2k2 - 4*q1(vb)*q1tq2t*P8*m2*q2k1 + 8*
     & q1(vb)*q1tk1t*q2tk1t*P8*m2 - 8*q1(vb)*q1tk1t*q2tk2t*P6*m2 - 
     & 4*q1(vb)*q1tk1t*t2*P6*m2 + 4*q1(vb)*q1tk1t*t2*P8*m2 + 8*q1(vb)
     & *q1tk2t*q2tk1t*P6*m2 + 8*q1(vb)*q1tk2t*q2tk1t*P8*m1 - 4*
     & q1(vb)*q1tk2t*t2*P6*m1 + 4*q1(vb)*q1tk2t*t2*P8*m1 + 4*q1(vb)*
     & q2tk1t*t1*P6*m2 + 4*q1(vb)*q2tk2t*t1*P6*m1 - 4*q2(vb)*q1tq2t*
     & q1tk1t*P5*m2 - 4*q2(vb)*q1tq2t*q1tk1t*P8*m2 - 4*q2(vb)*
     & q1tq2t*q1tk2t*P5*m1 - 4*q2(vb)*q1tq2t*q1tk2t*P8*m1 + 4*
     & q2(vb)*q1tq2t*q2tk1t*P5*m2 + 4*q2(vb)*q1tq2t*q2tk2t*P5*m1 - 
     & 4*q2(vb)*q1tq2t*P5*m1*q1k2 - 4*q2(vb)*q1tq2t*P5*m1*q2k2 - 4*
     & q2(vb)*q1tq2t*P5*m2*q1k1 - 4*q2(vb)*q1tq2t*P5*m2*q2k1 - 4*
     & q2(vb)*q1tq2t*P8*m1*q1k2
      A52(1,vb) = A52(1,vb) - 4*q2(vb)*q1tq2t*P8*m2*q1k1 - 8*q2(vb)*
     & q1tk1t*q2tk2t*P5*m1 + 4*q2(vb)*q1tk1t*t2*P5*m2 + 8*q2(vb)*
     & q1tk2t*q2tk1t*P5*m1 + 8*q2(vb)*q1tk2t*q2tk1t*P8*m2 + 8*
     & q2(vb)*q1tk2t*q2tk2t*P8*m1 + 4*q2(vb)*q1tk2t*t2*P5*m1 - 4*
     & q2(vb)*q2tk1t*t1*P5*m2 + 4*q2(vb)*q2tk1t*t1*P8*m2 - 4*q2(vb)*
     & q2tk2t*t1*P5*m1 + 4*q2(vb)*q2tk2t*t1*P8*m1 - 4*q1t(vb)*t2*P5*
     & m1*q1k2 - 4*q1t(vb)*t2*P5*m2*q1k1 - 4*q1t(vb)*t2*P6*m1*q1k2 - 4*
     & q1t(vb)*t2*P6*m2*q1k1 + 4*q1t(vb)*t2*P8*m1*q1k2 + 4*q1t(vb)*t2*
     & P8*m2*q1k1 - 4*q1t(vb)*q1tq2t*P5*m1*q2k2 - 4*q1t(vb)*q1tq2t*P5
     & *m2*q2k1 + 8*q1t(vb)*q1tq2t*P6*m1*q1k2 + 4*q1t(vb)*q1tq2t*P6*
     & m1*q2k2 + 8*q1t(vb)*q1tq2t*P6*m2*q1k1 + 4*q1t(vb)*q1tq2t*P6*m2
     & *q2k1 + 4*q1t(vb)*q1tq2t*P8*m1*q2k2 + 4*q1t(vb)*q1tq2t*P8*m2*
     & q2k1 - 4*q1t(vb)*q2tk1t*t1*P6*m2 + 4*q1t(vb)*q2tk1t*t2*P5*m2
     &  - 8*q1t(vb)*q2tk1t*P5*m1*q1k2 + 4*q1t(vb)*q2tk1t*P5*m2*q1q2
     &  - 8*q1t(vb)*q2tk1t*P5*m2*q1k1 - 8*q1t(vb)*q2tk1t*P5*m2*q2k1
     &  - 4*q1t(vb)*q2tk1t*P6*m2*q1q2
      A52(1,vb) = A52(1,vb) + 8*q1t(vb)*q2tk1t*P6*m2*q1k2 + 8*
     & q1t(vb)*q2tk1t*P8*m1*q1k2 - 4*q1t(vb)*q2tk1t*P8*m2*q1q2 + 8*
     & q1t(vb)*q2tk1t*P8*m2*q1k1 - 4*q1t(vb)*q2tk2t*t1*P6*m1 + 4*
     & q1t(vb)*q2tk2t*t2*P5*m1 + 4*q1t(vb)*q2tk2t*P5*m1*q1q2 - 8*
     & q1t(vb)*q2tk2t*P5*m1*q2k1 - 4*q1t(vb)*q2tk2t*P6*m1*q1q2 - 8*
     & q1t(vb)*q2tk2t*P6*m2*q1k1 - 4*q1t(vb)*q2tk2t*P8*m1*q1q2 - 4*
     & q2t(vb)*t1*P5*m1*q2k2 - 4*q2t(vb)*t1*P5*m2*q2k1 - 4*q2t(vb)*t1*
     & P6*m1*q2k2 - 4*q2t(vb)*t1*P6*m2*q2k1 + 4*q2t(vb)*t1*P8*m1*q2k2
     &  + 4*q2t(vb)*t1*P8*m2*q2k1 + 4*q2t(vb)*q1tq2t*P5*m1*q1k2 + 8*
     & q2t(vb)*q1tq2t*P5*m1*q2k2 + 4*q2t(vb)*q1tq2t*P5*m2*q1k1 + 8*
     & q2t(vb)*q1tq2t*P5*m2*q2k1 - 4*q2t(vb)*q1tq2t*P6*m1*q1k2 - 4*
     & q2t(vb)*q1tq2t*P6*m2*q1k1 + 4*q2t(vb)*q1tq2t*P8*m1*q1k2 + 4*
     & q2t(vb)*q1tq2t*P8*m2*q1k1 + 4*q2t(vb)*q1tk1t*t1*P6*m2 - 4*
     & q2t(vb)*q1tk1t*t2*P5*m2 - 8*q2t(vb)*q1tk1t*P5*m1*q2k2 - 4*
     & q2t(vb)*q1tk1t*P5*m2*q1q2 + 4*q2t(vb)*q1tk1t*P6*m2*q1q2 - 8*
     & q2t(vb)*q1tk1t*P6*m2*q1k2
      A52(1,vb) = A52(1,vb) - 4*q2t(vb)*q1tk1t*P8*m2*q1q2 + 4*
     & q2t(vb)*q1tk2t*t1*P6*m1 - 4*q2t(vb)*q1tk2t*t2*P5*m1 - 4*
     & q2t(vb)*q1tk2t*P5*m1*q1q2 + 8*q2t(vb)*q1tk2t*P5*m1*q2k1 + 4*
     & q2t(vb)*q1tk2t*P6*m1*q1q2 - 8*q2t(vb)*q1tk2t*P6*m1*q1k2 - 8*
     & q2t(vb)*q1tk2t*P6*m1*q2k2 - 8*q2t(vb)*q1tk2t*P6*m2*q2k1 - 4*
     & q2t(vb)*q1tk2t*P8*m1*q1q2 + 8*q2t(vb)*q1tk2t*P8*m1*q2k2 + 8*
     & q2t(vb)*q1tk2t*P8*m2*q2k1

      A52(2,vb) =
     &  - 4*k1(vb)*q2U*q1tq2t*q1tk2t*P5 + 4*k1(vb)*q2U*q1tq2t*
     & q1tk2t*P6 - 4*k1(vb)*q2U*q1tq2t*q1tk2t*P8 + 4*k1(vb)*q2U*
     & q1tq2t*q2tk2t*P5 - 4*k1(vb)*q2U*q1tq2t*q2tk2t*P6 - 4*k1(vb)
     & *q2U*q1tq2t*q2tk2t*P8 - 4*k1(vb)*q2U*q1tq2t*P5*m1*m2 - 4*
     & k1(vb)*q2U*q1tq2t*P5*q1k2 - 4*k1(vb)*q2U*q1tq2t*P5*q2k2 + 4*
     & k1(vb)*q2U*q1tq2t*P5*k1k2 + 4*k1(vb)*q2U*q1tq2t*P6*m1*m2 - 4
     & *k1(vb)*q2U*q1tq2t*P6*q1k2 - 4*k1(vb)*q2U*q1tq2t*P6*q2k2 - 4
     & *k1(vb)*q2U*q1tq2t*P6*k1k2 - 4*k1(vb)*q2U*q1tq2t*P8*q1k2 - 4
     & *k1(vb)*q2U*q1tq2t*P8*q2k2 - 4*k1(vb)*q2U*q1tk1t*q2tk2t*P5
     &  - 4*k1(vb)*q2U*q1tk1t*q2tk2t*P6 + 4*k1(vb)*q2U*q1tk2t*
     & q2tk1t*P5 + 4*k1(vb)*q2U*q1tk2t*q2tk1t*P6 + 8*k1(vb)*q2U*
     & q1tk2t*q2tk1t*P8 + 8*k1(vb)*q2U*q1tk2t*q2tk2t*P8 + 4*k1(vb)
     & *q2U*q1tk2t*t2*P5 - 4*k1(vb)*q2U*q1tk2t*t2*P6 + 4*k1(vb)*
     & q2U*q1tk2t*t2*P8 - 4*k1(vb)*q2U*q2tk2t*t1*P5 + 4*k1(vb)*q2U
     & *q2tk2t*t1*P6
      A52(2,vb) = A52(2,vb) + 4*k1(vb)*q2U*q2tk2t*t1*P8 - 4*k1(vb)*
     & q1tq2t*q1tU*P5*q2k2 + 4*k1(vb)*q1tq2t*q1tU*P6*m1*m2 - 8*
     & k1(vb)*q1tq2t*q1tU*P6*q1k2 - 4*k1(vb)*q1tq2t*q1tU*P6*q2k2 - 
     & 4*k1(vb)*q1tq2t*q1tU*P6*k1k2 - 4*k1(vb)*q1tq2t*q1tU*P8*q2k2
     &  + 4*k1(vb)*q1tq2t*q2tU*P5*m1*m2 + 4*k1(vb)*q1tq2t*q2tU*P5*
     & q1k2 + 8*k1(vb)*q1tq2t*q2tU*P5*q2k2 - 4*k1(vb)*q1tq2t*q2tU*
     & P5*k1k2 + 4*k1(vb)*q1tq2t*q2tU*P6*q1k2 + 4*k1(vb)*q1tq2t*
     & q2tU*P8*q1k2 - 4*k1(vb)*q1tk1t*q2tU*P5*q2k2 + 4*k1(vb)*
     & q1tk1t*q2tU*P6*q1k2 - 4*k1(vb)*q1tk2t*q2tU*t1*P6 - 4*k1(vb)*
     & q1tk2t*q2tU*t2*P5 - 4*k1(vb)*q1tk2t*q2tU*P5*q1q2 + 4*k1(vb)*
     & q1tk2t*q2tU*P5*q2k1 - 8*k1(vb)*q1tk2t*q2tU*P6*m1*m2 - 4*
     & k1(vb)*q1tk2t*q2tU*P6*q1q2 - 4*k1(vb)*q1tk2t*q2tU*P6*q1k1 + 
     & 8*k1(vb)*q1tk2t*q2tU*P6*q1k2 + 8*k1(vb)*q1tk2t*q2tU*P6*q2k2
     &  + 8*k1(vb)*q1tk2t*q2tU*P6*k1k2 - 4*k1(vb)*q1tk2t*q2tU*P8*
     & q1q2 + 8*k1(vb)*q1tk2t*q2tU*P8*q2k2 - 8*k1(vb)*q1tU*q2tk1t*
     & P5*m1*m2
      A52(2,vb) = A52(2,vb) - 8*k1(vb)*q1tU*q2tk1t*P5*q1k2 - 4*
     & k1(vb)*q1tU*q2tk1t*P5*q2k2 + 8*k1(vb)*q1tU*q2tk1t*P5*k1k2 - 
     & 4*k1(vb)*q1tU*q2tk1t*P6*q1k2 - 8*k1(vb)*q1tU*q2tk1t*P8*q1k2
     &  + 4*k1(vb)*q1tU*q2tk2t*t1*P6 + 4*k1(vb)*q1tU*q2tk2t*t2*P5
     &  + 4*k1(vb)*q1tU*q2tk2t*P5*q1q2 - 4*k1(vb)*q1tU*q2tk2t*P5*
     & q2k1 + 4*k1(vb)*q1tU*q2tk2t*P6*q1q2 + 4*k1(vb)*q1tU*q2tk2t*
     & P6*q1k1 + 4*k1(vb)*q1tU*q2tk2t*P8*q1q2 - 4*k1(vb)*q1tU*t2*P5*
     & m1*m2 - 4*k1(vb)*q1tU*t2*P5*q1k2 + 4*k1(vb)*q1tU*t2*P5*k1k2 + 
     & 4*k1(vb)*q1tU*t2*P6*q1k2 - 4*k1(vb)*q1tU*t2*P8*q1k2 - 4*k1(vb)
     & *q2tU*t1*P5*q2k2 - 4*k1(vb)*q2tU*t1*P6*m1*m2 + 4*k1(vb)*q2tU*
     & t1*P6*q2k2 + 4*k1(vb)*q2tU*t1*P6*k1k2 + 4*k1(vb)*q2tU*t1*P8*
     & q2k2 + 4*k2(vb)*q2U*q1tq2t*q1tk1t*P5 - 4*k2(vb)*q2U*q1tq2t*
     & q1tk1t*P6 + 4*k2(vb)*q2U*q1tq2t*q1tk1t*P8 - 4*k2(vb)*q2U*
     & q1tq2t*q2tk1t*P5 + 4*k2(vb)*q2U*q1tq2t*q2tk1t*P6 + 4*k2(vb)
     & *q2U*q1tq2t*q2tk1t*P8 - 4*k2(vb)*q2U*q1tq2t*P5*m1*m2 + 4*
     & k2(vb)*q2U*q1tq2t*P5*q1k1
      A52(2,vb) = A52(2,vb) + 4*k2(vb)*q2U*q1tq2t*P5*q2k1 + 4*
     & k2(vb)*q2U*q1tq2t*P5*k1k2 + 4*k2(vb)*q2U*q1tq2t*P6*m1*m2 + 4
     & *k2(vb)*q2U*q1tq2t*P6*q1k1 + 4*k2(vb)*q2U*q1tq2t*P6*q2k1 - 4
     & *k2(vb)*q2U*q1tq2t*P6*k1k2 + 4*k2(vb)*q2U*q1tq2t*P8*q1k1 + 4
     & *k2(vb)*q2U*q1tq2t*P8*q2k1 - 8*k2(vb)*q2U*q1tk1t*q2tk1t*P8
     &  + 4*k2(vb)*q2U*q1tk1t*q2tk2t*P5 + 4*k2(vb)*q2U*q1tk1t*
     & q2tk2t*P6 - 4*k2(vb)*q2U*q1tk1t*t2*P5 + 4*k2(vb)*q2U*q1tk1t
     & *t2*P6 - 4*k2(vb)*q2U*q1tk1t*t2*P8 - 4*k2(vb)*q2U*q1tk2t*
     & q2tk1t*P5 - 4*k2(vb)*q2U*q1tk2t*q2tk1t*P6 - 8*k2(vb)*q2U*
     & q1tk2t*q2tk1t*P8 + 4*k2(vb)*q2U*q2tk1t*t1*P5 - 4*k2(vb)*q2U
     & *q2tk1t*t1*P6 - 4*k2(vb)*q2U*q2tk1t*t1*P8 + 4*k2(vb)*q1tq2t*
     & q1tU*P5*q2k1 + 4*k2(vb)*q1tq2t*q1tU*P6*m1*m2 + 8*k2(vb)*
     & q1tq2t*q1tU*P6*q1k1 + 4*k2(vb)*q1tq2t*q1tU*P6*q2k1 - 4*
     & k2(vb)*q1tq2t*q1tU*P6*k1k2 + 4*k2(vb)*q1tq2t*q1tU*P8*q2k1 + 
     & 4*k2(vb)*q1tq2t*q2tU*P5*m1*m2 - 4*k2(vb)*q1tq2t*q2tU*P5*q1k1
     &  - 8*k2(vb)*q1tq2t*q2tU*P5*q2k1
      A52(2,vb) = A52(2,vb) - 4*k2(vb)*q1tq2t*q2tU*P5*k1k2 - 4*
     & k2(vb)*q1tq2t*q2tU*P6*q1k1 - 4*k2(vb)*q1tq2t*q2tU*P8*q1k1 + 
     & 4*k2(vb)*q1tk1t*q2tU*t1*P6 + 4*k2(vb)*q1tk1t*q2tU*t2*P5 + 4*
     & k2(vb)*q1tk1t*q2tU*P5*q1q2 + 4*k2(vb)*q1tk1t*q2tU*P5*q2k2 + 
     & 4*k2(vb)*q1tk1t*q2tU*P6*q1q2 - 4*k2(vb)*q1tk1t*q2tU*P6*q1k2
     &  + 4*k2(vb)*q1tk1t*q2tU*P8*q1q2 - 4*k2(vb)*q1tk2t*q2tU*P5*
     & q2k1 - 8*k2(vb)*q1tk2t*q2tU*P6*m1*m2 - 4*k2(vb)*q1tk2t*q2tU*
     & P6*q1k1 - 8*k2(vb)*q1tk2t*q2tU*P6*q2k1 + 8*k2(vb)*q1tk2t*
     & q2tU*P6*k1k2 - 8*k2(vb)*q1tk2t*q2tU*P8*q2k1 - 4*k2(vb)*q1tU*
     & q2tk1t*t1*P6 - 4*k2(vb)*q1tU*q2tk1t*t2*P5 - 8*k2(vb)*q1tU*
     & q2tk1t*P5*m1*m2 - 4*k2(vb)*q1tU*q2tk1t*P5*q1q2 + 8*k2(vb)*
     & q1tU*q2tk1t*P5*q1k1 + 8*k2(vb)*q1tU*q2tk1t*P5*q2k1 - 4*
     & k2(vb)*q1tU*q2tk1t*P5*q2k2 + 8*k2(vb)*q1tU*q2tk1t*P5*k1k2 - 
     & 4*k2(vb)*q1tU*q2tk1t*P6*q1q2 + 4*k2(vb)*q1tU*q2tk1t*P6*q1k2
     &  - 4*k2(vb)*q1tU*q2tk1t*P8*q1q2 + 8*k2(vb)*q1tU*q2tk1t*P8*
     & q1k1
      A52(2,vb) = A52(2,vb) + 4*k2(vb)*q1tU*q2tk2t*P5*q2k1 - 4*
     & k2(vb)*q1tU*q2tk2t*P6*q1k1 - 4*k2(vb)*q1tU*t2*P5*m1*m2 + 4*
     & k2(vb)*q1tU*t2*P5*q1k1 + 4*k2(vb)*q1tU*t2*P5*k1k2 - 4*k2(vb)*
     & q1tU*t2*P6*q1k1 + 4*k2(vb)*q1tU*t2*P8*q1k1 + 4*k2(vb)*q2tU*t1
     & *P5*q2k1 - 4*k2(vb)*q2tU*t1*P6*m1*m2 - 4*k2(vb)*q2tU*t1*P6*
     & q2k1 + 4*k2(vb)*q2tU*t1*P6*k1k2 - 4*k2(vb)*q2tU*t1*P8*q2k1 - 4
     & *kw(vb)*q2U*q1tq2t*P5*m1*m2 + 4*kw(vb)*q2U*q1tq2t*P5*k1k2 + 
     & 4*kw(vb)*q2U*q1tq2t*P6*m1*m2 - 4*kw(vb)*q2U*q1tq2t*P6*k1k2
     &  + 4*kw(vb)*q2U*q1tk1t*q2tk2t*P5 - 4*kw(vb)*q2U*q1tk1t*
     & q2tk2t*P6 - 4*kw(vb)*q2U*q1tk2t*q2tk1t*P5 + 4*kw(vb)*q2U*
     & q1tk2t*q2tk1t*P6 + 4*kw(vb)*q1tq2t*q1tU*P6*m1*m2 - 4*kw(vb)*
     & q1tq2t*q1tU*P6*k1k2 + 4*kw(vb)*q1tq2t*q2tU*P5*m1*m2 - 4*
     & kw(vb)*q1tq2t*q2tU*P5*k1k2 + 4*kw(vb)*q1tk1t*q2tU*P5*q2k2 + 
     & 4*kw(vb)*q1tk1t*q2tU*P6*q1k2 - 4*kw(vb)*q1tk2t*q2tU*P5*q2k1
     &  - 8*kw(vb)*q1tk2t*q2tU*P6*m1*m2 - 4*kw(vb)*q1tk2t*q2tU*P6*
     & q1k1
      A52(2,vb) = A52(2,vb) + 8*kw(vb)*q1tk2t*q2tU*P6*k1k2 - 8*
     & kw(vb)*q1tU*q2tk1t*P5*m1*m2 - 4*kw(vb)*q1tU*q2tk1t*P5*q2k2
     &  + 8*kw(vb)*q1tU*q2tk1t*P5*k1k2 - 4*kw(vb)*q1tU*q2tk1t*P6*
     & q1k2 + 4*kw(vb)*q1tU*q2tk2t*P5*q2k1 + 4*kw(vb)*q1tU*q2tk2t*
     & P6*q1k1 - 4*kw(vb)*q1tU*t2*P5*m1*m2 + 4*kw(vb)*q1tU*t2*P5*k1k2
     &  - 4*kw(vb)*q2tU*t1*P6*m1*m2 + 4*kw(vb)*q2tU*t1*P6*k1k2 + 4*
     & q1(vb)*q2U*q1tq2t*P8*m1*m2 - 4*q1(vb)*q2U*q1tq2t*P8*k1k2 - 4
     & *q1(vb)*q2U*q1tk1t*q2tk2t*P8 + 4*q1(vb)*q2U*q1tk2t*q2tk1t*
     & P8 + 4*q1(vb)*q1tq2t*q1tU*P6*m1*m2 - 4*q1(vb)*q1tq2t*q1tU*P6
     & *k1k2 - 4*q1(vb)*q1tq2t*q2tU*P6*m1*m2 + 4*q1(vb)*q1tq2t*q2tU
     & *P6*k1k2 - 4*q1(vb)*q1tq2t*q2tU*P8*m1*m2 + 4*q1(vb)*q1tq2t*
     & q2tU*P8*k1k2 - 4*q1(vb)*q1tk1t*q2tU*P6*q1k2 - 4*q1(vb)*
     & q1tk1t*q2tU*P6*q2k2 - 4*q1(vb)*q1tk1t*q2tU*P8*q2k2 + 8*
     & q1(vb)*q1tk2t*q2tU*P6*m1*m2 + 4*q1(vb)*q1tk2t*q2tU*P6*q1k1
     &  + 4*q1(vb)*q1tk2t*q2tU*P6*q2k1 - 8*q1(vb)*q1tk2t*q2tU*P6*
     & k1k2
      A52(2,vb) = A52(2,vb) + 4*q1(vb)*q1tk2t*q2tU*P8*q2k1 + 4*
     & q1(vb)*q1tU*q2tk1t*P6*q1k2 + 4*q1(vb)*q1tU*q2tk1t*P6*q2k2 + 
     & 8*q1(vb)*q1tU*q2tk1t*P8*m1*m2 + 4*q1(vb)*q1tU*q2tk1t*P8*q2k2
     &  - 8*q1(vb)*q1tU*q2tk1t*P8*k1k2 - 8*q1(vb)*q1tU*q2tk2t*P6*m1
     & *m2 - 4*q1(vb)*q1tU*q2tk2t*P6*q1k1 - 4*q1(vb)*q1tU*q2tk2t*P6
     & *q2k1 + 8*q1(vb)*q1tU*q2tk2t*P6*k1k2 - 4*q1(vb)*q1tU*q2tk2t*
     & P8*q2k1 - 4*q1(vb)*q1tU*t2*P6*m1*m2 + 4*q1(vb)*q1tU*t2*P6*k1k2
     &  + 4*q1(vb)*q1tU*t2*P8*m1*m2 - 4*q1(vb)*q1tU*t2*P8*k1k2 + 4*
     & q1(vb)*q2tU*t1*P6*m1*m2 - 4*q1(vb)*q2tU*t1*P6*k1k2 - 4*q2(vb)*
     & q2U*q1tq2t*P8*m1*m2 + 4*q2(vb)*q2U*q1tq2t*P8*k1k2 + 4*q2(vb)
     & *q2U*q1tk1t*q2tk2t*P8 - 4*q2(vb)*q2U*q1tk2t*q2tk1t*P8 - 4*
     & q2(vb)*q1tq2t*q1tU*P5*m1*m2 + 4*q2(vb)*q1tq2t*q1tU*P5*k1k2
     &  - 4*q2(vb)*q1tq2t*q1tU*P8*m1*m2 + 4*q2(vb)*q1tq2t*q1tU*P8*
     & k1k2 + 4*q2(vb)*q1tq2t*q2tU*P5*m1*m2 - 4*q2(vb)*q1tq2t*q2tU*
     & P5*k1k2 - 8*q2(vb)*q1tk1t*q2tU*P5*m1*m2 - 4*q2(vb)*q1tk1t*
     & q2tU*P5*q1k2
      A52(2,vb) = A52(2,vb) - 4*q2(vb)*q1tk1t*q2tU*P5*q2k2 + 8*
     & q2(vb)*q1tk1t*q2tU*P5*k1k2 - 4*q2(vb)*q1tk1t*q2tU*P8*q1k2 + 
     & 4*q2(vb)*q1tk2t*q2tU*P5*q1k1 + 4*q2(vb)*q1tk2t*q2tU*P5*q2k1
     &  + 8*q2(vb)*q1tk2t*q2tU*P8*m1*m2 + 4*q2(vb)*q1tk2t*q2tU*P8*
     & q1k1 - 8*q2(vb)*q1tk2t*q2tU*P8*k1k2 + 8*q2(vb)*q1tU*q2tk1t*
     & P5*m1*m2 + 4*q2(vb)*q1tU*q2tk1t*P5*q1k2 + 4*q2(vb)*q1tU*
     & q2tk1t*P5*q2k2 - 8*q2(vb)*q1tU*q2tk1t*P5*k1k2 + 4*q2(vb)*
     & q1tU*q2tk1t*P8*q1k2 - 4*q2(vb)*q1tU*q2tk2t*P5*q1k1 - 4*
     & q2(vb)*q1tU*q2tk2t*P5*q2k1 - 4*q2(vb)*q1tU*q2tk2t*P8*q1k1 + 
     & 4*q2(vb)*q1tU*t2*P5*m1*m2 - 4*q2(vb)*q1tU*t2*P5*k1k2 - 4*
     & q2(vb)*q2tU*t1*P5*m1*m2 + 4*q2(vb)*q2tU*t1*P5*k1k2 + 4*q2(vb)*
     & q2tU*t1*P8*m1*m2 - 4*q2(vb)*q2tU*t1*P8*k1k2 + 4*q1t(vb)*q2U*
     & q1tq2t*P5*m1*m2 - 4*q1t(vb)*q2U*q1tq2t*P5*k1k2 + 4*q1t(vb)*
     & q2U*q1tq2t*P6*m1*m2 - 4*q1t(vb)*q2U*q1tq2t*P6*k1k2 - 4*
     & q1t(vb)*q2U*q1tq2t*P8*m1*m2 + 4*q1t(vb)*q2U*q1tq2t*P8*k1k2
     &  - 8*q1t(vb)*q2U*q2tk1t*P5*m1*m2
      A52(2,vb) = A52(2,vb) - 4*q1t(vb)*q2U*q2tk1t*P5*q1k2 - 4*
     & q1t(vb)*q2U*q2tk1t*P5*q2k2 + 8*q1t(vb)*q2U*q2tk1t*P5*k1k2 + 
     & 4*q1t(vb)*q2U*q2tk1t*P6*q1k2 + 4*q1t(vb)*q2U*q2tk1t*P6*q2k2
     &  + 8*q1t(vb)*q2U*q2tk1t*P8*m1*m2 + 4*q1t(vb)*q2U*q2tk1t*P8*
     & q1k2 + 4*q1t(vb)*q2U*q2tk1t*P8*q2k2 - 8*q1t(vb)*q2U*q2tk1t*
     & P8*k1k2 + 4*q1t(vb)*q2U*q2tk2t*P5*q1k1 + 4*q1t(vb)*q2U*
     & q2tk2t*P5*q2k1 - 8*q1t(vb)*q2U*q2tk2t*P6*m1*m2 - 4*q1t(vb)*
     & q2U*q2tk2t*P6*q1k1 - 4*q1t(vb)*q2U*q2tk2t*P6*q2k1 + 8*
     & q1t(vb)*q2U*q2tk2t*P6*k1k2 - 4*q1t(vb)*q2U*q2tk2t*P8*q1k1 - 
     & 4*q1t(vb)*q2U*q2tk2t*P8*q2k1 - 4*q1t(vb)*q2U*t2*P5*m1*m2 + 4*
     & q1t(vb)*q2U*t2*P5*k1k2 - 4*q1t(vb)*q2U*t2*P6*m1*m2 + 4*q1t(vb)
     & *q2U*t2*P6*k1k2 + 4*q1t(vb)*q2U*t2*P8*m1*m2 - 4*q1t(vb)*q2U*
     & t2*P8*k1k2 - 4*q1t(vb)*q2tU*t1*P6*m1*m2 + 4*q1t(vb)*q2tU*t1*P6
     & *k1k2 + 4*q1t(vb)*q2tU*t2*P5*m1*m2 - 4*q1t(vb)*q2tU*t2*P5*k1k2
     &  + 4*q1t(vb)*q2tU*P5*m1*m2*q1q2 - 8*q1t(vb)*q2tU*P5*m1*m2*q2k1
     &  - 4*q1t(vb)*q2tU*P5*q1q2*k1k2
      A52(2,vb) = A52(2,vb) + 4*q1t(vb)*q2tU*P5*q1k1*q2k2 - 4*
     & q1t(vb)*q2tU*P5*q1k2*q2k1 + 8*q1t(vb)*q2tU*P5*q2k1*k1k2 - 4*
     & q1t(vb)*q2tU*P6*m1*m2*q1q2 + 8*q1t(vb)*q2tU*P6*m1*m2*q1k2 + 4*
     & q1t(vb)*q2tU*P6*q1q2*k1k2 - 4*q1t(vb)*q2tU*P6*q1k1*q2k2 + 4*
     & q1t(vb)*q2tU*P6*q1k2*q2k1 - 8*q1t(vb)*q2tU*P6*q1k2*k1k2 - 4*
     & q1t(vb)*q2tU*P8*m1*m2*q1q2 + 4*q1t(vb)*q2tU*P8*q1q2*k1k2 - 4*
     & q1t(vb)*q2tU*P8*q1k1*q2k2 + 4*q1t(vb)*q2tU*P8*q1k2*q2k1 - 4*
     & q2t(vb)*q2U*q1tq2t*P5*m1*m2 + 4*q2t(vb)*q2U*q1tq2t*P5*k1k2
     &  - 4*q2t(vb)*q2U*q1tq2t*P6*m1*m2 + 4*q2t(vb)*q2U*q1tq2t*P6*
     & k1k2 + 4*q2t(vb)*q2U*q1tq2t*P8*m1*m2 - 4*q2t(vb)*q2U*q1tq2t*
     & P8*k1k2 + 8*q2t(vb)*q2U*q1tk1t*P5*m1*m2 + 4*q2t(vb)*q2U*
     & q1tk1t*P5*q1k2 + 4*q2t(vb)*q2U*q1tk1t*P5*q2k2 - 8*q2t(vb)*
     & q2U*q1tk1t*P5*k1k2 - 4*q2t(vb)*q2U*q1tk1t*P6*q1k2 - 4*
     & q2t(vb)*q2U*q1tk1t*P6*q2k2 + 4*q2t(vb)*q2U*q1tk1t*P8*q1k2 + 
     & 4*q2t(vb)*q2U*q1tk1t*P8*q2k2 - 4*q2t(vb)*q2U*q1tk2t*P5*q1k1
     &  - 4*q2t(vb)*q2U*q1tk2t*P5*q2k1
      A52(2,vb) = A52(2,vb) + 8*q2t(vb)*q2U*q1tk2t*P6*m1*m2 + 4*
     & q2t(vb)*q2U*q1tk2t*P6*q1k1 + 4*q2t(vb)*q2U*q1tk2t*P6*q2k1 - 
     & 8*q2t(vb)*q2U*q1tk2t*P6*k1k2 - 8*q2t(vb)*q2U*q1tk2t*P8*m1*m2
     &  - 4*q2t(vb)*q2U*q1tk2t*P8*q1k1 - 4*q2t(vb)*q2U*q1tk2t*P8*
     & q2k1 + 8*q2t(vb)*q2U*q1tk2t*P8*k1k2 + 4*q2t(vb)*q2U*t1*P5*m1*
     & m2 - 4*q2t(vb)*q2U*t1*P5*k1k2 + 4*q2t(vb)*q2U*t1*P6*m1*m2 - 4*
     & q2t(vb)*q2U*t1*P6*k1k2 - 4*q2t(vb)*q2U*t1*P8*m1*m2 + 4*q2t(vb)
     & *q2U*t1*P8*k1k2 + 4*q2t(vb)*q1tU*t1*P6*m1*m2 - 4*q2t(vb)*q1tU
     & *t1*P6*k1k2 - 4*q2t(vb)*q1tU*t2*P5*m1*m2 + 4*q2t(vb)*q1tU*t2*
     & P5*k1k2 - 4*q2t(vb)*q1tU*P5*m1*m2*q1q2 + 8*q2t(vb)*q1tU*P5*m1*
     & m2*q2k1 + 4*q2t(vb)*q1tU*P5*q1q2*k1k2 - 4*q2t(vb)*q1tU*P5*q1k1
     & *q2k2 + 4*q2t(vb)*q1tU*P5*q1k2*q2k1 - 8*q2t(vb)*q1tU*P5*q2k1*
     & k1k2 + 4*q2t(vb)*q1tU*P6*m1*m2*q1q2 - 8*q2t(vb)*q1tU*P6*m1*m2*
     & q1k2 - 4*q2t(vb)*q1tU*P6*q1q2*k1k2 + 4*q2t(vb)*q1tU*P6*q1k1*
     & q2k2 - 4*q2t(vb)*q1tU*P6*q1k2*q2k1 + 8*q2t(vb)*q1tU*P6*q1k2*
     & k1k2
      A52(2,vb) = A52(2,vb) - 4*q2t(vb)*q1tU*P8*m1*m2*q1q2 + 4*
     & q2t(vb)*q1tU*P8*q1q2*k1k2 - 4*q2t(vb)*q1tU*P8*q1k1*q2k2 + 4*
     & q2t(vb)*q1tU*P8*q1k2*q2k1 - 4*U(vb)*t1*t2*P5*m1*m2 + 4*U(vb)*t1
     & *t2*P5*k1k2 - 4*U(vb)*t1*t2*P6*m1*m2 + 4*U(vb)*t1*t2*P6*k1k2 + 4
     & *U(vb)*t1*t2*P8*m1*m2 - 4*U(vb)*t1*t2*P8*k1k2 - 4*U(vb)*q1tq2t*
     & q1tk1t*P5*q2k2 - 8*U(vb)*q1tq2t*q1tk1t*P6*q1k2 - 4*U(vb)*
     & q1tq2t*q1tk1t*P6*q2k2 - 4*U(vb)*q1tq2t*q1tk1t*P8*q2k2 + 4*
     & U(vb)*q1tq2t*q1tk2t*P5*q2k1 + 16*U(vb)*q1tq2t*q1tk2t*P6*m1*
     & m2 + 8*U(vb)*q1tq2t*q1tk2t*P6*q1k1 + 4*U(vb)*q1tq2t*q1tk2t*
     & P6*q2k1 - 16*U(vb)*q1tq2t*q1tk2t*P6*k1k2 + 4*U(vb)*q1tq2t*
     & q1tk2t*P8*q2k1 + 16*U(vb)*q1tq2t*q2tk1t*P5*m1*m2 + 4*U(vb)*
     & q1tq2t*q2tk1t*P5*q1k2 + 8*U(vb)*q1tq2t*q2tk1t*P5*q2k2 - 16*
     & U(vb)*q1tq2t*q2tk1t*P5*k1k2 + 4*U(vb)*q1tq2t*q2tk1t*P6*q1k2
     &  + 4*U(vb)*q1tq2t*q2tk1t*P8*q1k2 - 4*U(vb)*q1tq2t*q2tk2t*P5*
     & q1k1 - 8*U(vb)*q1tq2t*q2tk2t*P5*q2k1 - 4*U(vb)*q1tq2t*q2tk2t
     & *P6*q1k1
      A52(2,vb) = A52(2,vb) - 4*U(vb)*q1tq2t*q2tk2t*P8*q1k1 + 4*
     & U(vb)*q1tq2t*t1*P6*m1*m2 - 4*U(vb)*q1tq2t*t1*P6*k1k2 + 4*U(vb)
     & *q1tq2t*t2*P5*m1*m2 - 4*U(vb)*q1tq2t*t2*P5*k1k2 - 4*U(vb)*
     & q1tq2t*P5*m1*m2*q1q2 + 8*U(vb)*q1tq2t*P5*m1*m2*q2k1 + 4*U(vb)*
     & q1tq2t*P5*q1q2*k1k2 - 4*U(vb)*q1tq2t*P5*q1k1*q2k2 + 4*U(vb)*
     & q1tq2t*P5*q1k2*q2k1 - 8*U(vb)*q1tq2t*P5*q2k1*k1k2 - 4*U(vb)*
     & q1tq2t*P6*m1*m2*q1q2 + 8*U(vb)*q1tq2t*P6*m1*m2*q1k2 + 4*U(vb)*
     & q1tq2t*P6*q1q2*k1k2 - 4*U(vb)*q1tq2t*P6*q1k1*q2k2 + 4*U(vb)*
     & q1tq2t*P6*q1k2*q2k1 - 8*U(vb)*q1tq2t*P6*q1k2*k1k2 - 4*U(vb)*
     & q1tq2t*P8*m1*m2*q1q2 + 4*U(vb)*q1tq2t*P8*q1q2*k1k2 - 4*U(vb)*
     & q1tq2t*P8*q1k1*q2k2 + 4*U(vb)*q1tq2t*P8*q1k2*q2k1 - 4*U(vb)*
     & q1tq2t**2*P5*m1*m2 + 4*U(vb)*q1tq2t**2*P5*k1k2 - 4*U(vb)*
     & q1tq2t**2*P6*m1*m2 + 4*U(vb)*q1tq2t**2*P6*k1k2 - 4*U(vb)*
     & q1tq2t**2*P8*m1*m2 + 4*U(vb)*q1tq2t**2*P8*k1k2 - 16*U(vb)*
     & q1tk1t*q2tk1t*P5*m1*m2 - 8*U(vb)*q1tk1t*q2tk1t*P5*q1k2 - 8*
     & U(vb)*q1tk1t*q2tk1t*P5*q2k2
      A52(2,vb) = A52(2,vb) + 16*U(vb)*q1tk1t*q2tk1t*P5*k1k2 - 8*
     & U(vb)*q1tk1t*q2tk1t*P8*q1k2 + 4*U(vb)*q1tk1t*q2tk2t*t1*P6 + 
     & 4*U(vb)*q1tk1t*q2tk2t*t2*P5 + 4*U(vb)*q1tk1t*q2tk2t*P5*q1q2
     &  + 4*U(vb)*q1tk1t*q2tk2t*P6*q1q2 + 4*U(vb)*q1tk1t*q2tk2t*P8*
     & q1q2 - 8*U(vb)*q1tk1t*t2*P5*m1*m2 - 4*U(vb)*q1tk1t*t2*P5*q1k2
     &  + 8*U(vb)*q1tk1t*t2*P5*k1k2 + 4*U(vb)*q1tk1t*t2*P6*q1k2 - 4*
     & U(vb)*q1tk1t*t2*P8*q1k2 - 4*U(vb)*q1tk2t*q2tk1t*t1*P6 - 4*
     & U(vb)*q1tk2t*q2tk1t*t2*P5 - 4*U(vb)*q1tk2t*q2tk1t*P5*q1q2 + 
     & 8*U(vb)*q1tk2t*q2tk1t*P5*q1k1 + 8*U(vb)*q1tk2t*q2tk1t*P5*
     & q2k1 - 4*U(vb)*q1tk2t*q2tk1t*P6*q1q2 + 8*U(vb)*q1tk2t*q2tk1t
     & *P6*q1k2 + 8*U(vb)*q1tk2t*q2tk1t*P6*q2k2 + 16*U(vb)*q1tk2t*
     & q2tk1t*P8*m1*m2 - 4*U(vb)*q1tk2t*q2tk1t*P8*q1q2 + 8*U(vb)*
     & q1tk2t*q2tk1t*P8*q1k1 + 8*U(vb)*q1tk2t*q2tk1t*P8*q2k2 - 16*
     & U(vb)*q1tk2t*q2tk1t*P8*k1k2 - 16*U(vb)*q1tk2t*q2tk2t*P6*m1*
     & m2 - 8*U(vb)*q1tk2t*q2tk2t*P6*q1k1 - 8*U(vb)*q1tk2t*q2tk2t*
     & P6*q2k1
      A52(2,vb) = A52(2,vb) + 16*U(vb)*q1tk2t*q2tk2t*P6*k1k2 - 8*
     & U(vb)*q1tk2t*q2tk2t*P8*q2k1 + 4*U(vb)*q1tk2t*t2*P5*q1k1 - 8*
     & U(vb)*q1tk2t*t2*P6*m1*m2 - 4*U(vb)*q1tk2t*t2*P6*q1k1 + 8*U(vb)
     & *q1tk2t*t2*P6*k1k2 + 8*U(vb)*q1tk2t*t2*P8*m1*m2 + 4*U(vb)*
     & q1tk2t*t2*P8*q1k1 - 8*U(vb)*q1tk2t*t2*P8*k1k2 - 8*U(vb)*
     & q2tk1t*t1*P5*m1*m2 - 4*U(vb)*q2tk1t*t1*P5*q2k2 + 8*U(vb)*
     & q2tk1t*t1*P5*k1k2 + 4*U(vb)*q2tk1t*t1*P6*q2k2 + 8*U(vb)*
     & q2tk1t*t1*P8*m1*m2 + 4*U(vb)*q2tk1t*t1*P8*q2k2 - 8*U(vb)*
     & q2tk1t*t1*P8*k1k2 + 4*U(vb)*q2tk2t*t1*P5*q2k1 - 8*U(vb)*
     & q2tk2t*t1*P6*m1*m2 - 4*U(vb)*q2tk2t*t1*P6*q2k1 + 8*U(vb)*
     & q2tk2t*t1*P6*k1k2 - 4*U(vb)*q2tk2t*t1*P8*q2k1

      A52(3,vb) =
     &  - 8*k1(vb)*q1L*q1tq2t*q1tk2t*P6 + 4*k1(vb)*q1L*q1tq2t*
     & q2tk2t*P5 + 4*k1(vb)*q1L*q1tq2t*q2tk2t*P6 + 4*k1(vb)*q1L*
     & q1tq2t*q2tk2t*P8 + 4*k1(vb)*q1L*q1tq2t*P5*q2k2 - 4*k1(vb)*
     & q1L*q1tq2t*P6*m1*m2 + 4*k1(vb)*q1L*q1tq2t*P6*q2k2 + 4*k1(vb)
     & *q1L*q1tq2t*P6*k1k2 + 4*k1(vb)*q1L*q1tq2t*P8*q2k2 + 4*k1(vb)
     & *q1L*q1tk1t*q2tk2t*P6 - 8*k1(vb)*q1L*q1tk2t*q2tk1t*P5 - 4*
     & k1(vb)*q1L*q1tk2t*q2tk1t*P6 - 8*k1(vb)*q1L*q1tk2t*q2tk1t*
     & P8 + 8*k1(vb)*q1L*q1tk2t*q2tk2t*P6 - 4*k1(vb)*q1L*q1tk2t*t2
     & *P5 + 4*k1(vb)*q1L*q1tk2t*t2*P6 - 4*k1(vb)*q1L*q1tk2t*t2*P8
     &  - 4*k1(vb)*q2L*q1tq2t*q1tk2t*P5 - 4*k1(vb)*q2L*q1tq2t*
     & q1tk2t*P6 - 4*k1(vb)*q2L*q1tq2t*q1tk2t*P8 + 8*k1(vb)*q2L*
     & q1tq2t*q2tk2t*P5 - 4*k1(vb)*q2L*q1tq2t*P5*m1*m2 - 4*k1(vb)*
     & q2L*q1tq2t*P5*q1k2 + 4*k1(vb)*q2L*q1tq2t*P5*k1k2 - 4*k1(vb)*
     & q2L*q1tq2t*P6*q1k2 - 4*k1(vb)*q2L*q1tq2t*P8*q1k2 - 4*k1(vb)*
     & q2L*q1tk1t*q2tk2t*P5
      A52(3,vb) = A52(3,vb) - 4*k1(vb)*q2L*q1tk2t*q2tk1t*P5 + 8*
     & k1(vb)*q2L*q1tk2t*q2tk2t*P6 + 8*k1(vb)*q2L*q1tk2t*q2tk2t*
     & P8 - 4*k1(vb)*q2L*q2tk2t*t1*P5 + 4*k1(vb)*q2L*q2tk2t*t1*P6
     &  + 4*k1(vb)*q2L*q2tk2t*t1*P8 - 4*k1(vb)*q1tq2t*q1tL*P5*q2k2
     &  + 4*k1(vb)*q1tq2t*q1tL*P6*m1*m2 - 8*k1(vb)*q1tq2t*q1tL*P6*
     & q1k2 - 4*k1(vb)*q1tq2t*q1tL*P6*q2k2 - 4*k1(vb)*q1tq2t*q1tL*
     & P6*k1k2 - 4*k1(vb)*q1tq2t*q1tL*P8*q2k2 + 4*k1(vb)*q1tq2t*
     & q2tL*P5*m1*m2 + 4*k1(vb)*q1tq2t*q2tL*P5*q1k2 + 8*k1(vb)*
     & q1tq2t*q2tL*P5*q2k2 - 4*k1(vb)*q1tq2t*q2tL*P5*k1k2 + 4*
     & k1(vb)*q1tq2t*q2tL*P6*q1k2 + 4*k1(vb)*q1tq2t*q2tL*P8*q1k2 - 
     & 4*k1(vb)*q1tk1t*q2tL*P5*q2k2 + 4*k1(vb)*q1tk1t*q2tL*P6*q1k2
     &  - 4*k1(vb)*q1tk2t*q2tL*t1*P6 - 4*k1(vb)*q1tk2t*q2tL*t2*P5
     &  - 4*k1(vb)*q1tk2t*q2tL*P5*q1q2 + 4*k1(vb)*q1tk2t*q2tL*P5*
     & q2k1 - 8*k1(vb)*q1tk2t*q2tL*P6*m1*m2 - 4*k1(vb)*q1tk2t*q2tL*
     & P6*q1q2 - 4*k1(vb)*q1tk2t*q2tL*P6*q1k1 + 8*k1(vb)*q1tk2t*
     & q2tL*P6*q1k2
      A52(3,vb) = A52(3,vb) + 8*k1(vb)*q1tk2t*q2tL*P6*q2k2 + 8*
     & k1(vb)*q1tk2t*q2tL*P6*k1k2 - 4*k1(vb)*q1tk2t*q2tL*P8*q1q2 + 
     & 8*k1(vb)*q1tk2t*q2tL*P8*q2k2 - 8*k1(vb)*q1tL*q2tk1t*P5*m1*m2
     &  - 8*k1(vb)*q1tL*q2tk1t*P5*q1k2 - 4*k1(vb)*q1tL*q2tk1t*P5*
     & q2k2 + 8*k1(vb)*q1tL*q2tk1t*P5*k1k2 - 4*k1(vb)*q1tL*q2tk1t*
     & P6*q1k2 - 8*k1(vb)*q1tL*q2tk1t*P8*q1k2 + 4*k1(vb)*q1tL*
     & q2tk2t*t1*P6 + 4*k1(vb)*q1tL*q2tk2t*t2*P5 + 4*k1(vb)*q1tL*
     & q2tk2t*P5*q1q2 - 4*k1(vb)*q1tL*q2tk2t*P5*q2k1 + 4*k1(vb)*
     & q1tL*q2tk2t*P6*q1q2 + 4*k1(vb)*q1tL*q2tk2t*P6*q1k1 + 4*
     & k1(vb)*q1tL*q2tk2t*P8*q1q2 - 4*k1(vb)*q1tL*t2*P5*m1*m2 - 4*
     & k1(vb)*q1tL*t2*P5*q1k2 + 4*k1(vb)*q1tL*t2*P5*k1k2 + 4*k1(vb)*
     & q1tL*t2*P6*q1k2 - 4*k1(vb)*q1tL*t2*P8*q1k2 - 4*k1(vb)*q2tL*t1
     & *P5*q2k2 - 4*k1(vb)*q2tL*t1*P6*m1*m2 + 4*k1(vb)*q2tL*t1*P6*
     & q2k2 + 4*k1(vb)*q2tL*t1*P6*k1k2 + 4*k1(vb)*q2tL*t1*P8*q2k2 + 8
     & *k2(vb)*q1L*q1tq2t*q1tk1t*P6 - 4*k2(vb)*q1L*q1tq2t*q2tk1t*
     & P5
      A52(3,vb) = A52(3,vb) - 4*k2(vb)*q1L*q1tq2t*q2tk1t*P6 - 4*
     & k2(vb)*q1L*q1tq2t*q2tk1t*P8 - 4*k2(vb)*q1L*q1tq2t*P5*q2k1
     &  - 4*k2(vb)*q1L*q1tq2t*P6*m1*m2 - 4*k2(vb)*q1L*q1tq2t*P6*
     & q2k1 + 4*k2(vb)*q1L*q1tq2t*P6*k1k2 - 4*k2(vb)*q1L*q1tq2t*P8*
     & q2k1 + 8*k2(vb)*q1L*q1tk1t*q2tk1t*P5 + 8*k2(vb)*q1L*q1tk1t*
     & q2tk1t*P8 - 4*k2(vb)*q1L*q1tk1t*q2tk2t*P6 + 4*k2(vb)*q1L*
     & q1tk1t*t2*P5 - 4*k2(vb)*q1L*q1tk1t*t2*P6 + 4*k2(vb)*q1L*
     & q1tk1t*t2*P8 - 4*k2(vb)*q1L*q1tk2t*q2tk1t*P6 + 4*k2(vb)*q2L
     & *q1tq2t*q1tk1t*P5 + 4*k2(vb)*q2L*q1tq2t*q1tk1t*P6 + 4*
     & k2(vb)*q2L*q1tq2t*q1tk1t*P8 - 8*k2(vb)*q2L*q1tq2t*q2tk1t*
     & P5 - 4*k2(vb)*q2L*q1tq2t*P5*m1*m2 + 4*k2(vb)*q2L*q1tq2t*P5*
     & q1k1 + 4*k2(vb)*q2L*q1tq2t*P5*k1k2 + 4*k2(vb)*q2L*q1tq2t*P6*
     & q1k1 + 4*k2(vb)*q2L*q1tq2t*P8*q1k1 + 8*k2(vb)*q2L*q1tk1t*
     & q2tk1t*P5 + 4*k2(vb)*q2L*q1tk1t*q2tk2t*P5 - 4*k2(vb)*q2L*
     & q1tk2t*q2tk1t*P5 - 8*k2(vb)*q2L*q1tk2t*q2tk1t*P6 - 8*k2(vb)
     & *q2L*q1tk2t*q2tk1t*P8
      A52(3,vb) = A52(3,vb) + 4*k2(vb)*q2L*q2tk1t*t1*P5 - 4*k2(vb)*
     & q2L*q2tk1t*t1*P6 - 4*k2(vb)*q2L*q2tk1t*t1*P8 + 4*k2(vb)*
     & q1tq2t*q1tL*P5*q2k1 + 4*k2(vb)*q1tq2t*q1tL*P6*m1*m2 + 8*
     & k2(vb)*q1tq2t*q1tL*P6*q1k1 + 4*k2(vb)*q1tq2t*q1tL*P6*q2k1 - 
     & 4*k2(vb)*q1tq2t*q1tL*P6*k1k2 + 4*k2(vb)*q1tq2t*q1tL*P8*q2k1
     &  + 4*k2(vb)*q1tq2t*q2tL*P5*m1*m2 - 4*k2(vb)*q1tq2t*q2tL*P5*
     & q1k1 - 8*k2(vb)*q1tq2t*q2tL*P5*q2k1 - 4*k2(vb)*q1tq2t*q2tL*
     & P5*k1k2 - 4*k2(vb)*q1tq2t*q2tL*P6*q1k1 - 4*k2(vb)*q1tq2t*
     & q2tL*P8*q1k1 + 4*k2(vb)*q1tk1t*q2tL*t1*P6 + 4*k2(vb)*q1tk1t*
     & q2tL*t2*P5 + 4*k2(vb)*q1tk1t*q2tL*P5*q1q2 + 4*k2(vb)*q1tk1t*
     & q2tL*P5*q2k2 + 4*k2(vb)*q1tk1t*q2tL*P6*q1q2 - 4*k2(vb)*
     & q1tk1t*q2tL*P6*q1k2 + 4*k2(vb)*q1tk1t*q2tL*P8*q1q2 - 4*
     & k2(vb)*q1tk2t*q2tL*P5*q2k1 - 8*k2(vb)*q1tk2t*q2tL*P6*m1*m2
     &  - 4*k2(vb)*q1tk2t*q2tL*P6*q1k1 - 8*k2(vb)*q1tk2t*q2tL*P6*
     & q2k1 + 8*k2(vb)*q1tk2t*q2tL*P6*k1k2 - 8*k2(vb)*q1tk2t*q2tL*
     & P8*q2k1
      A52(3,vb) = A52(3,vb) - 4*k2(vb)*q1tL*q2tk1t*t1*P6 - 4*k2(vb)
     & *q1tL*q2tk1t*t2*P5 - 8*k2(vb)*q1tL*q2tk1t*P5*m1*m2 - 4*
     & k2(vb)*q1tL*q2tk1t*P5*q1q2 + 8*k2(vb)*q1tL*q2tk1t*P5*q1k1 + 
     & 8*k2(vb)*q1tL*q2tk1t*P5*q2k1 - 4*k2(vb)*q1tL*q2tk1t*P5*q2k2
     &  + 8*k2(vb)*q1tL*q2tk1t*P5*k1k2 - 4*k2(vb)*q1tL*q2tk1t*P6*
     & q1q2 + 4*k2(vb)*q1tL*q2tk1t*P6*q1k2 - 4*k2(vb)*q1tL*q2tk1t*
     & P8*q1q2 + 8*k2(vb)*q1tL*q2tk1t*P8*q1k1 + 4*k2(vb)*q1tL*
     & q2tk2t*P5*q2k1 - 4*k2(vb)*q1tL*q2tk2t*P6*q1k1 - 4*k2(vb)*
     & q1tL*t2*P5*m1*m2 + 4*k2(vb)*q1tL*t2*P5*q1k1 + 4*k2(vb)*q1tL*
     & t2*P5*k1k2 - 4*k2(vb)*q1tL*t2*P6*q1k1 + 4*k2(vb)*q1tL*t2*P8*
     & q1k1 + 4*k2(vb)*q2tL*t1*P5*q2k1 - 4*k2(vb)*q2tL*t1*P6*m1*m2 - 
     & 4*k2(vb)*q2tL*t1*P6*q2k1 + 4*k2(vb)*q2tL*t1*P6*k1k2 - 4*k2(vb)
     & *q2tL*t1*P8*q2k1 - 4*kw(vb)*q1L*q1tq2t*P6*m1*m2 + 4*kw(vb)*
     & q1L*q1tq2t*P6*k1k2 + 4*kw(vb)*q1L*q1tk1t*q2tk2t*P6 - 4*
     & kw(vb)*q1L*q1tk2t*q2tk1t*P6 - 4*kw(vb)*q2L*q1tq2t*P5*m1*m2
     &  + 4*kw(vb)*q2L*q1tq2t*P5*k1k2
      A52(3,vb) = A52(3,vb) + 4*kw(vb)*q2L*q1tk1t*q2tk2t*P5 - 4*
     & kw(vb)*q2L*q1tk2t*q2tk1t*P5 + 4*kw(vb)*q1tq2t*q1tL*P6*m1*m2
     &  - 4*kw(vb)*q1tq2t*q1tL*P6*k1k2 + 4*kw(vb)*q1tq2t*q2tL*P5*m1
     & *m2 - 4*kw(vb)*q1tq2t*q2tL*P5*k1k2 + 4*kw(vb)*q1tk1t*q2tL*P5
     & *q2k2 + 4*kw(vb)*q1tk1t*q2tL*P6*q1k2 - 4*kw(vb)*q1tk2t*q2tL*
     & P5*q2k1 - 8*kw(vb)*q1tk2t*q2tL*P6*m1*m2 - 4*kw(vb)*q1tk2t*
     & q2tL*P6*q1k1 + 8*kw(vb)*q1tk2t*q2tL*P6*k1k2 - 8*kw(vb)*q1tL*
     & q2tk1t*P5*m1*m2 - 4*kw(vb)*q1tL*q2tk1t*P5*q2k2 + 8*kw(vb)*
     & q1tL*q2tk1t*P5*k1k2 - 4*kw(vb)*q1tL*q2tk1t*P6*q1k2 + 4*
     & kw(vb)*q1tL*q2tk2t*P5*q2k1 + 4*kw(vb)*q1tL*q2tk2t*P6*q1k1 - 
     & 4*kw(vb)*q1tL*t2*P5*m1*m2 + 4*kw(vb)*q1tL*t2*P5*k1k2 - 4*
     & kw(vb)*q2tL*t1*P6*m1*m2 + 4*kw(vb)*q2tL*t1*P6*k1k2 + 4*q1(vb)*
     & q1L*q1tq2t*P6*m1*m2 - 4*q1(vb)*q1L*q1tq2t*P6*k1k2 - 4*q1(vb)
     & *q1L*q1tk1t*q2tk2t*P6 + 4*q1(vb)*q1L*q1tk2t*q2tk1t*P6 + 4*
     & q1(vb)*q2L*q1tq2t*P6*m1*m2 - 4*q1(vb)*q2L*q1tq2t*P6*k1k2 + 4
     & *q1(vb)*q2L*q1tq2t*P8*m1*m2
      A52(3,vb) = A52(3,vb) - 4*q1(vb)*q2L*q1tq2t*P8*k1k2 - 4*
     & q1(vb)*q2L*q1tk1t*q2tk2t*P6 - 4*q1(vb)*q2L*q1tk1t*q2tk2t*
     & P8 + 4*q1(vb)*q2L*q1tk2t*q2tk1t*P6 + 4*q1(vb)*q2L*q1tk2t*
     & q2tk1t*P8 + 4*q1(vb)*q1tq2t*q1tL*P6*m1*m2 - 4*q1(vb)*q1tq2t*
     & q1tL*P6*k1k2 - 4*q1(vb)*q1tq2t*q2tL*P6*m1*m2 + 4*q1(vb)*
     & q1tq2t*q2tL*P6*k1k2 - 4*q1(vb)*q1tq2t*q2tL*P8*m1*m2 + 4*
     & q1(vb)*q1tq2t*q2tL*P8*k1k2 - 4*q1(vb)*q1tk1t*q2tL*P6*q1k2 - 
     & 4*q1(vb)*q1tk1t*q2tL*P6*q2k2 - 4*q1(vb)*q1tk1t*q2tL*P8*q2k2
     &  + 8*q1(vb)*q1tk2t*q2tL*P6*m1*m2 + 4*q1(vb)*q1tk2t*q2tL*P6*
     & q1k1 + 4*q1(vb)*q1tk2t*q2tL*P6*q2k1 - 8*q1(vb)*q1tk2t*q2tL*
     & P6*k1k2 + 4*q1(vb)*q1tk2t*q2tL*P8*q2k1 + 4*q1(vb)*q1tL*
     & q2tk1t*P6*q1k2 + 4*q1(vb)*q1tL*q2tk1t*P6*q2k2 + 8*q1(vb)*
     & q1tL*q2tk1t*P8*m1*m2 + 4*q1(vb)*q1tL*q2tk1t*P8*q2k2 - 8*
     & q1(vb)*q1tL*q2tk1t*P8*k1k2 - 8*q1(vb)*q1tL*q2tk2t*P6*m1*m2
     &  - 4*q1(vb)*q1tL*q2tk2t*P6*q1k1 - 4*q1(vb)*q1tL*q2tk2t*P6*
     & q2k1
      A52(3,vb) = A52(3,vb) + 8*q1(vb)*q1tL*q2tk2t*P6*k1k2 - 4*
     & q1(vb)*q1tL*q2tk2t*P8*q2k1 - 4*q1(vb)*q1tL*t2*P6*m1*m2 + 4*
     & q1(vb)*q1tL*t2*P6*k1k2 + 4*q1(vb)*q1tL*t2*P8*m1*m2 - 4*q1(vb)*
     & q1tL*t2*P8*k1k2 + 4*q1(vb)*q2tL*t1*P6*m1*m2 - 4*q1(vb)*q2tL*
     & t1*P6*k1k2 + 4*q2(vb)*q1L*q1tq2t*P5*m1*m2 - 4*q2(vb)*q1L*
     & q1tq2t*P5*k1k2 + 4*q2(vb)*q1L*q1tq2t*P8*m1*m2 - 4*q2(vb)*q1L
     & *q1tq2t*P8*k1k2 - 4*q2(vb)*q1L*q1tk1t*q2tk2t*P5 - 4*q2(vb)*
     & q1L*q1tk1t*q2tk2t*P8 + 4*q2(vb)*q1L*q1tk2t*q2tk1t*P5 + 4*
     & q2(vb)*q1L*q1tk2t*q2tk1t*P8 + 4*q2(vb)*q2L*q1tq2t*P5*m1*m2
     &  - 4*q2(vb)*q2L*q1tq2t*P5*k1k2 - 4*q2(vb)*q2L*q1tk1t*q2tk2t
     & *P5 + 4*q2(vb)*q2L*q1tk2t*q2tk1t*P5 - 4*q2(vb)*q1tq2t*q1tL*
     & P5*m1*m2 + 4*q2(vb)*q1tq2t*q1tL*P5*k1k2 - 4*q2(vb)*q1tq2t*
     & q1tL*P8*m1*m2 + 4*q2(vb)*q1tq2t*q1tL*P8*k1k2 + 4*q2(vb)*
     & q1tq2t*q2tL*P5*m1*m2 - 4*q2(vb)*q1tq2t*q2tL*P5*k1k2 - 8*
     & q2(vb)*q1tk1t*q2tL*P5*m1*m2 - 4*q2(vb)*q1tk1t*q2tL*P5*q1k2
     &  - 4*q2(vb)*q1tk1t*q2tL*P5*q2k2
      A52(3,vb) = A52(3,vb) + 8*q2(vb)*q1tk1t*q2tL*P5*k1k2 - 4*
     & q2(vb)*q1tk1t*q2tL*P8*q1k2 + 4*q2(vb)*q1tk2t*q2tL*P5*q1k1 + 
     & 4*q2(vb)*q1tk2t*q2tL*P5*q2k1 + 8*q2(vb)*q1tk2t*q2tL*P8*m1*m2
     &  + 4*q2(vb)*q1tk2t*q2tL*P8*q1k1 - 8*q2(vb)*q1tk2t*q2tL*P8*
     & k1k2 + 8*q2(vb)*q1tL*q2tk1t*P5*m1*m2 + 4*q2(vb)*q1tL*q2tk1t*
     & P5*q1k2 + 4*q2(vb)*q1tL*q2tk1t*P5*q2k2 - 8*q2(vb)*q1tL*
     & q2tk1t*P5*k1k2 + 4*q2(vb)*q1tL*q2tk1t*P8*q1k2 - 4*q2(vb)*
     & q1tL*q2tk2t*P5*q1k1 - 4*q2(vb)*q1tL*q2tk2t*P5*q2k1 - 4*
     & q2(vb)*q1tL*q2tk2t*P8*q1k1 + 4*q2(vb)*q1tL*t2*P5*m1*m2 - 4*
     & q2(vb)*q1tL*t2*P5*k1k2 - 4*q2(vb)*q2tL*t1*P5*m1*m2 + 4*q2(vb)*
     & q2tL*t1*P5*k1k2 + 4*q2(vb)*q2tL*t1*P8*m1*m2 - 4*q2(vb)*q2tL*
     & t1*P8*k1k2 - 8*q1t(vb)*q1L*q1tq2t*P6*m1*m2 + 8*q1t(vb)*q1L*
     & q1tq2t*P6*k1k2 + 8*q1t(vb)*q1L*q2tk1t*P5*m1*m2 + 4*q1t(vb)*
     & q1L*q2tk1t*P5*q2k2 - 8*q1t(vb)*q1L*q2tk1t*P5*k1k2 - 4*
     & q1t(vb)*q1L*q2tk1t*P6*q2k2 - 8*q1t(vb)*q1L*q2tk1t*P8*m1*m2
     &  - 4*q1t(vb)*q1L*q2tk1t*P8*q2k2
      A52(3,vb) = A52(3,vb) + 8*q1t(vb)*q1L*q2tk1t*P8*k1k2 - 4*
     & q1t(vb)*q1L*q2tk2t*P5*q2k1 + 8*q1t(vb)*q1L*q2tk2t*P6*m1*m2
     &  + 4*q1t(vb)*q1L*q2tk2t*P6*q2k1 - 8*q1t(vb)*q1L*q2tk2t*P6*
     & k1k2 + 4*q1t(vb)*q1L*q2tk2t*P8*q2k1 + 4*q1t(vb)*q1L*t2*P5*m1*
     & m2 - 4*q1t(vb)*q1L*t2*P5*k1k2 + 4*q1t(vb)*q1L*t2*P6*m1*m2 - 4*
     & q1t(vb)*q1L*t2*P6*k1k2 - 4*q1t(vb)*q1L*t2*P8*m1*m2 + 4*q1t(vb)
     & *q1L*t2*P8*k1k2 + 4*q1t(vb)*q2L*q1tq2t*P5*m1*m2 - 4*q1t(vb)*
     & q2L*q1tq2t*P5*k1k2 - 4*q1t(vb)*q2L*q1tq2t*P6*m1*m2 + 4*
     & q1t(vb)*q2L*q1tq2t*P6*k1k2 - 4*q1t(vb)*q2L*q1tq2t*P8*m1*m2
     &  + 4*q1t(vb)*q2L*q1tq2t*P8*k1k2 - 4*q1t(vb)*q2L*q2tk1t*P5*
     & q1k2 + 4*q1t(vb)*q2L*q2tk1t*P6*q1k2 + 4*q1t(vb)*q2L*q2tk1t*
     & P8*q1k2 + 4*q1t(vb)*q2L*q2tk2t*P5*q1k1 - 4*q1t(vb)*q2L*
     & q2tk2t*P6*q1k1 - 4*q1t(vb)*q2L*q2tk2t*P8*q1k1 - 4*q1t(vb)*
     & q2tL*t1*P6*m1*m2 + 4*q1t(vb)*q2tL*t1*P6*k1k2 + 4*q1t(vb)*q2tL
     & *t2*P5*m1*m2 - 4*q1t(vb)*q2tL*t2*P5*k1k2 + 4*q1t(vb)*q2tL*P5*
     & m1*m2*q1q2
      A52(3,vb) = A52(3,vb) - 8*q1t(vb)*q2tL*P5*m1*m2*q2k1 - 4*
     & q1t(vb)*q2tL*P5*q1q2*k1k2 + 4*q1t(vb)*q2tL*P5*q1k1*q2k2 - 4*
     & q1t(vb)*q2tL*P5*q1k2*q2k1 + 8*q1t(vb)*q2tL*P5*q2k1*k1k2 - 4*
     & q1t(vb)*q2tL*P6*m1*m2*q1q2 + 8*q1t(vb)*q2tL*P6*m1*m2*q1k2 + 4*
     & q1t(vb)*q2tL*P6*q1q2*k1k2 - 4*q1t(vb)*q2tL*P6*q1k1*q2k2 + 4*
     & q1t(vb)*q2tL*P6*q1k2*q2k1 - 8*q1t(vb)*q2tL*P6*q1k2*k1k2 - 4*
     & q1t(vb)*q2tL*P8*m1*m2*q1q2 + 4*q1t(vb)*q2tL*P8*q1q2*k1k2 - 4*
     & q1t(vb)*q2tL*P8*q1k1*q2k2 + 4*q1t(vb)*q2tL*P8*q1k2*q2k1 - 4*
     & q2t(vb)*q1L*q1tq2t*P5*m1*m2 + 4*q2t(vb)*q1L*q1tq2t*P5*k1k2
     &  + 4*q2t(vb)*q1L*q1tq2t*P6*m1*m2 - 4*q2t(vb)*q1L*q1tq2t*P6*
     & k1k2 - 4*q2t(vb)*q1L*q1tq2t*P8*m1*m2 + 4*q2t(vb)*q1L*q1tq2t*
     & P8*k1k2 - 4*q2t(vb)*q1L*q1tk1t*P5*q2k2 + 4*q2t(vb)*q1L*
     & q1tk1t*P6*q2k2 - 4*q2t(vb)*q1L*q1tk1t*P8*q2k2 + 4*q2t(vb)*
     & q1L*q1tk2t*P5*q2k1 - 4*q2t(vb)*q1L*q1tk2t*P6*q2k1 + 4*
     & q2t(vb)*q1L*q1tk2t*P8*q2k1 - 8*q2t(vb)*q2L*q1tq2t*P5*m1*m2
     &  + 8*q2t(vb)*q2L*q1tq2t*P5*k1k2
      A52(3,vb) = A52(3,vb) + 8*q2t(vb)*q2L*q1tk1t*P5*m1*m2 + 4*
     & q2t(vb)*q2L*q1tk1t*P5*q1k2 - 8*q2t(vb)*q2L*q1tk1t*P5*k1k2 - 
     & 4*q2t(vb)*q2L*q1tk1t*P6*q1k2 + 4*q2t(vb)*q2L*q1tk1t*P8*q1k2
     &  - 4*q2t(vb)*q2L*q1tk2t*P5*q1k1 + 8*q2t(vb)*q2L*q1tk2t*P6*m1
     & *m2 + 4*q2t(vb)*q2L*q1tk2t*P6*q1k1 - 8*q2t(vb)*q2L*q1tk2t*P6
     & *k1k2 - 8*q2t(vb)*q2L*q1tk2t*P8*m1*m2 - 4*q2t(vb)*q2L*q1tk2t
     & *P8*q1k1 + 8*q2t(vb)*q2L*q1tk2t*P8*k1k2 + 4*q2t(vb)*q2L*t1*P5
     & *m1*m2 - 4*q2t(vb)*q2L*t1*P5*k1k2 + 4*q2t(vb)*q2L*t1*P6*m1*m2
     &  - 4*q2t(vb)*q2L*t1*P6*k1k2 - 4*q2t(vb)*q2L*t1*P8*m1*m2 + 4*
     & q2t(vb)*q2L*t1*P8*k1k2 + 4*q2t(vb)*q1tL*t1*P6*m1*m2 - 4*
     & q2t(vb)*q1tL*t1*P6*k1k2 - 4*q2t(vb)*q1tL*t2*P5*m1*m2 + 4*
     & q2t(vb)*q1tL*t2*P5*k1k2 - 4*q2t(vb)*q1tL*P5*m1*m2*q1q2 + 8*
     & q2t(vb)*q1tL*P5*m1*m2*q2k1 + 4*q2t(vb)*q1tL*P5*q1q2*k1k2 - 4*
     & q2t(vb)*q1tL*P5*q1k1*q2k2 + 4*q2t(vb)*q1tL*P5*q1k2*q2k1 - 8*
     & q2t(vb)*q1tL*P5*q2k1*k1k2 + 4*q2t(vb)*q1tL*P6*m1*m2*q1q2 - 8*
     & q2t(vb)*q1tL*P6*m1*m2*q1k2
      A52(3,vb) = A52(3,vb) - 4*q2t(vb)*q1tL*P6*q1q2*k1k2 + 4*
     & q2t(vb)*q1tL*P6*q1k1*q2k2 - 4*q2t(vb)*q1tL*P6*q1k2*q2k1 + 8*
     & q2t(vb)*q1tL*P6*q1k2*k1k2 - 4*q2t(vb)*q1tL*P8*m1*m2*q1q2 + 4*
     & q2t(vb)*q1tL*P8*q1q2*k1k2 - 4*q2t(vb)*q1tL*P8*q1k1*q2k2 + 4*
     & q2t(vb)*q1tL*P8*q1k2*q2k1 - 4*L(vb)*t1*t2*P5*m1*m2 + 4*L(vb)*t1
     & *t2*P5*k1k2 - 4*L(vb)*t1*t2*P6*m1*m2 + 4*L(vb)*t1*t2*P6*k1k2 + 4
     & *L(vb)*t1*t2*P8*m1*m2 - 4*L(vb)*t1*t2*P8*k1k2 - 4*L(vb)*q1tq2t*
     & q1tk1t*P5*q2k2 - 8*L(vb)*q1tq2t*q1tk1t*P6*q1k2 - 4*L(vb)*
     & q1tq2t*q1tk1t*P6*q2k2 - 4*L(vb)*q1tq2t*q1tk1t*P8*q2k2 + 4*
     & L(vb)*q1tq2t*q1tk2t*P5*q2k1 + 16*L(vb)*q1tq2t*q1tk2t*P6*m1*
     & m2 + 8*L(vb)*q1tq2t*q1tk2t*P6*q1k1 + 4*L(vb)*q1tq2t*q1tk2t*
     & P6*q2k1 - 16*L(vb)*q1tq2t*q1tk2t*P6*k1k2 + 4*L(vb)*q1tq2t*
     & q1tk2t*P8*q2k1 + 16*L(vb)*q1tq2t*q2tk1t*P5*m1*m2 + 4*L(vb)*
     & q1tq2t*q2tk1t*P5*q1k2 + 8*L(vb)*q1tq2t*q2tk1t*P5*q2k2 - 16*
     & L(vb)*q1tq2t*q2tk1t*P5*k1k2 + 4*L(vb)*q1tq2t*q2tk1t*P6*q1k2
     &  + 4*L(vb)*q1tq2t*q2tk1t*P8*q1k2
      A52(3,vb) = A52(3,vb) - 4*L(vb)*q1tq2t*q2tk2t*P5*q1k1 - 8*
     & L(vb)*q1tq2t*q2tk2t*P5*q2k1 - 4*L(vb)*q1tq2t*q2tk2t*P6*q1k1
     &  - 4*L(vb)*q1tq2t*q2tk2t*P8*q1k1 + 4*L(vb)*q1tq2t*t1*P6*m1*m2
     &  - 4*L(vb)*q1tq2t*t1*P6*k1k2 + 4*L(vb)*q1tq2t*t2*P5*m1*m2 - 4*
     & L(vb)*q1tq2t*t2*P5*k1k2 - 4*L(vb)*q1tq2t*P5*m1*m2*q1q2 + 8*
     & L(vb)*q1tq2t*P5*m1*m2*q2k1 + 4*L(vb)*q1tq2t*P5*q1q2*k1k2 - 4*
     & L(vb)*q1tq2t*P5*q1k1*q2k2 + 4*L(vb)*q1tq2t*P5*q1k2*q2k1 - 8*
     & L(vb)*q1tq2t*P5*q2k1*k1k2 - 4*L(vb)*q1tq2t*P6*m1*m2*q1q2 + 8*
     & L(vb)*q1tq2t*P6*m1*m2*q1k2 + 4*L(vb)*q1tq2t*P6*q1q2*k1k2 - 4*
     & L(vb)*q1tq2t*P6*q1k1*q2k2 + 4*L(vb)*q1tq2t*P6*q1k2*q2k1 - 8*
     & L(vb)*q1tq2t*P6*q1k2*k1k2 - 4*L(vb)*q1tq2t*P8*m1*m2*q1q2 + 4*
     & L(vb)*q1tq2t*P8*q1q2*k1k2 - 4*L(vb)*q1tq2t*P8*q1k1*q2k2 + 4*
     & L(vb)*q1tq2t*P8*q1k2*q2k1 - 4*L(vb)*q1tq2t**2*P5*m1*m2 + 4*
     & L(vb)*q1tq2t**2*P5*k1k2 - 4*L(vb)*q1tq2t**2*P6*m1*m2 + 4*L(vb)
     & *q1tq2t**2*P6*k1k2 - 4*L(vb)*q1tq2t**2*P8*m1*m2 + 4*L(vb)*
     & q1tq2t**2*P8*k1k2
      A52(3,vb) = A52(3,vb) - 16*L(vb)*q1tk1t*q2tk1t*P5*m1*m2 - 8*
     & L(vb)*q1tk1t*q2tk1t*P5*q1k2 - 8*L(vb)*q1tk1t*q2tk1t*P5*q2k2
     &  + 16*L(vb)*q1tk1t*q2tk1t*P5*k1k2 - 8*L(vb)*q1tk1t*q2tk1t*P8
     & *q1k2 + 4*L(vb)*q1tk1t*q2tk2t*t1*P6 + 4*L(vb)*q1tk1t*q2tk2t*
     & t2*P5 + 4*L(vb)*q1tk1t*q2tk2t*P5*q1q2 + 4*L(vb)*q1tk1t*
     & q2tk2t*P6*q1q2 + 4*L(vb)*q1tk1t*q2tk2t*P8*q1q2 - 8*L(vb)*
     & q1tk1t*t2*P5*m1*m2 - 4*L(vb)*q1tk1t*t2*P5*q1k2 + 8*L(vb)*
     & q1tk1t*t2*P5*k1k2 + 4*L(vb)*q1tk1t*t2*P6*q1k2 - 4*L(vb)*
     & q1tk1t*t2*P8*q1k2 - 4*L(vb)*q1tk2t*q2tk1t*t1*P6 - 4*L(vb)*
     & q1tk2t*q2tk1t*t2*P5 - 4*L(vb)*q1tk2t*q2tk1t*P5*q1q2 + 8*
     & L(vb)*q1tk2t*q2tk1t*P5*q1k1 + 8*L(vb)*q1tk2t*q2tk1t*P5*q2k1
     &  - 4*L(vb)*q1tk2t*q2tk1t*P6*q1q2 + 8*L(vb)*q1tk2t*q2tk1t*P6*
     & q1k2 + 8*L(vb)*q1tk2t*q2tk1t*P6*q2k2 + 16*L(vb)*q1tk2t*
     & q2tk1t*P8*m1*m2 - 4*L(vb)*q1tk2t*q2tk1t*P8*q1q2 + 8*L(vb)*
     & q1tk2t*q2tk1t*P8*q1k1 + 8*L(vb)*q1tk2t*q2tk1t*P8*q2k2 - 16*
     & L(vb)*q1tk2t*q2tk1t*P8*k1k2
      A52(3,vb) = A52(3,vb) - 16*L(vb)*q1tk2t*q2tk2t*P6*m1*m2 - 8*
     & L(vb)*q1tk2t*q2tk2t*P6*q1k1 - 8*L(vb)*q1tk2t*q2tk2t*P6*q2k1
     &  + 16*L(vb)*q1tk2t*q2tk2t*P6*k1k2 - 8*L(vb)*q1tk2t*q2tk2t*P8
     & *q2k1 + 4*L(vb)*q1tk2t*t2*P5*q1k1 - 8*L(vb)*q1tk2t*t2*P6*m1*m2
     &  - 4*L(vb)*q1tk2t*t2*P6*q1k1 + 8*L(vb)*q1tk2t*t2*P6*k1k2 + 8*
     & L(vb)*q1tk2t*t2*P8*m1*m2 + 4*L(vb)*q1tk2t*t2*P8*q1k1 - 8*L(vb)
     & *q1tk2t*t2*P8*k1k2 - 8*L(vb)*q2tk1t*t1*P5*m1*m2 - 4*L(vb)*
     & q2tk1t*t1*P5*q2k2 + 8*L(vb)*q2tk1t*t1*P5*k1k2 + 4*L(vb)*
     & q2tk1t*t1*P6*q2k2 + 8*L(vb)*q2tk1t*t1*P8*m1*m2 + 4*L(vb)*
     & q2tk1t*t1*P8*q2k2 - 8*L(vb)*q2tk1t*t1*P8*k1k2 + 4*L(vb)*
     & q2tk2t*t1*P5*q2k1 - 8*L(vb)*q2tk2t*t1*P6*m1*m2 - 4*L(vb)*
     & q2tk2t*t1*P6*q2k1 + 8*L(vb)*q2tk2t*t1*P6*k1k2 - 4*L(vb)*
     & q2tk2t*t1*P8*q2k1

      A52(4,vb) =
     &  - 4*k1(vb)*q1L*q2U*q1tq2t*P5*m2 - 4*k1(vb)*q1L*q2U*q1tq2t
     & *P6*m2 - 4*k1(vb)*q1L*q2U*q1tq2t*P8*m2 - 8*k1(vb)*q1L*
     & q1tq2t*q1tU*P6*m2 + 4*k1(vb)*q1L*q1tq2t*q2tU*P5*m2 + 4*
     & k1(vb)*q1L*q1tq2t*q2tU*P6*m2 + 4*k1(vb)*q1L*q1tq2t*q2tU*P8
     & *m2 + 4*k1(vb)*q1L*q1tk1t*q2tU*P6*m2 + 4*k1(vb)*q1L*q1tk2t*
     & q2tU*P6*m1 - 8*k1(vb)*q1L*q1tU*q2tk1t*P5*m2 - 4*k1(vb)*q1L*
     & q1tU*q2tk1t*P6*m2 - 8*k1(vb)*q1L*q1tU*q2tk1t*P8*m2 - 4*
     & k1(vb)*q1L*q1tU*q2tk2t*P6*m1 + 8*k1(vb)*q1L*q1tU*q2tk2t*P6
     & *m2 - 4*k1(vb)*q1L*q1tU*t2*P5*m2 + 4*k1(vb)*q1L*q1tU*t2*P6*
     & m2 - 4*k1(vb)*q1L*q1tU*t2*P8*m2 - 4*k1(vb)*q2U*q2L*q1tq2t*
     & P5*m2 - 4*k1(vb)*q2U*q2L*q1tq2t*P6*m2 - 4*k1(vb)*q2U*q2L*
     & q1tq2t*P8*m2 + 4*k1(vb)*q2U*q1tq2t*q1tL*P5*m2 - 4*k1(vb)*
     & q2U*q1tq2t*q1tL*P6*m2 + 4*k1(vb)*q2U*q1tq2t*q1tL*P8*m2 - 4
     & *k1(vb)*q2U*q1tq2t*q2tL*P5*m2 + 4*k1(vb)*q2U*q1tq2t*q2tL*
     & P6*m2
      A52(4,vb) = A52(4,vb) + 4*k1(vb)*q2U*q1tq2t*q2tL*P8*m2 + 4*
     & k1(vb)*q2U*q1tk1t*q2tL*P5*m2 + 4*k1(vb)*q2U*q1tk1t*q2tL*P6
     & *m2 - 4*k1(vb)*q2U*q1tk2t*q2tL*P5*m1 + 4*k1(vb)*q2U*q1tk2t*
     & q2tL*P6*m1 - 8*k1(vb)*q2U*q1tk2t*q2tL*P6*m2 - 8*k1(vb)*q2U*
     & q1tk2t*q2tL*P8*m2 - 4*k1(vb)*q2U*q1tL*q2tk1t*P5*m2 - 4*
     & k1(vb)*q2U*q1tL*q2tk1t*P6*m2 - 8*k1(vb)*q2U*q1tL*q2tk1t*P8
     & *m2 + 4*k1(vb)*q2U*q1tL*q2tk2t*P5*m1 - 4*k1(vb)*q2U*q1tL*
     & q2tk2t*P6*m1 + 8*k1(vb)*q2U*q1tL*q2tk2t*P6*m2 - 4*k1(vb)*
     & q2U*q1tL*t2*P5*m2 + 4*k1(vb)*q2U*q1tL*t2*P6*m2 - 4*k1(vb)*
     & q2U*q1tL*t2*P8*m2 + 4*k1(vb)*q2U*q2tL*t1*P5*m2 - 4*k1(vb)*
     & q2U*q2tL*t1*P6*m2 - 4*k1(vb)*q2U*q2tL*t1*P8*m2 - 4*k1(vb)*
     & q2L*q1tq2t*q1tU*P5*m2 - 4*k1(vb)*q2L*q1tq2t*q1tU*P6*m2 - 4
     & *k1(vb)*q2L*q1tq2t*q1tU*P8*m2 + 8*k1(vb)*q2L*q1tq2t*q2tU*
     & P5*m2 - 4*k1(vb)*q2L*q1tk1t*q2tU*P5*m2 + 4*k1(vb)*q2L*
     & q1tk2t*q2tU*P5*m1 + 8*k1(vb)*q2L*q1tk2t*q2tU*P6*m2 + 8*
     & k1(vb)*q2L*q1tk2t*q2tU*P8*m2
      A52(4,vb) = A52(4,vb) - 4*k1(vb)*q2L*q1tU*q2tk1t*P5*m2 - 4*
     & k1(vb)*q2L*q1tU*q2tk2t*P5*m1 - 4*k1(vb)*q2L*q2tU*t1*P5*m2
     &  + 4*k1(vb)*q2L*q2tU*t1*P6*m2 + 4*k1(vb)*q2L*q2tU*t1*P8*m2
     &  - 4*k1(vb)*q1tU*q2tL*t1*P6*m2 - 4*k1(vb)*q1tU*q2tL*t2*P5*m2
     &  - 4*k1(vb)*q1tU*q2tL*P5*m1*q2k2 - 4*k1(vb)*q1tU*q2tL*P5*m2*
     & q1q2 + 4*k1(vb)*q1tU*q2tL*P5*m2*q2k1 - 4*k1(vb)*q1tU*q2tL*P6
     & *m1*q1k2 - 4*k1(vb)*q1tU*q2tL*P6*m2*q1q2 - 4*k1(vb)*q1tU*
     & q2tL*P6*m2*q1k1 + 8*k1(vb)*q1tU*q2tL*P6*m2*q1k2 - 4*k1(vb)*
     & q1tU*q2tL*P8*m2*q1q2 + 4*k1(vb)*q1tL*q2tU*t1*P6*m2 + 4*
     & k1(vb)*q1tL*q2tU*t2*P5*m2 + 4*k1(vb)*q1tL*q2tU*P5*m1*q2k2 + 
     & 4*k1(vb)*q1tL*q2tU*P5*m2*q1q2 - 4*k1(vb)*q1tL*q2tU*P5*m2*
     & q2k1 + 4*k1(vb)*q1tL*q2tU*P6*m1*q1k2 + 4*k1(vb)*q1tL*q2tU*P6
     & *m2*q1q2 + 4*k1(vb)*q1tL*q2tU*P6*m2*q1k1 - 8*k1(vb)*q1tL*
     & q2tU*P6*m2*q1k2 + 4*k1(vb)*q1tL*q2tU*P8*m2*q1q2 - 4*k2(vb)*
     & q1L*q2U*q1tq2t*P5*m1 - 4*k2(vb)*q1L*q2U*q1tq2t*P6*m1 - 4*
     & k2(vb)*q1L*q2U*q1tq2t*P8*m1
      A52(4,vb) = A52(4,vb) - 8*k2(vb)*q1L*q1tq2t*q1tU*P6*m1 + 4*
     & k2(vb)*q1L*q1tq2t*q2tU*P5*m1 + 4*k2(vb)*q1L*q1tq2t*q2tU*P6
     & *m1 + 4*k2(vb)*q1L*q1tq2t*q2tU*P8*m1 + 4*k2(vb)*q1L*q1tk1t*
     & q2tU*P6*m2 + 4*k2(vb)*q1L*q1tk2t*q2tU*P6*m1 - 8*k2(vb)*q1L*
     & q1tU*q2tk1t*P5*m1 - 4*k2(vb)*q1L*q1tU*q2tk1t*P6*m2 - 8*
     & k2(vb)*q1L*q1tU*q2tk1t*P8*m1 + 4*k2(vb)*q1L*q1tU*q2tk2t*P6
     & *m1 - 4*k2(vb)*q1L*q1tU*t2*P5*m1 + 4*k2(vb)*q1L*q1tU*t2*P6*
     & m1 - 4*k2(vb)*q1L*q1tU*t2*P8*m1 - 4*k2(vb)*q2U*q2L*q1tq2t*
     & P5*m1 - 4*k2(vb)*q2U*q2L*q1tq2t*P6*m1 - 4*k2(vb)*q2U*q2L*
     & q1tq2t*P8*m1 + 4*k2(vb)*q2U*q1tq2t*q1tL*P5*m1 - 4*k2(vb)*
     & q2U*q1tq2t*q1tL*P6*m1 + 4*k2(vb)*q2U*q1tq2t*q1tL*P8*m1 - 4
     & *k2(vb)*q2U*q1tq2t*q2tL*P5*m1 + 4*k2(vb)*q2U*q1tq2t*q2tL*
     & P6*m1 + 4*k2(vb)*q2U*q1tq2t*q2tL*P8*m1 + 8*k2(vb)*q2U*
     & q1tk1t*q2tL*P5*m1 - 4*k2(vb)*q2U*q1tk1t*q2tL*P5*m2 + 4*
     & k2(vb)*q2U*q1tk1t*q2tL*P6*m2 - 4*k2(vb)*q2U*q1tk2t*q2tL*P5
     & *m1
      A52(4,vb) = A52(4,vb) - 4*k2(vb)*q2U*q1tk2t*q2tL*P6*m1 - 8*
     & k2(vb)*q2U*q1tk2t*q2tL*P8*m1 - 8*k2(vb)*q2U*q1tL*q2tk1t*P5
     & *m1 + 4*k2(vb)*q2U*q1tL*q2tk1t*P5*m2 - 4*k2(vb)*q2U*q1tL*
     & q2tk1t*P6*m2 - 8*k2(vb)*q2U*q1tL*q2tk1t*P8*m1 + 4*k2(vb)*
     & q2U*q1tL*q2tk2t*P5*m1 + 4*k2(vb)*q2U*q1tL*q2tk2t*P6*m1 - 4
     & *k2(vb)*q2U*q1tL*t2*P5*m1 + 4*k2(vb)*q2U*q1tL*t2*P6*m1 - 4*
     & k2(vb)*q2U*q1tL*t2*P8*m1 + 4*k2(vb)*q2U*q2tL*t1*P5*m1 - 4*
     & k2(vb)*q2U*q2tL*t1*P6*m1 - 4*k2(vb)*q2U*q2tL*t1*P8*m1 - 4*
     & k2(vb)*q2L*q1tq2t*q1tU*P5*m1 - 4*k2(vb)*q2L*q1tq2t*q1tU*P6
     & *m1 - 4*k2(vb)*q2L*q1tq2t*q1tU*P8*m1 + 8*k2(vb)*q2L*q1tq2t*
     & q2tU*P5*m1 - 8*k2(vb)*q2L*q1tk1t*q2tU*P5*m1 + 4*k2(vb)*q2L*
     & q1tk1t*q2tU*P5*m2 + 4*k2(vb)*q2L*q1tk2t*q2tU*P5*m1 + 8*
     & k2(vb)*q2L*q1tk2t*q2tU*P6*m1 + 8*k2(vb)*q2L*q1tk2t*q2tU*P8
     & *m1 - 4*k2(vb)*q2L*q1tU*q2tk1t*P5*m2 - 4*k2(vb)*q2L*q1tU*
     & q2tk2t*P5*m1 - 4*k2(vb)*q2L*q2tU*t1*P5*m1 + 4*k2(vb)*q2L*
     & q2tU*t1*P6*m1
      A52(4,vb) = A52(4,vb) + 4*k2(vb)*q2L*q2tU*t1*P8*m1 - 4*k2(vb)
     & *q1tU*q2tL*t1*P6*m1 - 4*k2(vb)*q1tU*q2tL*t2*P5*m1 - 4*k2(vb)
     & *q1tU*q2tL*P5*m1*q1q2 + 8*k2(vb)*q1tU*q2tL*P5*m1*q2k1 - 4*
     & k2(vb)*q1tU*q2tL*P5*m1*q2k2 - 4*k2(vb)*q1tU*q2tL*P5*m2*q2k1
     &  - 4*k2(vb)*q1tU*q2tL*P6*m1*q1q2 + 4*k2(vb)*q1tU*q2tL*P6*m1*
     & q1k2 - 4*k2(vb)*q1tU*q2tL*P6*m2*q1k1 - 4*k2(vb)*q1tU*q2tL*P8
     & *m1*q1q2 + 4*k2(vb)*q1tL*q2tU*t1*P6*m1 + 4*k2(vb)*q1tL*q2tU*
     & t2*P5*m1 + 4*k2(vb)*q1tL*q2tU*P5*m1*q1q2 - 8*k2(vb)*q1tL*
     & q2tU*P5*m1*q2k1 + 4*k2(vb)*q1tL*q2tU*P5*m1*q2k2 + 4*k2(vb)*
     & q1tL*q2tU*P5*m2*q2k1 + 4*k2(vb)*q1tL*q2tU*P6*m1*q1q2 - 4*
     & k2(vb)*q1tL*q2tU*P6*m1*q1k2 + 4*k2(vb)*q1tL*q2tU*P6*m2*q1k1
     &  + 4*k2(vb)*q1tL*q2tU*P8*m1*q1q2 + 4*kw(vb)*q1L*q1tk1t*q2tU
     & *P6*m2 + 4*kw(vb)*q1L*q1tk2t*q2tU*P6*m1 - 4*kw(vb)*q1L*q1tU
     & *q2tk1t*P6*m2 - 4*kw(vb)*q1L*q1tU*q2tk2t*P6*m1 - 4*kw(vb)*
     & q2U*q1tk1t*q2tL*P5*m2 + 4*kw(vb)*q2U*q1tk1t*q2tL*P6*m2 - 4
     & *kw(vb)*q2U*q1tk2t*q2tL*P5*m1
      A52(4,vb) = A52(4,vb) + 4*kw(vb)*q2U*q1tk2t*q2tL*P6*m1 + 4*
     & kw(vb)*q2U*q1tL*q2tk1t*P5*m2 - 4*kw(vb)*q2U*q1tL*q2tk1t*P6
     & *m2 + 4*kw(vb)*q2U*q1tL*q2tk2t*P5*m1 - 4*kw(vb)*q2U*q1tL*
     & q2tk2t*P6*m1 + 4*kw(vb)*q2L*q1tk1t*q2tU*P5*m2 + 4*kw(vb)*
     & q2L*q1tk2t*q2tU*P5*m1 - 4*kw(vb)*q2L*q1tU*q2tk1t*P5*m2 - 4
     & *kw(vb)*q2L*q1tU*q2tk2t*P5*m1 - 4*kw(vb)*q1tU*q2tL*P5*m1*
     & q2k2 - 4*kw(vb)*q1tU*q2tL*P5*m2*q2k1 - 4*kw(vb)*q1tU*q2tL*P6
     & *m1*q1k2 - 4*kw(vb)*q1tU*q2tL*P6*m2*q1k1 + 4*kw(vb)*q1tL*
     & q2tU*P5*m1*q2k2 + 4*kw(vb)*q1tL*q2tU*P5*m2*q2k1 + 4*kw(vb)*
     & q1tL*q2tU*P6*m1*q1k2 + 4*kw(vb)*q1tL*q2tU*P6*m2*q1k1 - 4*
     & q1(vb)*q1L*q1tk1t*q2tU*P6*m2 - 4*q1(vb)*q1L*q1tk2t*q2tU*P6
     & *m1 + 4*q1(vb)*q1L*q1tU*q2tk1t*P6*m2 + 4*q1(vb)*q1L*q1tU*
     & q2tk2t*P6*m1 + 4*q1(vb)*q2U*q1tk1t*q2tL*P8*m2 + 4*q1(vb)*
     & q2U*q1tk2t*q2tL*P8*m1 - 4*q1(vb)*q2U*q1tL*q2tk1t*P8*m2 - 4
     & *q1(vb)*q2U*q1tL*q2tk2t*P8*m1 - 4*q1(vb)*q2L*q1tk1t*q2tU*
     & P6*m2
      A52(4,vb) = A52(4,vb) - 4*q1(vb)*q2L*q1tk1t*q2tU*P8*m2 - 4*
     & q1(vb)*q2L*q1tk2t*q2tU*P6*m1 - 4*q1(vb)*q2L*q1tk2t*q2tU*P8
     & *m1 + 4*q1(vb)*q2L*q1tU*q2tk1t*P6*m2 + 4*q1(vb)*q2L*q1tU*
     & q2tk1t*P8*m2 + 4*q1(vb)*q2L*q1tU*q2tk2t*P6*m1 + 4*q1(vb)*
     & q2L*q1tU*q2tk2t*P8*m1 + 4*q1(vb)*q1tU*q2tL*P6*m1*q1k2 + 4*
     & q1(vb)*q1tU*q2tL*P6*m1*q2k2 + 4*q1(vb)*q1tU*q2tL*P6*m2*q1k1
     &  + 4*q1(vb)*q1tU*q2tL*P6*m2*q2k1 + 4*q1(vb)*q1tU*q2tL*P8*m1*
     & q2k2 + 4*q1(vb)*q1tU*q2tL*P8*m2*q2k1 - 4*q1(vb)*q1tL*q2tU*P6
     & *m1*q1k2 - 4*q1(vb)*q1tL*q2tU*P6*m1*q2k2 - 4*q1(vb)*q1tL*
     & q2tU*P6*m2*q1k1 - 4*q1(vb)*q1tL*q2tU*P6*m2*q2k1 - 4*q1(vb)*
     & q1tL*q2tU*P8*m1*q2k2 - 4*q1(vb)*q1tL*q2tU*P8*m2*q2k1 - 4*
     & q2(vb)*q1L*q1tk1t*q2tU*P5*m2 - 4*q2(vb)*q1L*q1tk1t*q2tU*P8
     & *m2 - 4*q2(vb)*q1L*q1tk2t*q2tU*P5*m1 - 4*q2(vb)*q1L*q1tk2t*
     & q2tU*P8*m1 + 4*q2(vb)*q1L*q1tU*q2tk1t*P5*m2 + 4*q2(vb)*q1L*
     & q1tU*q2tk1t*P8*m2 + 4*q2(vb)*q1L*q1tU*q2tk2t*P5*m1 + 4*
     & q2(vb)*q1L*q1tU*q2tk2t*P8*m1
      A52(4,vb) = A52(4,vb) - 4*q2(vb)*q2U*q1tk1t*q2tL*P8*m2 - 4*
     & q2(vb)*q2U*q1tk2t*q2tL*P8*m1 + 4*q2(vb)*q2U*q1tL*q2tk1t*P8
     & *m2 + 4*q2(vb)*q2U*q1tL*q2tk2t*P8*m1 - 4*q2(vb)*q2L*q1tk1t*
     & q2tU*P5*m2 - 4*q2(vb)*q2L*q1tk2t*q2tU*P5*m1 + 4*q2(vb)*q2L*
     & q1tU*q2tk1t*P5*m2 + 4*q2(vb)*q2L*q1tU*q2tk2t*P5*m1 + 4*
     & q2(vb)*q1tU*q2tL*P5*m1*q1k2 + 4*q2(vb)*q1tU*q2tL*P5*m1*q2k2
     &  + 4*q2(vb)*q1tU*q2tL*P5*m2*q1k1 + 4*q2(vb)*q1tU*q2tL*P5*m2*
     & q2k1 + 4*q2(vb)*q1tU*q2tL*P8*m1*q1k2 + 4*q2(vb)*q1tU*q2tL*P8
     & *m2*q1k1 - 4*q2(vb)*q1tL*q2tU*P5*m1*q1k2 - 4*q2(vb)*q1tL*
     & q2tU*P5*m1*q2k2 - 4*q2(vb)*q1tL*q2tU*P5*m2*q1k1 - 4*q2(vb)*
     & q1tL*q2tU*P5*m2*q2k1 - 4*q2(vb)*q1tL*q2tU*P8*m1*q1k2 - 4*
     & q2(vb)*q1tL*q2tU*P8*m2*q1k1 - 4*q1t(vb)*q1L*q2U*q2tk1t*P5*
     & m2 + 4*q1t(vb)*q1L*q2U*q2tk1t*P6*m2 + 4*q1t(vb)*q1L*q2U*
     & q2tk1t*P8*m2 - 4*q1t(vb)*q1L*q2U*q2tk2t*P5*m1 + 4*q1t(vb)*
     & q1L*q2U*q2tk2t*P6*m1 + 4*q1t(vb)*q1L*q2U*q2tk2t*P8*m1 - 4*
     & q1t(vb)*q1L*q2tU*P5*m1*q2k2
      A52(4,vb) = A52(4,vb) - 4*q1t(vb)*q1L*q2tU*P5*m2*q2k1 + 4*
     & q1t(vb)*q1L*q2tU*P6*m1*q2k2 + 4*q1t(vb)*q1L*q2tU*P6*m2*q2k1
     &  + 4*q1t(vb)*q1L*q2tU*P8*m1*q2k2 + 4*q1t(vb)*q1L*q2tU*P8*m2*
     & q2k1 - 4*q1t(vb)*q2U*q2L*q2tk1t*P5*m2 + 4*q1t(vb)*q2U*q2L*
     & q2tk1t*P6*m2 + 4*q1t(vb)*q2U*q2L*q2tk1t*P8*m2 - 4*q1t(vb)*
     & q2U*q2L*q2tk2t*P5*m1 + 4*q1t(vb)*q2U*q2L*q2tk2t*P6*m1 + 4*
     & q1t(vb)*q2U*q2L*q2tk2t*P8*m1 - 4*q1t(vb)*q2U*q2tL*P5*m1*
     & q1k2 - 4*q1t(vb)*q2U*q2tL*P5*m1*q2k2 - 4*q1t(vb)*q2U*q2tL*P5
     & *m2*q1k1 - 4*q1t(vb)*q2U*q2tL*P5*m2*q2k1 + 4*q1t(vb)*q2U*
     & q2tL*P6*m1*q1k2 + 4*q1t(vb)*q2U*q2tL*P6*m1*q2k2 + 4*q1t(vb)*
     & q2U*q2tL*P6*m2*q1k1 + 4*q1t(vb)*q2U*q2tL*P6*m2*q2k1 + 4*
     & q1t(vb)*q2U*q2tL*P8*m1*q1k2 + 4*q1t(vb)*q2U*q2tL*P8*m1*q2k2
     &  + 4*q1t(vb)*q2U*q2tL*P8*m2*q1k1 + 4*q1t(vb)*q2U*q2tL*P8*m2*
     & q2k1 + 4*q1t(vb)*q2L*q2tU*P5*m1*q1k2 + 4*q1t(vb)*q2L*q2tU*P5
     & *m2*q1k1 - 4*q1t(vb)*q2L*q2tU*P6*m1*q1k2 - 4*q1t(vb)*q2L*
     & q2tU*P6*m2*q1k1
      A52(4,vb) = A52(4,vb) - 4*q1t(vb)*q2L*q2tU*P8*m1*q1k2 - 4*
     & q1t(vb)*q2L*q2tU*P8*m2*q1k1 + 4*q2t(vb)*q1L*q2U*q1tk1t*P5*
     & m2 - 4*q2t(vb)*q1L*q2U*q1tk1t*P6*m2 + 4*q2t(vb)*q1L*q2U*
     & q1tk1t*P8*m2 + 4*q2t(vb)*q1L*q2U*q1tk2t*P5*m1 - 4*q2t(vb)*
     & q1L*q2U*q1tk2t*P6*m1 + 4*q2t(vb)*q1L*q2U*q1tk2t*P8*m1 + 4*
     & q2t(vb)*q1L*q1tU*P5*m1*q2k2 + 4*q2t(vb)*q1L*q1tU*P5*m2*q2k1
     &  - 4*q2t(vb)*q1L*q1tU*P6*m1*q2k2 - 4*q2t(vb)*q1L*q1tU*P6*m2*
     & q2k1 + 4*q2t(vb)*q1L*q1tU*P8*m1*q2k2 + 4*q2t(vb)*q1L*q1tU*P8
     & *m2*q2k1 + 4*q2t(vb)*q2U*q2L*q1tk1t*P5*m2 - 4*q2t(vb)*q2U*
     & q2L*q1tk1t*P6*m2 + 4*q2t(vb)*q2U*q2L*q1tk1t*P8*m2 + 4*
     & q2t(vb)*q2U*q2L*q1tk2t*P5*m1 - 4*q2t(vb)*q2U*q2L*q1tk2t*P6
     & *m1 + 4*q2t(vb)*q2U*q2L*q1tk2t*P8*m1 + 4*q2t(vb)*q2U*q1tL*
     & P5*m1*q1k2 + 4*q2t(vb)*q2U*q1tL*P5*m1*q2k2 + 4*q2t(vb)*q2U*
     & q1tL*P5*m2*q1k1 + 4*q2t(vb)*q2U*q1tL*P5*m2*q2k1 - 4*q2t(vb)*
     & q2U*q1tL*P6*m1*q1k2 - 4*q2t(vb)*q2U*q1tL*P6*m1*q2k2 - 4*
     & q2t(vb)*q2U*q1tL*P6*m2*q1k1
      A52(4,vb) = A52(4,vb) - 4*q2t(vb)*q2U*q1tL*P6*m2*q2k1 + 4*
     & q2t(vb)*q2U*q1tL*P8*m1*q1k2 + 4*q2t(vb)*q2U*q1tL*P8*m1*q2k2
     &  + 4*q2t(vb)*q2U*q1tL*P8*m2*q1k1 + 4*q2t(vb)*q2U*q1tL*P8*m2*
     & q2k1 - 4*q2t(vb)*q2L*q1tU*P5*m1*q1k2 - 4*q2t(vb)*q2L*q1tU*P5
     & *m2*q1k1 + 4*q2t(vb)*q2L*q1tU*P6*m1*q1k2 + 4*q2t(vb)*q2L*
     & q1tU*P6*m2*q1k1 - 4*q2t(vb)*q2L*q1tU*P8*m1*q1k2 - 4*q2t(vb)*
     & q2L*q1tU*P8*m2*q1k1 - 8*U(vb)*q1L*q1tq2t*q1tk1t*P6*m2 - 8*
     & U(vb)*q1L*q1tq2t*q1tk2t*P6*m1 + 4*U(vb)*q1L*q1tq2t*q2tk1t*
     & P5*m2 + 4*U(vb)*q1L*q1tq2t*q2tk1t*P6*m2 + 4*U(vb)*q1L*
     & q1tq2t*q2tk1t*P8*m2 + 4*U(vb)*q1L*q1tq2t*q2tk2t*P5*m1 + 4*
     & U(vb)*q1L*q1tq2t*q2tk2t*P6*m1 + 4*U(vb)*q1L*q1tq2t*q2tk2t*
     & P8*m1 + 4*U(vb)*q1L*q1tq2t*P5*m1*q2k2 + 4*U(vb)*q1L*q1tq2t*
     & P5*m2*q2k1 + 4*U(vb)*q1L*q1tq2t*P6*m1*q2k2 + 4*U(vb)*q1L*
     & q1tq2t*P6*m2*q2k1 + 4*U(vb)*q1L*q1tq2t*P8*m1*q2k2 + 4*U(vb)*
     & q1L*q1tq2t*P8*m2*q2k1 - 8*U(vb)*q1L*q1tk1t*q2tk1t*P5*m2 - 8
     & *U(vb)*q1L*q1tk1t*q2tk1t*P8*m2
      A52(4,vb) = A52(4,vb) + 8*U(vb)*q1L*q1tk1t*q2tk2t*P6*m2 - 4*
     & U(vb)*q1L*q1tk1t*t2*P5*m2 + 4*U(vb)*q1L*q1tk1t*t2*P6*m2 - 4*
     & U(vb)*q1L*q1tk1t*t2*P8*m2 - 8*U(vb)*q1L*q1tk2t*q2tk1t*P5*m1
     &  - 8*U(vb)*q1L*q1tk2t*q2tk1t*P8*m1 + 8*U(vb)*q1L*q1tk2t*
     & q2tk2t*P6*m1 - 4*U(vb)*q1L*q1tk2t*t2*P5*m1 + 4*U(vb)*q1L*
     & q1tk2t*t2*P6*m1 - 4*U(vb)*q1L*q1tk2t*t2*P8*m1 - 4*U(vb)*q2L*
     & q1tq2t*q1tk1t*P5*m2 - 4*U(vb)*q2L*q1tq2t*q1tk1t*P6*m2 - 4*
     & U(vb)*q2L*q1tq2t*q1tk1t*P8*m2 - 4*U(vb)*q2L*q1tq2t*q1tk2t*
     & P5*m1 - 4*U(vb)*q2L*q1tq2t*q1tk2t*P6*m1 - 4*U(vb)*q2L*
     & q1tq2t*q1tk2t*P8*m1 + 8*U(vb)*q2L*q1tq2t*q2tk1t*P5*m2 + 8*
     & U(vb)*q2L*q1tq2t*q2tk2t*P5*m1 - 4*U(vb)*q2L*q1tq2t*P5*m1*
     & q1k2 - 4*U(vb)*q2L*q1tq2t*P5*m2*q1k1 - 4*U(vb)*q2L*q1tq2t*P6
     & *m1*q1k2 - 4*U(vb)*q2L*q1tq2t*P6*m2*q1k1 - 4*U(vb)*q2L*
     & q1tq2t*P8*m1*q1k2 - 4*U(vb)*q2L*q1tq2t*P8*m2*q1k1 - 8*U(vb)*
     & q2L*q1tk1t*q2tk1t*P5*m2 - 8*U(vb)*q2L*q1tk1t*q2tk2t*P5*m1
     &  + 8*U(vb)*q2L*q1tk2t*q2tk1t*P6*m2
      A52(4,vb) = A52(4,vb) + 8*U(vb)*q2L*q1tk2t*q2tk1t*P8*m2 + 8*
     & U(vb)*q2L*q1tk2t*q2tk2t*P6*m1 + 8*U(vb)*q2L*q1tk2t*q2tk2t*
     & P8*m1 - 4*U(vb)*q2L*q2tk1t*t1*P5*m2 + 4*U(vb)*q2L*q2tk1t*t1*
     & P6*m2 + 4*U(vb)*q2L*q2tk1t*t1*P8*m2 - 4*U(vb)*q2L*q2tk2t*t1*
     & P5*m1 + 4*U(vb)*q2L*q2tk2t*t1*P6*m1 + 4*U(vb)*q2L*q2tk2t*t1*
     & P8*m1 - 4*U(vb)*q1tq2t*q1tL*P5*m1*q2k2 - 4*U(vb)*q1tq2t*q1tL
     & *P5*m2*q2k1 - 8*U(vb)*q1tq2t*q1tL*P6*m1*q1k2 - 4*U(vb)*q1tq2t
     & *q1tL*P6*m1*q2k2 - 8*U(vb)*q1tq2t*q1tL*P6*m2*q1k1 - 4*U(vb)*
     & q1tq2t*q1tL*P6*m2*q2k1 - 4*U(vb)*q1tq2t*q1tL*P8*m1*q2k2 - 4*
     & U(vb)*q1tq2t*q1tL*P8*m2*q2k1 + 4*U(vb)*q1tq2t*q2tL*P5*m1*
     & q1k2 + 8*U(vb)*q1tq2t*q2tL*P5*m1*q2k2 + 4*U(vb)*q1tq2t*q2tL*
     & P5*m2*q1k1 + 8*U(vb)*q1tq2t*q2tL*P5*m2*q2k1 + 4*U(vb)*q1tq2t*
     & q2tL*P6*m1*q1k2 + 4*U(vb)*q1tq2t*q2tL*P6*m2*q1k1 + 4*U(vb)*
     & q1tq2t*q2tL*P8*m1*q1k2 + 4*U(vb)*q1tq2t*q2tL*P8*m2*q1k1 - 4*
     & U(vb)*q1tk1t*q2tL*t1*P6*m2 - 4*U(vb)*q1tk1t*q2tL*t2*P5*m2 - 
     & 8*U(vb)*q1tk1t*q2tL*P5*m1*q2k2
      A52(4,vb) = A52(4,vb) - 4*U(vb)*q1tk1t*q2tL*P5*m2*q1q2 - 4*
     & U(vb)*q1tk1t*q2tL*P6*m2*q1q2 + 8*U(vb)*q1tk1t*q2tL*P6*m2*
     & q1k2 - 4*U(vb)*q1tk1t*q2tL*P8*m2*q1q2 - 4*U(vb)*q1tk2t*q2tL*
     & t1*P6*m1 - 4*U(vb)*q1tk2t*q2tL*t2*P5*m1 - 4*U(vb)*q1tk2t*
     & q2tL*P5*m1*q1q2 + 8*U(vb)*q1tk2t*q2tL*P5*m1*q2k1 - 4*U(vb)*
     & q1tk2t*q2tL*P6*m1*q1q2 + 8*U(vb)*q1tk2t*q2tL*P6*m1*q1k2 + 8*
     & U(vb)*q1tk2t*q2tL*P6*m1*q2k2 + 8*U(vb)*q1tk2t*q2tL*P6*m2*
     & q2k1 - 4*U(vb)*q1tk2t*q2tL*P8*m1*q1q2 + 8*U(vb)*q1tk2t*q2tL*
     & P8*m1*q2k2 + 8*U(vb)*q1tk2t*q2tL*P8*m2*q2k1 + 4*U(vb)*q1tL*
     & q2tk1t*t1*P6*m2 + 4*U(vb)*q1tL*q2tk1t*t2*P5*m2 - 8*U(vb)*
     & q1tL*q2tk1t*P5*m1*q1k2 + 4*U(vb)*q1tL*q2tk1t*P5*m2*q1q2 - 8*
     & U(vb)*q1tL*q2tk1t*P5*m2*q1k1 - 8*U(vb)*q1tL*q2tk1t*P5*m2*
     & q2k1 + 4*U(vb)*q1tL*q2tk1t*P6*m2*q1q2 - 8*U(vb)*q1tL*q2tk1t*
     & P6*m2*q1k2 - 8*U(vb)*q1tL*q2tk1t*P8*m1*q1k2 + 4*U(vb)*q1tL*
     & q2tk1t*P8*m2*q1q2 - 8*U(vb)*q1tL*q2tk1t*P8*m2*q1k1 + 4*U(vb)*
     & q1tL*q2tk2t*t1*P6*m1
      A52(4,vb) = A52(4,vb) + 4*U(vb)*q1tL*q2tk2t*t2*P5*m1 + 4*
     & U(vb)*q1tL*q2tk2t*P5*m1*q1q2 - 8*U(vb)*q1tL*q2tk2t*P5*m1*
     & q2k1 + 4*U(vb)*q1tL*q2tk2t*P6*m1*q1q2 + 8*U(vb)*q1tL*q2tk2t*
     & P6*m2*q1k1 + 4*U(vb)*q1tL*q2tk2t*P8*m1*q1q2 - 4*U(vb)*q1tL*t2
     & *P5*m1*q1k2 - 4*U(vb)*q1tL*t2*P5*m2*q1k1 + 4*U(vb)*q1tL*t2*P6*
     & m1*q1k2 + 4*U(vb)*q1tL*t2*P6*m2*q1k1 - 4*U(vb)*q1tL*t2*P8*m1*
     & q1k2 - 4*U(vb)*q1tL*t2*P8*m2*q1k1 - 4*U(vb)*q2tL*t1*P5*m1*q2k2
     &  - 4*U(vb)*q2tL*t1*P5*m2*q2k1 + 4*U(vb)*q2tL*t1*P6*m1*q2k2 + 4
     & *U(vb)*q2tL*t1*P6*m2*q2k1 + 4*U(vb)*q2tL*t1*P8*m1*q2k2 + 4*
     & U(vb)*q2tL*t1*P8*m2*q2k1 + 4*L(vb)*q2U*q1tq2t*q1tk1t*P5*m2
     &  - 4*L(vb)*q2U*q1tq2t*q1tk1t*P6*m2 + 4*L(vb)*q2U*q1tq2t*
     & q1tk1t*P8*m2 + 4*L(vb)*q2U*q1tq2t*q1tk2t*P5*m1 - 4*L(vb)*
     & q2U*q1tq2t*q1tk2t*P6*m1 + 4*L(vb)*q2U*q1tq2t*q1tk2t*P8*m1
     &  - 4*L(vb)*q2U*q1tq2t*q2tk1t*P5*m2 + 4*L(vb)*q2U*q1tq2t*
     & q2tk1t*P6*m2 + 4*L(vb)*q2U*q1tq2t*q2tk1t*P8*m2 - 4*L(vb)*
     & q2U*q1tq2t*q2tk2t*P5*m1
      A52(4,vb) = A52(4,vb) + 4*L(vb)*q2U*q1tq2t*q2tk2t*P6*m1 + 4*
     & L(vb)*q2U*q1tq2t*q2tk2t*P8*m1 + 4*L(vb)*q2U*q1tq2t*P5*m1*
     & q1k2 + 4*L(vb)*q2U*q1tq2t*P5*m1*q2k2 + 4*L(vb)*q2U*q1tq2t*P5
     & *m2*q1k1 + 4*L(vb)*q2U*q1tq2t*P5*m2*q2k1 + 4*L(vb)*q2U*
     & q1tq2t*P6*m1*q1k2 + 4*L(vb)*q2U*q1tq2t*P6*m1*q2k2 + 4*L(vb)*
     & q2U*q1tq2t*P6*m2*q1k1 + 4*L(vb)*q2U*q1tq2t*P6*m2*q2k1 + 4*
     & L(vb)*q2U*q1tq2t*P8*m1*q1k2 + 4*L(vb)*q2U*q1tq2t*P8*m1*q2k2
     &  + 4*L(vb)*q2U*q1tq2t*P8*m2*q1k1 + 4*L(vb)*q2U*q1tq2t*P8*m2*
     & q2k1 - 8*L(vb)*q2U*q1tk1t*q2tk1t*P8*m2 + 8*L(vb)*q2U*q1tk1t
     & *q2tk2t*P5*m1 + 8*L(vb)*q2U*q1tk1t*q2tk2t*P6*m2 - 4*L(vb)*
     & q2U*q1tk1t*t2*P5*m2 + 4*L(vb)*q2U*q1tk1t*t2*P6*m2 - 4*L(vb)*
     & q2U*q1tk1t*t2*P8*m2 - 8*L(vb)*q2U*q1tk2t*q2tk1t*P5*m1 - 8*
     & L(vb)*q2U*q1tk2t*q2tk1t*P6*m2 - 8*L(vb)*q2U*q1tk2t*q2tk1t*
     & P8*m1 - 8*L(vb)*q2U*q1tk2t*q2tk1t*P8*m2 - 8*L(vb)*q2U*
     & q1tk2t*q2tk2t*P8*m1 - 4*L(vb)*q2U*q1tk2t*t2*P5*m1 + 4*L(vb)*
     & q2U*q1tk2t*t2*P6*m1
      A52(4,vb) = A52(4,vb) - 4*L(vb)*q2U*q1tk2t*t2*P8*m1 + 4*L(vb)
     & *q2U*q2tk1t*t1*P5*m2 - 4*L(vb)*q2U*q2tk1t*t1*P6*m2 - 4*L(vb)
     & *q2U*q2tk1t*t1*P8*m2 + 4*L(vb)*q2U*q2tk2t*t1*P5*m1 - 4*L(vb)
     & *q2U*q2tk2t*t1*P6*m1 - 4*L(vb)*q2U*q2tk2t*t1*P8*m1 + 4*L(vb)
     & *q1tq2t*q1tU*P5*m1*q2k2 + 4*L(vb)*q1tq2t*q1tU*P5*m2*q2k1 + 8
     & *L(vb)*q1tq2t*q1tU*P6*m1*q1k2 + 4*L(vb)*q1tq2t*q1tU*P6*m1*
     & q2k2 + 8*L(vb)*q1tq2t*q1tU*P6*m2*q1k1 + 4*L(vb)*q1tq2t*q1tU*
     & P6*m2*q2k1 + 4*L(vb)*q1tq2t*q1tU*P8*m1*q2k2 + 4*L(vb)*q1tq2t*
     & q1tU*P8*m2*q2k1 - 4*L(vb)*q1tq2t*q2tU*P5*m1*q1k2 - 8*L(vb)*
     & q1tq2t*q2tU*P5*m1*q2k2 - 4*L(vb)*q1tq2t*q2tU*P5*m2*q1k1 - 8*
     & L(vb)*q1tq2t*q2tU*P5*m2*q2k1 - 4*L(vb)*q1tq2t*q2tU*P6*m1*
     & q1k2 - 4*L(vb)*q1tq2t*q2tU*P6*m2*q1k1 - 4*L(vb)*q1tq2t*q2tU*
     & P8*m1*q1k2 - 4*L(vb)*q1tq2t*q2tU*P8*m2*q1k1 + 4*L(vb)*q1tk1t*
     & q2tU*t1*P6*m2 + 4*L(vb)*q1tk1t*q2tU*t2*P5*m2 + 8*L(vb)*
     & q1tk1t*q2tU*P5*m1*q2k2 + 4*L(vb)*q1tk1t*q2tU*P5*m2*q1q2 + 4*
     & L(vb)*q1tk1t*q2tU*P6*m2*q1q2
      A52(4,vb) = A52(4,vb) - 8*L(vb)*q1tk1t*q2tU*P6*m2*q1k2 + 4*
     & L(vb)*q1tk1t*q2tU*P8*m2*q1q2 + 4*L(vb)*q1tk2t*q2tU*t1*P6*m1
     &  + 4*L(vb)*q1tk2t*q2tU*t2*P5*m1 + 4*L(vb)*q1tk2t*q2tU*P5*m1*
     & q1q2 - 8*L(vb)*q1tk2t*q2tU*P5*m1*q2k1 + 4*L(vb)*q1tk2t*q2tU*
     & P6*m1*q1q2 - 8*L(vb)*q1tk2t*q2tU*P6*m1*q1k2 - 8*L(vb)*q1tk2t*
     & q2tU*P6*m1*q2k2 - 8*L(vb)*q1tk2t*q2tU*P6*m2*q2k1 + 4*L(vb)*
     & q1tk2t*q2tU*P8*m1*q1q2 - 8*L(vb)*q1tk2t*q2tU*P8*m1*q2k2 - 8*
     & L(vb)*q1tk2t*q2tU*P8*m2*q2k1 - 4*L(vb)*q1tU*q2tk1t*t1*P6*m2
     &  - 4*L(vb)*q1tU*q2tk1t*t2*P5*m2 + 8*L(vb)*q1tU*q2tk1t*P5*m1*
     & q1k2 - 4*L(vb)*q1tU*q2tk1t*P5*m2*q1q2 + 8*L(vb)*q1tU*q2tk1t*
     & P5*m2*q1k1 + 8*L(vb)*q1tU*q2tk1t*P5*m2*q2k1 - 4*L(vb)*q1tU*
     & q2tk1t*P6*m2*q1q2 + 8*L(vb)*q1tU*q2tk1t*P6*m2*q1k2 + 8*L(vb)*
     & q1tU*q2tk1t*P8*m1*q1k2 - 4*L(vb)*q1tU*q2tk1t*P8*m2*q1q2 + 8*
     & L(vb)*q1tU*q2tk1t*P8*m2*q1k1 - 4*L(vb)*q1tU*q2tk2t*t1*P6*m1
     &  - 4*L(vb)*q1tU*q2tk2t*t2*P5*m1 - 4*L(vb)*q1tU*q2tk2t*P5*m1*
     & q1q2
      A52(4,vb) = A52(4,vb) + 8*L(vb)*q1tU*q2tk2t*P5*m1*q2k1 - 4*
     & L(vb)*q1tU*q2tk2t*P6*m1*q1q2 - 8*L(vb)*q1tU*q2tk2t*P6*m2*
     & q1k1 - 4*L(vb)*q1tU*q2tk2t*P8*m1*q1q2 + 4*L(vb)*q1tU*t2*P5*m1
     & *q1k2 + 4*L(vb)*q1tU*t2*P5*m2*q1k1 - 4*L(vb)*q1tU*t2*P6*m1*
     & q1k2 - 4*L(vb)*q1tU*t2*P6*m2*q1k1 + 4*L(vb)*q1tU*t2*P8*m1*q1k2
     &  + 4*L(vb)*q1tU*t2*P8*m2*q1k1 + 4*L(vb)*q2tU*t1*P5*m1*q2k2 + 4
     & *L(vb)*q2tU*t1*P5*m2*q2k1 - 4*L(vb)*q2tU*t1*P6*m1*q2k2 - 4*
     & L(vb)*q2tU*t1*P6*m2*q2k1 - 4*L(vb)*q2tU*t1*P8*m1*q2k2 - 4*
     & L(vb)*q2tU*t1*P8*m2*q2k1

      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccc

      do kk=1,4
        do j=1,4
          do i=1,4
	    A1=-2*A8(kk,i)-AS1(kk,i)+AS2(kk,i)
	    B1=-2*A8(kk,j)-AS1(kk,j)+AS2(kk,j)
	    A2=AS1(kk,i)+AS2(kk,i)
	    B2=AS1(kk,j)+AS2(kk,j)
	    A15=-2*A58(kk,i)-A51(kk,i)+A52(kk,i)
	    B15=-2*A58(kk,j)-A51(kk,j)+A52(kk,j)
	    A25=A51(kk,i)+A52(kk,i)
	    B25=A51(kk,j)+A52(kk,j)
            M(kk,i,j)=(3.*A1*B1+7./3.*A2*B2)/dabs(N0(kk))
            M5(kk,i,j)=(3.*A15*B15+7./3.*A25*B25)/dabs(N5(kk))
          enddo
        enddo
      enddo

c      write(6,*) '*'
c      write(6,*) N0
c      write(6,*) N5

      do i=1,4
        do j=1,4
          d(i,j)=0.
	enddo
      enddo

      dd(4,4)=-1.
      dd(3,3)=1.
      dd(2,2)=1.
      dd(1,1)=1.

c      write(6,*) d(1,1),d(2,2)

      do i=1,4
        do j=1,4
          d(i,j)=dd(i,j)+dd(i,i)*dd(j,j)*P(c,i)*P(c,j)/mw**2
	enddo
      enddo

      ME=0.

      do kk=1,4
        do i=1,4
          do j=1,4
            ME=ME+vq**2*M(kk,i,j)*d(i,j)+aq**2*M5(kk,i,j)*d(i,j)
	  enddo
	enddo
      enddo

      return
      end
