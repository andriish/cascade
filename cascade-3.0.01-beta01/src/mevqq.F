      SUBROUTINE mevqq(WT1)
      IMPLICIT None
C
c
C    GLUON GLUON ----> Q Q_BAR +W 
C
c use the momenta ala CCH
c
c   P1 (Kin1)
c   ---------------------------
c      \
c    K1 \ KIA1     / P4 (KF1)
c        \        /
c         ----------------P5(W/Z)
c        /        \
c    K2 / KIA2     \ P3 (KF2)
c      / 
c   --------------------------
c   P2 (Kin2)
c
c

C

#include "caslujets.inc"
#include "caspara1.inc"
#include "caspara.inc"
#include "casluco.inc"
#include "casludat2.inc"
#include "caspara5.inc"
#include "caspart.inc"
#include "casgki.inc"
#include "casshwr.inc"
#include "casmult.inc"
#include "casscal.inc"
      DOUBLE PRECISION SIN2W,XMW2
      COMMON/ELWEAK/SIN2W,XMW2
	DOUBLE PRECISION XMW
      Integer IGENFL,IMEPHI
      COMMON/GENWEI/IGENFL,IMEPHI
      DOUBLE PRECISION K1K2,K1P2,K1P3,K1P4,K2P1,K2P3,K2P4,K12,K22
      DOUBLE PRECISION P1P2,P1P3,P1P4,P2P3,P2P4,M2,PT2Q
      Double Precision k1p3t,k12t
      Double Precision XPQ1(-6:6),XPQ2(-6:6),SFSCAL1,SFSCAL2 
      Double Precision STMAX1,STMAX2
      common/myscale/STMAX1,STMAX2
      DOUBLE PRECISION DSHN(4),DSH(4),dshnf(4),pfin(4)
      Double Precision DScale,alph_em
      Double Precision scale
      Double Precision CUT
c      DATA CUT /1.d-10/
c      DATA CUT /1.d-3/
      DATA CUT /1.d-6/ ! as smallx
      Double Precision Sfcut
      DATA Sfcut /1.d-3/ 
c      DATA Sfcut /1.d-9/ 
	
      Integer I,KIA1,KIA2,KF1,KF2,Kin1,Kin2
      Integer KPZ,ud
      Double Precision QF2,alpha_s,qg2,qg1,sh,th,uh,wt1,sum2,wtg
      Double Precision cab,cabB,SH_NEW,pcm_fin,cos_fin,phi,phio,test
      Double Precision sh_old,shf,check,dphi,dtheta
      Double Precision E1,E2,Aab,Anab,yph,yph_bar,xsi_hard(2),aL,epsilon
      Double Precision k1q,k1q1,k1p,kq,pq,pq1,plus,minus
      Double Precision X1,X2
      Double Precision gw2,gs2,als,XME,XM1,XM2
      Double Precision pcm_new,alphas,nq

      double Precision Aas

      Double Precision dot,dot1,dcasrn
      Double Precision sL,rt,ru
      Double Precision Del,Nc
      Double Precision rap1,rap2,rapp,ptstm1,ptstm2,zz1,zz2
      Integer npfin,kpfl
      Integer IQBAR
      Integer  ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
	Double Precision ME,XMTOT,rntest
      Double Precision Vqq
      Integer KW1,KW2
      
	Double Precision QQg,lambda
	Double Precision Qscal
	Integer Iqqbar
	Common/GLUDAT2/Qscal,Iqqbar
c      Data IQBAR/1/
c      Data IQBAR/2/
      Double Precision PT2
      Integer Nacc

      Data Nc/3.d0/
      Logical ltest
c      data ltest/.true./
      data ltest/.false./
      Logical first
      data first/.true./
c      SAVE STMAX1,STMAX2, 
      Save XPQ2, phio   

c     call pylist(1)
c KPZ=23 for Z
c KPZ=24 for W-+
c Z case:
c ud=0   for u-type
c ud=1   for d-type
      KPZ=K(NF1+1,2)
      M2 = P(NF1+1,5)**2
      ud=0
      
      IQBAR = Iqqbar
c      IQbar = 2
      NPFIN = N
c	call pylist(1)
      WT1 = 0.d0
c      write(6,*) ' mevqq ',am(1),am(2)
c	call pylist(1)
c boost  cms
      CALL PYROBO(0,N,0.D0,0.D0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
      KIA1=N+1
      KIA2=N+2
      KF1=N+3
      KF2=N+5
      Kin1=N+6
      Kin2=N+7
      N=N+7
      DO I=1,4
         P(KIA1,I) = P(NIA1,I)
         P(KIA2,I) = P(NIA2,I)
         P(KF1,I) = P(NF1,I)
         P(KF1+1,I) = P(NF1+1,I)
         P(KF2,I) = P(NF2,I)
         K(KF1+1,I) = K(NF1+1,I)
         DSH(I)= P(NIA1,I) + P(NIA2,I)
      ENDDO
      K(KIA1,1)=13
      K(KIA2,1)=13
      K(KF1,1)=13
      K(KF2,1)=13
      K(KF1+1,1)=13

      K(Kin1,1)=13
      K(Kin2,1)=13
c define incoming particles as massless, to ensure proper cancellation in
c matrix-element
      P(Kin1,4)=Dabs(P(1,3))
      P(Kin1,5)=0
      P(Kin2,4)=Dabs(P(2,3))
      P(Kin2,5)=0
      Do I=1,3
         P(Kin1,I)=P(1,I)
         P(Kin2,I)=P(2,I)
      enddo
c
      SH_old = DOT(DSH,DSH)

      if(ltest) then
         SH_NEW = DOT(DSH,DSH)
         goto 10
      endif
c      call pylist(1)
      CALL PYROBO(0,0,0.D0,0.D0,-DSH(1)/DSH(4), -DSH(2)/DSH(4), -DSH(3)
     +/DSH(4))
      DPHI = PYANGL(P(NIA1,1),P(NIA1,2))
      call PYROBO(0,0,0.D0,-dphi,0.d0,0.d0,0.d0)
      DTHETA = PYANGL(P(NIA1,3),P(NIA1,1))
      call PYROBO(0,0,-DTHETA,0.D0,0.d0,0.d0,0.d0)

      K12 = DOT1(KIA1,KIA1)
      K22 = DOT1(KIA2,KIA2)
c      IF(ICCFM.EQ.0) Then
c	   If(ires(1).eq.1) then
c            P(KIA1,1) = 0D0
c            P(KIA1,2) = 0D0
c            P(KIA1,5) = 0
c            K12 = 0
c         endif
c	   If(ires(2).eq.1) then
c            P(KIA2,1) = 0D0
c            P(KIA2,2) = 0D0
c            P(KIA2,5) = 0
c            K22 = 0
c         endif
c
c      ENDIF
c      P(KIA1,3) = dsqrt(((sh_old-K12-K22)**2 - 4d0*K12*K22)/4d0/sh_old)
c      P(KIA1,4) = 0.5d0*(sh_old+(-k22+k12))/sqrt(sh_old)
c      P(KIA2,3) =-dsqrt(((sh_old-K12-K22)**2 - 4d0*K12*K22)/4d0/sh_old)
c      P(KIA2,4) = 0.5d0*(sh_old-(-k22+k12))/sqrt(sh_old) 

c      write(6,*) ' Kia1,2 = ',KIA1,KIA2,' KF1,2 = ',KF1,KF2,KF1+1
c      call pylist(1)
            
      do i=1,4
      dshnf(i) = P(kf1,i)+P(kf2,i)+P(kf1+1,i)
c      dshnf(i) = P(kia1,i)+P(kia2,i) 
      enddo
      shf = dot(dshnf,dshnf)

      check = abs(P(kf1,1)+P(kf1+1,1)+P(kf2,1)+
     &  P(kf1,2)+P(kf1+1,2)+P(kf2,2) +
     &  P(kf1,3)+P(kf1+1,3)+P(kf2,3) +
     &  P(kf1,4)+P(kf1+1,4)+P(kf2,4)-P(KIA2,4)-P(KIA1,4))
      if(check.gt.1e-5) then 
      write(6,*) 'shat ',sqrt(sh_old),sqrt(shf),check
      write(6,*) ' mevqq momentum checks ',P(kf1,1)+P(kf1+1,1)+P(kf2,1),
     &  P(kf1,2)+P(kf1+1,2)+P(kf2,2) +
     &  P(kf1,3)+P(kf1+1,3)+P(kf2,3) +
     &  P(kf1,4)+P(kf1+1,4)+P(kf2,4)-P(KIA2,4)-P(KIA1,4)
      endif
ctest
      PCM_FIN=SQRT(P(NF1,4)**2 - AM(1)**2)
      COS_FIN=P(NF1,3)/PCM_FIN
cend     
      call PYROBO(0,0,DTHETA,0.D0,0.d0,0.d0,0.d0)
      call PYROBO(0,0,0.D0,dphi,0.d0,0.d0,0.d0)

      CALL PYROBO(0,0,0.D0,0.D0,DSH(1)/DSH(4) ,DSH(2)/DSH(4), DSH(3)/
     +DSH(4))

cccc do not perform this rotation
      goto 345      
ctest
cc      P(KIA1,3) = 0.5d0*(P(KIA1,3)-P(KIA1,4))
cc      P(KIA1,4) = DABS(P(KIA1,3))
cc      P(KIA2,3) = 0.5d0*(P(KIA2,3)+P(KIA2,4))
cc      P(KIA2,4) = DABS(P(KIA2,3))
cc       DO I=1,4
cc         DSHN(I)= P(KIA1,I) + P(KIA2,I)
cc      ENDDO
cc      SH_NEW = DOT(DSHN,DSHN)
cc      write(6,*) ' sh_new ', sh_new, sh_old
cc      if(sh_new.lt.0) goto 20
cc     
cc      CALL PYROBO(0,0,0.D0,0.D0,-DSHN(1)/DSHN(4), -DSHN(2)/DSHN(4),
cc     +-DSHN(3)/DSHN(4))
cc      
cc      P(KF1,4) = 0.5d0*DSQRT(SH_NEW)
cc      if(P(KF1,4)**2.le. AM(1)**2) goto 20
cc      PCM_NEW=SQRT(P(KF1,4)**2 - AM(1)**2)
cc      P(KF1,3) =PCM_NEW*COS_FIN
cc      IF(IMEPHI.EQ.0) THEN
cc         phi=cos(2.d0*PI*dcasrn())
cc         phio = phi
cc      ELSEIF(IMEPHI.EQ.1) THEN
cc         phi=phio
cc      ENDIF
cc      P(KF1,1) = DSQRT(PCM_NEW**2-P(KF1,3)**2)*cos(phi)
cc      P(KF1,2) = DSQRT(PCM_NEW**2-P(KF1,3)**2)*sin(phi) 
cc      DO I=1,4
cc         P(KF2,I) = P(KIA1,I) + P(KIA2,I) - P(KF1,I)
cc      ENDDO
cc
cc      CALL PYROBO(0,0,0.D0,0.D0,DSHN(1)/DSHN(4) ,DSHN(2)/DSHN(4),
cc     +DSHN(3)/DSHN(4))
     
345   Continue     
   10 continue
      K12 = DOT1(KIA1,KIA1)
      K22 = DOT1(KIA2,KIA2)
c check on scale:
      PT2Q = (P(KF1,1)**2 +P(KF1,2)**2+P(KF2,1)**2 +P(KF2,2)**2)/2.
      SCALE = PT2Q+M2
      IF(IQ2.EQ.1) THEN
         SCALE = 4.D0*M2
      ELSEIF(IQ2.EQ.2) THEN
         SCALE = SH_NEW
      ELSEIF(IQ2.EQ.3) THEN
         SCALE = M2 + PT2Q
      ELSEIF(IQ2.EQ.6) THEN
         SCALE = DABS(K12)+DABS(K22)
      ELSE
         WRITE(6,*) ' NO VALID Q2 SCALE. STOP'
         STOP
      ENDIF
	
c	write(6,*) ' scales stmax = ',stmax,' sqrt(scale) = ',sqrt(scale)
c	stmax = sqrt(scale)*scalfa
      Dscale = scale
      alpha_s = ALPHAS(sqrt(dscale))
      ALPH_EM = ALPH
c	alph_em = 1./137.
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)
      gw2=alph_em*4.D0*pi/SIN2W
      gs2=alpha_s*4.D0*pi
      if(IABS(KPZ).eq.23.and.IPRO.eq.503) then
        ud = 0 ! up type quarks: u
        XM1 = ME(KIA1,KIA2,KF1,KF2,kpz,ud)
        ud = 1 ! down type quarks: d,s
        XM2 = ME(KIA1,KIA2,KF1,KF2,kpz,ud)
        XME =4.D0/256D0*gw2*gs2**2*(XM1+2.*XM2)/K12/K22
        xmtot = XM1+2.*XM2
        rntest = dcasrn()
        if((xm1+xm2)/xmtot.le.rntest) then
          KPFL = 3
        elseif(xm1/xmtot.le.rntest) then
          KPFL = 1
        else
          KPFL = 2
        endif  
        if(dcasrn().gt.0.5) then
          K(NF1,2)= KPFL
          K(NF2,2)= -KPFL
        else
          K(NF1,2)= -KPFL
          K(NF2,2)= KPFL
        endif    
      elseif(IABS(KPZ).eq.23.and.IPRO.ge.504) then
        IF(iabs(K(NF1,2)).eq.4.or.iabs(K(NF1,2)).eq.6) then 
          ud=0
          elseif(iabs(K(NF1,2)).eq.5) then
          ud=1
          endif
        XM1 = ME(KIA1,KIA2,KF1,KF2,kpz,ud)
        XME =4.D0/256D0*gw2*gs2**2*XM1/K12/K22
      elseif(IABS(KPZ).eq.24) then
c        write(6,*) ' flavors ',K(NF1+1,2),K(NF1,2),K(NF2,2)
        KW1=Abs(K(NF1,2))
        KW2=Abs(K(NF2,2))
        Vqq=VCKM((KW1+1)/2,(KW2+1)/2)
cc        write(6,*) ' mevqq ',KW1,KW2,Vqq,(KW1+1)/2,(KW2+1)/2
        XM1 = ME(KIA1,KIA2,KF1,KF2,kpz,ud)
        XME =4.D0/256D0*gw2*gs2**2*XM1/K12/K22 * Vqq
      endif

      IF(FIRST) THEN
	   XG2SF = 0.001
	   SFSCAL2 = 2.
	   STMAX2 = 5.
         call cauniglu(K(2,2),XG2SF,SFSCAL2,STMAX2,XPQ2)
         IQBAR = Iqqbar
         IF(IQBAR.EQ.1) THEN
            write(6,*) ' mevqq - scale: qq-pair '
         ELSE
            write(6,*) ' mevqq - scale: q or q_bar '
         ENDIF
	      write(6,*) ' mevqq - scale for alphas ',Qscal
         FIRST=.FALSE.
      ENDIF

c	write(6,*) K12,k22,K(nia1,2),nia1,K(nia1,2),nia2
c	call pylist(1)
      QG1 = P(NIA1,1)**2+ P(NIA1,2)**2
      QG2 = P(NIA2,1)**2+ P(NIA2,2)**2
      SFSCAL1 = DABS(QG1)
      SFSCAL2 = DABS(QG2)
c      write(6,*) 'k12 ',k12,' k22 ',k22
c      write(6,*) 'qg1 ',qg1,' qg2 ',qg2
c new scale
      PT2Q = (P(KF1,1)**2 +P(KF1,2)**2+P(KF2,1)**2 +P(KF2,2)**2)/2.
cc      PT2Q = (P(KF1,1)**2 +P(KF1,2)**2)
      xg_bar =  -qg2/(1.d0-xg)/SSS
      yph = dble(YY)
      yph_bar =  -qg1/(1.d0 - yph)/SSS

	if(sfscal1.le.sfcut) sfscal1=sfcut 
	if(sfscal2.le.sfcut) sfscal2=sfcut 
      xsi_hard(1) = (xg+yph_bar)/(yph+xg_bar)	
      xsi_hard(2) = (yph+xg_bar)/(xg+yph_bar)
      if(xsi_hard(2).le.0) then
c            write(6,*) 'mevqq xsi ',xsi_hard(2)
c            write(6,*) 'k12 ',k12,' k22 ',k22
c            write(6,*) ' new xg,xg_bar ',xg,xg_bar
c            write(6,*) ' new yph,yph_bar ',yph,yph_bar
      endif

c which of the quarks is closer to the proton:
c      rapp = PYP(NPFIN,17)
      rap1 = rapp-PYP(NF1,17)
      rap2 = rapp-PYP(NF2,17)
c      zz = q2/yph/sss/xg
      zz1 = dabs(qg1)/yph/sss/xg
      zz2 = dabs(qg2)/yph/sss/xg
c      write(6,*) ' old zz1,2 ',zz1,zz2
      zz1 = 1. - (P(NF1,4)+ P(NF1,3))/2./P(2,4)/yph
      zz2 = 1. - (P(NF2,4)+ P(NF2,3))/2./P(2,4)/xg
c      write(6,*) ' new zz1,2 ',zz1,zz2
c      if(rap1.le.rap2) then
c         ptstm=dsqrt(dabs(qg1)+P(NF1,5)**2)
c      else
c         ptstm=dsqrt(dabs(qg2)+P(NF2,5)**2)
c      endif
* new
      ptstm1=dsqrt(dabs(qg1)+P(NF1,5)**2)
      ptstm2=dsqrt(dabs(qg2)+P(NF2,5)**2)

      If(IQBAR.eq.1) then
         stmax1 = yph*sqrt(xsi_hard(1)*SSS)
         stmax2 = xg*sqrt(xsi_hard(2)*SSS)
      elseif(IQBAR.EQ.2) then
         stmax1 = ptstm1/(1.-zz1)
         stmax2 = ptstm2/(1.-zz2)
      else
         write(6,*) ' mevqq. no valid scale selected ',IQBAR
      endif
      stmax1 = stmax1 *scalfaf
      stmax2 = stmax2 *scalfaf
	
c      stmax1 = 2 * stmax1	
c      stmax2 = 2 * stmax2
c            write(6,*) 'mevqq xsi ',xsi_hard(2),' stmax ',stmax
c            write(6,*) ' new xg,xg_bar ',xg,xg_bar
c            write(6,*) ' new yph,yph_bar ',yph,yph_bar
c      write(6,*) ' old newstmax2 ',xg*sqrt(xsi_hard(2)*SSS),stmax2
      CALL PYROBO(0,N,0.D0,0.D0,-CM(1)/CM(4),-CM(2)/CM(4),-CM(3)/CM(4)
     +)

c      write(6,*) ' mevqq new alphas ',alpha_s,pt2q
       
      SUM2 = XME
      if(stmax1.ge.sqrt(sss)) then
            XPQ1(0) = 0.d0
      elseif(stmax2.ge.sqrt(sss)) then
            XPQ2(0) = 0.d0
      else
            call cauniglu(K(1,2),XG1SF,SFSCAL1,STMAX1,XPQ1)
            call cauniglu(K(2,2),XG2SF,SFSCAL2,STMAX2,XPQ2)
      endif
      WTG= dble(XPQ1(0)*XPQ2(0))


C new cut for on-shell DGLAP
	If(Iccfm.eq.0) then
c	   If(ires(1).eq.1.and.(Qg1.gt.(PT2Q+M2))) SUM2 = 0 
	   If(ires(1).eq.1.and.(Qg1.gt.scale*scalfa)) SUM2 = 0 
	   If(ires(1).eq.1.and.(Qg1.gt.STMAX1**2)) SUM2 = 0 
c	   If(ires(2).eq.1.and.(Qg2.gt.(PT2Q+M2))) SUM2 = 0 
	   If(ires(2).eq.1.and.(Qg2.gt.scale*scalfa)) SUM2 = 0 
	   If(ires(2).eq.1.and.(Qg2.gt.STMAX2**2)) SUM2 = 0 
	Endif

c      Write(6,*) ' mevqq Cab,Aab,CabB,Anab ',Cab,Aab,CabB,Anab
c      Write(6,*) ' mevqq Wtg ',Wtg,XG1SF,XG2SF
c      Write(6,*) ' mevqq Cab,Aab,CabB,Anab ',Cab,Aab,CabB,Anab
      IF(SUM2.LT.0D0) THEN
c     Write(6,*) ' mevqq Cab,Aab,CabB,Anab ',Cab,Aab,CabB,Anab
c     write(6,*) ' SUM2 ',SUM2
         if(sum2.ne.sum2) then
            write(6,*) ' mevqq SUM2,xsi_hard ',SUM2,
     &        xsi_hard(1),xsi_hard(2)
            Write(6,*) ' mevqq Cab,Aab,CabB,Anab ',Cab,Aab,CabB,
     +      Anab
            write(6,*) ' mevqq k12,k22 ',k12,k22
            write(6,*) ' mevqq th,uh ',th,uh
         endif
         SUM2 = 0D0
      ENDIF
c      write(6,*) sum2,wtg,xme
      SUM2 = SUM2*WTG
      
      
c apply pt cuts here
      PT2 = PYP(KF1,9)
      CALL CUTG(PT2,NACC)
      IF(NACC.EQ.0) SUM2 = 0
      PT2 = PYP(KF2,9)
      CALL CUTG(PT2,NACC)
      IF(NACC.EQ.0) SUM2 = 0

      

      WT1 = SUM2

      if(sum2.ne.sum2) then

         write(6,*) ' mevqq xpq(0),alpha_s',xpq1(0),xpq2(0),
     +      alpha_s
         write(6,*) ' mevqq SUM2,xsi_hard ',SUM2,
     &      xsi_hard(1),xsi_hard(2)
         Write(6,*) ' mevqq Cab,Aab,CabB,Anab ',Cab,Aab,CabB,Anab
         write(6,*) ' mevqq k12,k22 ',k12,k22
         write(6,*) ' mevqq th,uh ',th,uh
         write(6,*) ' P(KIA1) ',(P(KIA1,I),I=1,5)
         write(6,*) ' P(KIA2) ',(P(KIA2,I),I=1,5)
         write(6,*) ' P(KF1) ',(P(KF1,I),I=1,5)
         write(6,*) ' P(KF2) ',(P(KF2,I),I=1,5)
         call pylist(1)
         stop

      endif
      CALL PYEDIT(13)
      RETURN
   20 WT1 = 0.d0
      CALL PYEDIT(13)
      RETURN
      END
