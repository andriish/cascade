      SUBROUTINE meoffsh(WT1)
      IMPLICIT None
C
c
c use the momenta ala CCH
c
c   P1 (Kin1)
c   ---------------------------
c      \
c    K1 \ KIA1     / P4 (KF1)
c        \        /
c         --------
c        /        \
c    K2 / KIA2     \ P3 (KF2)
c      / 
c   --------------------------
c   P2 (Kin2)
c
c
C updated 16.Feb 2013
C ME checked for gg->gg and qg->qg (including prefactors)
C

#include "caslujets.inc"
#include "caspara1.inc"
#include "caspara.inc"
#include "casluco.inc"
#include "caspara5.inc"
#include "caspart.inc"
#include "casgki.inc"
#include "casshwr.inc"
#include "casmult.inc"
#include "casscal.inc"
      Integer IGENFL,IMEPHI
      COMMON/GENWEI/IGENFL,IMEPHI
      DOUBLE PRECISION K1K2,K1P2,K1P3,K1P4,K2P1,K2P3,K2P4,K12,K22
      DOUBLE PRECISION P1P2,P1P3,P1P4,P2P3,P2P4,P3P4,M2,PT2Q
      Double Precision k1p3t,k12t,qq1
      Double Precision XPQ1(-6:6),XPQ2(-6:6),SFSCAL1,SFSCAL2 
      Double Precision STMAX1,STMAX2
      common/myscale/STMAX1,STMAX2
      DOUBLE PRECISION DSHN(4),DSH(4),dshnf(4),pfin(4)
      Double Precision DScale
      Double Precision scale
      Double Precision CUT
c      DATA CUT /1.d-10/
c      DATA CUT /1.d-3/
      DATA CUT /1.d-6/ ! as smallx
      Double Precision Sfcut
      DATA Sfcut /1.d-3/ 
c      DATA Sfcut /1.d-9/ 
	
      Integer I,KIA1,KIA2,KF1,KF2,Kin1,Kin2
      Double Precision QF2,alpha_s,qg2,qg1,sh,th,uh,wt1,sum2
      Double Precision wtg,wtq1,wtq2
      Double Precision cab,cabB,SH_NEW,pcm_fin,cos_fin,phi,phio,test
      Double Precision sh_old,shf,check,dphi,dtheta
      Double Precision E1,E2,yph,yph_bar,xsi_hard(2),epsilon
      Double Precision k1q,k1q1,k1p,kq,pq,pq1,plus,minus
      Double Precision X1,X2
      Double Precision alph_em,pcm_new,alphas,nq
      Double Precision SumAab,SumAnab,SumaL
      Double Precision SumB1ab,SumB1nab,SumB2ab,SumB2nab
      Double Precision SumC1ab,SumC1nab,SumC2ab,SumC2nab
      Double Precision SUMA,SUMB,SUMC,SUMC1,SUMC2
      Double Precision A1,A2,A3,Asum
      double Precision XPQC1,XPQC2,QSUM
	Double Precision rntest,Rleg,Rlego
      Integer NFLT
      Integer ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON /CASCOL/ ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      double Precision Aas
      Integer IRPAT,IRPBT,IRPCT,IRPDT,IRPET,IRPFT,IRPGT

      Double Precision dot,dot1,dcasrn
      Double Precision sL,rt,ru
      Double Precision Del,Nc
      Double Precision rap1,rap2,rapp,ptstm1,ptstm2,zz1,zz2,QT2
      Integer npfin
      Integer IQBAR
      Integer  ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
	Double Precision MEQG
      
	Double Precision QQg,lambda
	Double Precision Qscal
	Integer Iqqbar
	Common/GLUDAT2/Qscal,Iqqbar
c      Data IQBAR/1/
c      Data IQBAR/2/
      Double Precision PT2
      Integer Nacc

      Integer Icoll
      Common/Colltest/Icoll

      Data Nc/3.d0/
      Logical ltest
c      data ltest/.true./
      data ltest/.false./
      Logical first
      data first/.true./
c      SAVE STMAX1,STMAX2, 
      Save XPQ2, phio   

c     call pylist(1)
      IQBAR = Iqqbar
      IRESPRO = 0
      IRPAT = IRPA
      IRPBT = IRPB
      IRPCT = IRPC
      IRPDT = IRPD
      IRPET = IRPE
      IRPGT = IRPG
c      IQbar = 2
      IF(INTER.EQ.2) THEN
         write(6,*) ' INTER = 2 not implemented for meoffsh'
         STOP
      ENDIF
      
c      write(6,*) IPRO,IRPA,IRPB,IRPC
      IF(IPRO.EQ.11) THEN
c	   write(6,*) kpa
         QF2 = DFLOAT(PYCHGE(KPA))**2/9.D0
         NQ = 1d0
         IRPBT = 0
         IRPCT = 0
         IRPDT = 0
         IRPET = 0
      ELSE
         QF2=2.D0/3.D0
         NQ = 3d0
      ENDIF
      NPFIN = N
c	call pylist(1)
      SH = 2.D0 * AM(1)**2 + 2.D0 * DOT1(NF1,NF2)
      WT1 = 0.d0
      IF(SH.LE.0.) RETURN
c boost  cms
      CALL PYROBO(0,N,0.D0,0.D0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
      KIA1=N+1
      KIA2=N+2
      KF1=N+3
      KF2=N+4
      Kin1=N+5
      Kin2=N+6
      N=N+6
      DO I=1,4
         K(KIA1,I) = K(NIA1,I)
         K(KIA2,I) = K(NIA2,I)
         K(KF1,I) = K(NF1,I)
         K(KF2,I) = K(NF2,I)
         P(KIA1,I) = P(NIA1,I)
         P(KIA2,I) = P(NIA2,I)
         P(KF1,I) = P(NF1,I)
         P(KF2,I) = P(NF2,I)
         DSH(I)= P(NIA1,I) + P(NIA2,I)
      ENDDO
      K(KIA1,1)=13
      K(KIA2,1)=13
      K(KF1,1)=13
      K(KF2,1)=13

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

      K12 = DOT1(KIA1,KIA1)
      K22 = DOT1(KIA2,KIA2)
      IF(ICCFM.EQ.0) Then
	   If(ires(1).eq.1) then
            P(KIA1,1) = 0D0
            P(KIA1,2) = 0D0
            P(KIA1,5) = 0
            K12 = 0
         endif
	   If(ires(2).eq.1) then
            P(KIA2,1) = 0D0
            P(KIA2,2) = 0D0
            P(KIA2,5) = 0
            K22 = 0
         endif

      Else
         CALL PYROBO(0,0,0.D0,0.D0,-DSH(1)/DSH(4), -DSH(2)/DSH(4), -DSH(3)
     +   /DSH(4))
         DPHI = PYANGL(P(KIA1,1),P(KIA1,2))
         call PYROBO(0,0,0.D0,-dphi,0.d0,0.d0,0.d0)
         DTHETA = PYANGL(P(KIA1,3),P(KIA1,1))
         call PYROBO(0,0,-DTHETA,0.D0,0.d0,0.d0,0.d0)

         P(KIA1,3) = dsqrt(((sh_old-K12-K22)**2 - 4d0*K12*K22)/4d0/sh_old)
         P(KIA1,4) = 0.5d0*(sh_old+(-k22+k12))/sqrt(sh_old)
         P(KIA2,3) =-dsqrt(((sh_old-K12-K22)**2 - 4d0*K12*K22)/4d0/sh_old)
         P(KIA2,4) = 0.5d0*(sh_old-(-k22+k12))/sqrt(sh_old) 


            
         do i=1,4
            dshnf(i) = P(kf1,i)+P(kf2,i)
         enddo
         shf = dot(dshnf,dshnf)

         check = abs(P(kf1,1)+P(kf2,1)+
     &     P(kf1,2)+P(kf2,2) +
     &     P(kf1,3)+P(kf2,3) +
     &     P(kf1,4)+P(kf2,4)-P(KIA2,4)-P(KIA1,4))
         if(check.gt.1e-5) then 
            write(6,*) 'shat ',sqrt(sh_old),sqrt(sh_new),sqrt(shf),check
            write(6,*) ' meoffsh momentum checks ',P(kf1,1)+P(kf2,1),
     &      P(kf1,2)+P(kf2,2), 
     &      P(kf1,3)+P(kf2,3), 
     &      P(kf1,4)+P(kf2,4)-P(KIA2,4)-P(KIA1,4)
         endif

         PCM_FIN=SQRT(P(NF1,4)**2 - AM(1)**2)
         COS_FIN=P(NF1,3)/PCM_FIN

         call PYROBO(0,0,DTHETA,0.D0,0.d0,0.d0,0.d0)
         call PYROBO(0,0,0.D0,dphi,0.d0,0.d0,0.d0)

         CALL PYROBO(0,0,0.D0,0.D0,DSH(1)/DSH(4) ,DSH(2)/DSH(4), DSH(3)/
     +   DSH(4))
      Endif

      P(KIA1,3) = 0.5d0*(P(KIA1,3)-P(KIA1,4))
      P(KIA1,4) = DABS(P(KIA1,3))
      P(KIA2,3) = 0.5d0*(P(KIA2,3)+P(KIA2,4))
      P(KIA2,4) = DABS(P(KIA2,3))
       DO I=1,4
         DSHN(I)= P(KIA1,I) + P(KIA2,I)
      ENDDO
      SH_NEW = DOT(DSHN,DSHN)
      if(sh_new.lt.0) goto 20
      IF(P(NF1,5)-P(NF2,5).lt.0.01) then
         CALL PYROBO(0,0,0.D0,0.D0,-DSHN(1)/DSHN(4), -DSHN(2)/DSHN(4),
     +   -DSHN(3)/DSHN(4))

         P(KF1,4) = 0.5d0*DSQRT(SH_NEW)
         if(P(KF1,4)**2.le. AM(1)**2) goto 20
         PCM_NEW=SQRT(P(KF1,4)**2 - AM(1)**2)
         P(KF1,3) =PCM_NEW*COS_FIN
         IF(IMEPHI.EQ.0) THEN
            phi=2.d0*PI*dcasrn()
            phio = phi
         ELSEIF(IMEPHI.EQ.1) THEN
            phi=phio
         ENDIF
         P(KF1,1) = DSQRT(PCM_NEW**2-P(KF1,3)**2)*cos(phi)
         P(KF1,2) = DSQRT(PCM_NEW**2-P(KF1,3)**2)*sin(phi) 
         DO I=1,4
            P(KF2,I) = P(KIA1,I) + P(KIA2,I) - P(KF1,I)
         ENDDO
c	Call pylist(1)
         CALL PYROBO(0,0,0.D0,0.D0,DSHN(1)/DSHN(4) ,DSHN(2)/DSHN(4),
     +   DSHN(3)/DSHN(4))
      Endif
     
      K22 = DOT1(KIA2,KIA2)

	IF(Ires(1).ne.0.and.Ires(2).ne.0.and.IPRO.eq.10.and.
     &	(IRPBT.eq.1.or.IRPCT.eq.1).and.iccfm.ne.0) then

         do i=1,4
           pfin(i)=P(KF1,i)+P(KF2,i)
         end do
         plus=(pfin(4)-pfin(3))/2D0
         minus=(pfin(4)+pfin(3))/2D0
         P(KIA2,4)=minus
         P(KIA2,3)=minus
         P(KIA1,4)=plus
         P(KIA1,3)=-plus
c  select which leg should be on-shell
         IF(IMEPHI.EQ.0) THEN
            Rleg=dcasrn()
            Rlego = Rleg
         ELSEIF(IMEPHI.EQ.1) THEN
            Rleg=Rlego
         ENDIF
            
         if(Rleg.gt.0.5) then
           P(KIA2,2)=0D0
           P(KIA2,1)=0D0
           P(KIA1,2)=pfin(2)
           P(KIA1,1)=pfin(1)
         else
           P(KIA1,2)=0D0
           P(KIA1,1)=0D0
           P(KIA2,2)=pfin(2)
           P(KIA2,1)=pfin(1)
         endif
         
	else
c      write(6,*) ' should not enter ',IRES(1),IRPBT,IRPCT
         if(Iccfm.eq.0) then
            IF(IMEPHI.EQ.0) THEN
               Rleg=dcasrn()
               Rlego = Rleg
            ELSEIF(IMEPHI.EQ.1) THEN
               Rleg=Rlego
            ENDIF
         endif
	ENDIF

   10 continue

      K1K2 = DOT1(KIA1,KIA2)
      K1P2 = DOT1(KIA1,Kin2)
      K1P3 = DOT1(KIA1,KF2)
      K1P4 = DOT1(KIA1,KF1) 
      K2P1 = DOT1(KIA2,Kin1)
      K2P3 = DOT1(KIA2,KF2)
      K2P4 = DOT1(KIA2,KF1)
      K12 = DOT1(KIA1,KIA1)
      K22 = DOT1(KIA2,KIA2)
      P1P2 = DOT1(Kin1,Kin2)
      P1P3 = DOT1(Kin1,KF2)
      P1P4 = DOT1(Kin1,KF1)
      P2P3 = DOT1(Kin2,KF2)
      P2P4 = DOT1(Kin2,KF1)
      
      M2 = AM(1)**2

      TH = AM(1)**2 - 2.D0 * DOT1(KIA1,KF1) + K12
      UH = AM(1)**2 - 2.D0 * DOT1(KIA2,KF1) + K22
      SH = K12 + 2.D0 * DOT1(KIA1,KIA2) + K22
      
      test=Sh+th+uh-2*AM(1)**2-K12-K22	
	
      if(abs(test).gt.1e-4) then
         write(6,*) ' sh+th+uh ,m**2 ',Sh_new+th+uh,2*AM(1)**2+K12+
     +   K22,AM(1)
         write(6,*) ' sh+th+uh-m**2 ',Sh_new+th+uh-2*AM(1)**2-K12-
     +   K22
      endif
C Catani npb 366 (1991) 135
      IF(M2.gt.0.001) then
        E1 = -K12/M2
        E2 = -K22/M2
        else
        E1 = -K12*CUT
        E2 = -K22*CUT
      Endif
      SumAL=0d0
      SumAab = 0d0
      SumAnab = 0d0
      SumC1ab = 0d0
      SumC1nab = 0d0
      SumB1ab = 0d0
      SumB1nab = 0d0
      SumB2ab = 0d0
      SumB2nab = 0d0
      SumC2ab = 0d0
      SumC2nab = 0d0
      IF(E1.LT.CUT.AND.E2.LT.CUT) THEN
C all incoming partons are on mass shell
C parton 1 on shell, parton 2 on shell
C eq B23,B24
         SumAab = K1K2/K1P4/K1P3*((K1P3**2 + K1P4**2)/K1K2 + 2D0*M2*
     +   (1D0 - M2*K1K2/2D0/K1P4/K1P3))
         SumAnab = -2D0/K1K2*((K1P3**2 + K1P4**2)/K1K2 + 2D0*M2*(1D0 -
     +   M2*K1K2/2D0/K1P4/K1P3))
         p3p4 = DOT1(KF1,KF2)
c	for g*g -> gg	   
C labels as in paper JHEP 092009121 (note: p3 <-> p4 in code) 
	   SumB1nab=(k1k2/k1p2)**2*(p2p3**4+k1p2**4+p2p4**4)*
     &   (p3p4*k1p2+k1p3*p2p4+p2p3*k1p4)/(p2p3*k1p3*p3p4*k1p2*p2p4*k1p4)
c exchange momenta
	   SumB2nab=(k1k2/k2p1)**2*(p1p4**4+k2p1**4+p1p3**4)*
     &   (p3p4*k2p1+k2p4*p1p3+p1p4*k2p3)/(p1p4*k2p4*p3p4*k2p1*p1p3*k2p3)
     
c divide by 2, since later we add SumB1nab and SumB2nab
         SumB1nab=SumB1nab/2.
         SumB2nab=SumB2nab/2.
c         write(6,*) ' SumB1ab,SumB1nab ',SumB1ab,SumB1nab,SumB2ab,SumB2nab
         
c	for g*q -> gq
C labels as in paper JHEP 092009121 (note: p3 <-> p4 in code) 
 	   SumC1ab=(k1k2/k1p2)**2*(k1p2**2+p2p4**2)/(k1p3*p3p4)
	   SumC1nab= SumC1ab*((p3p4*k1p2+k1p3*p2p4)/(k1p4*p2p3)-1.)        
c exchange momenta
         SumC2ab=(k1k2/k2p1)**2*(k2p1**2+p1p3**2)/(k2p4*p3p4)
         SumC2nab= SumC2ab*((p3p4*k2p1+k2p4*p1p3)/(k2p3*p1p4)-1.)    
c         write(6,*) ' SumC1ab,SumC1nab ',SumC1ab,SumC1nab,SumC2ab,SumC2nab
         
c         write(6,*) ' SumC1ab,SumC1nab ',SumC1ab,SumC1nab
      ELSEIF(E1.LT.CUT.AND.E2.GT.CUT) THEN
C photoproduction, gluon off mass
C parton 1 on shell, parton 2 off shell
C eq. B21, B22
         SumAab = (K1K2/K1P2)**2*((P2P3**2+P2P4**2)/K1P4/K1P3 - 2D0*M2/
     +   K22*(P2P4/K1P3-P2P3/K1P4)**2)
         SumAnab = 0D0
         p3p4 = DOT1(KF1,KF2)
c	for g*g -> gg	   
C labels as in paper JHEP 092009121 (note: p3 <-> p4 in code)
	   SumB1nab=(k1k2/k1p2)**2*(p2p3**4+k1p2**4+p2p4**4)*
     &   (p3p4*k1p2+k1p3*p2p4+p2p3*k1p4)/(p2p3*k1p3*p3p4*k1p2*p2p4*k1p4)

c	for g*q -> gq
C labels as in paper JHEP 092009121 (note: p3 <-> p4 in code)
 	   SumC1ab=(k1k2/k1p2)**2*(k1p2**2+p2p4**2)/(k1p3*p3p4)
	   SumC1nab= SumC1ab*((p3p4*k1p2+k1p3*p2p4)/(k1p4*p2p3)-1.)        
c         write(6,*) ' E2 SumC1ab,SumC1nab ',SumC1ab,SumC1nab
         
      ELSEIF(E1.GT.CUT.AND.E2.LT.CUT) THEN
C DIS gluon onshell
C parton 1 off shell, parton 2 on shell
C eq. B21, B22
         SumAab = (K1K2/K2P1)**2*((P1P3**2+P1P4**2)/K2P4/K2P3 - 2D0*M2/
     +   K12*(P1P4/K2P3-P1P3/K2P4)**2)
         SumAnab = 0D0
c NIA1 = gluon
C NIA2 = quark
C NF1 = gluon
C NF2 = quark
         p3p4 = DOT1(KF1,KF2)
c	for g*g -> gg	   
C labels as in paper JHEP 092009121 (note: p3 <-> p4 in code)
	   SumB2nab=(k1k2/k2p1)**2*(p1p4**4+k2p1**4+p1p3**4)*
     &   (p3p4*k2p1+k2p4*p1p3+p1p4*k2p3)/(p1p4*k2p4*p3p4*k2p1*p1p3*k2p3)
     
c	for g*q -> gq (or gamma* q -> gq)
C labels as in paper JHEP 092009121 (note: p3 <-> p4 in code)
	   SumC2ab=(k1k2/k2p1)**2*(k2p1**2+p1p3**2)/(k2p4*p3p4)
	   SumC2nab= SumC2ab*((p3p4*k2p1+k2p4*p1p3)/(k2p3*p1p4)-1.)    
c         write(6,*) ' E1 SumC1ab,SumC1nab ',SumC1ab,SumC1nab
      ELSEIF(E1.GT.CUT.AND.E2.GT.CUT) THEN
C DIS, gluon off mass
C eq B13,B14
C parton 1 off shell, parton 2 off shell
         SumAab = 4D0*(K1P2*K2P1/P1P2**2)**2*(P1P2**2/(TH-M2)/(UH-M2) -
     +   1D0/K12/K22*(P1P2+2D0*P1P4*P2P3/(TH-M2)+ 2D0*P1P3*P2P4/(UH-
     +   M2))**2)
         Del = 2D0*P1P2*(2D0*P1P3*P2P4/P1P2 - 2D0*P1P4*P2P3/P1P2 - 
     +   K12*P1P4/K2P1 + K22*P2P4/K1P2 + K1P4-K2P4)
         SumAnab = (K1P2*K2P1/P1P2**2)**2*(4D0*P1P2**2*( -1D0/(TH-M2)/
     +   (UH-M2)-1D0/SH*(1D0/(TH-M2)-1D0/(UH-M2))* (P1P3/K2P1-P2P3/
     +   K1P2)+1D0/SH*P1P2/K2P1/K1P2)+2D0/K12/K22* (P1P2+4D0*P1P4
     +   *P2P3/(TH-M2)-Del/SH)* (P1P2+4D0*P1P3*P2P4/(UH-M2)+Del/SH))
         SumaL = SumAab*
     +    (1.d0 -1d0/3d0 - 8d0/3d0/SumAab*k2p1**2/p1p2**2*(
     +   (p2p4**2+p2p3**2)/(TH-M2)/(UH-M2) - (2d0*M2 + k12)/k22*
     +   (p2p3/(TH-M2) - p2p4/(UH-M2))**2))
      ENDIF


      IF(FIRST) THEN
	   XG2SF = 0.001
	   SFSCAL2 = 2.
	   STMAX2 = 5.
         call cauniglu(K(2,2),XG2SF,SFSCAL2,STMAX2,XPQ2)
         IQBAR = Iqqbar
         IF(IQBAR.EQ.1) THEN
            write(6,*) ' meoffsh - scale: qq-pair '
         ELSEIF(IQBAR.EQ.2) THEN
            write(6,*) ' meoffsh - scale: q or q_bar '
         ELSEIF(IQBAR.EQ.3) THEN
            write(6,*) ' meoffsh - scale: sh+qt2 '
         ENDIF
	      write(6,*) ' meoffsh - renorm-scale factor:',Scalfa
            write(6,*) ' meoffsh - fact-scale factor: ',scalfaf
         FIRST=.FALSE.
      ENDIF

      QG1 = P(NIA1,1)**2+ P(NIA1,2)**2
      QG2 = P(NIA2,1)**2+ P(NIA2,2)**2
      SFSCAL1 = DABS(QG1)
      SFSCAL2 = DABS(QG2)
c new scale
      PT2Q = (P(KF1,1)**2 +P(KF1,2)**2+P(KF2,1)**2 +P(KF2,2)**2)/2.
      xg_bar =  -qg2/(1.d0-xg)/SSS
      yph = dble(YY)
      yph_bar =  -qg1/(1.d0 - yph)/SSS


      If(Icoll.eq.1) then      
         xg_bar = 0.  ! hannes hack
         yph_bar = 0. ! hannes hack
      endif

	if(sfscal1.le.sfcut) sfscal1=sfcut 
	if(sfscal2.le.sfcut) sfscal2=sfcut 
      xsi_hard(1) = (xg+yph_bar)/(yph+xg_bar)	
      xsi_hard(2) = (yph+xg_bar)/(xg+yph_bar)

c which of the quarks is closer to the proton:
      rap1 = rapp-PYP(NF1,17)
      rap2 = rapp-PYP(NF2,17)
      zz1 = dabs(qg1)/yph/sss/xg
      zz2 = dabs(qg2)/yph/sss/xg
      zz1 = 1. - (P(NF1,4)+ P(NF1,3))/2./P(2,4)/yph
      zz2 = 1. - (P(NF2,4)+ P(NF2,3))/2./P(2,4)/xg

      ptstm1=dsqrt(dabs(qg1)+P(NF1,5)**2)
      ptstm2=dsqrt(dabs(qg2)+P(NF2,5)**2)

      If(IQBAR.eq.1) then
         stmax1 = yph*sqrt(xsi_hard(1)*SSS)
         stmax2 = xg*sqrt(xsi_hard(2)*SSS)
      elseif(IQBAR.EQ.2) then
         stmax1 = ptstm1/(1.-zz1)
         stmax2 = ptstm2/(1.-zz2)
      elseif(IQBAR.EQ.3) then
         QT2 = (P(KIA1,1)+P(KIA2,1))**2+(P(KIA1,2) +P(KIA2,2))**2
         stmax1 = sqrt(sh + QT2)
         stmax2 = stmax1
      else
         write(6,*) ' meoffsh. no valid scale selected ',IQBAR
      endif
	
      CALL PYROBO(0,N,0.D0,0.D0,-CM(1)/CM(4),-CM(2)/CM(4),-CM(3)/CM(4)
     +)

c new scale
c check on scale:
      SCALE = PT2Q+M2
      IF(IQ2.EQ.1) THEN
         SCALE = 4.D0*M2
      ELSEIF(IQ2.EQ.2) THEN
         SCALE = SH_NEW
      ELSEIF(IQ2.EQ.3) THEN
         SCALE = 4D0*M2 + PT2Q
      ELSEIF(IQ2.EQ.4) THEN
         SCALE = Q2
      ELSEIF(IQ2.EQ.5) THEN
         SCALE = Q2 + PT2Q + 4.D0*M2
      ELSEIF(IQ2.EQ.6) THEN
         SCALE = DABS(K12)+DABS(K22)
      ELSEIF(IQ2.EQ.7) THEN
         Scale = max(stmax1,stmax2)**2
      ELSE
         WRITE(6,*) ' NO VALID Q2 SCALE. STOP'
         STOP
      ENDIF
	
c multiplication with Qscal done in alphas.F      
      Dscale = scale
      alpha_s = ALPHAS(sqrt(dscale))
	QQg = Qg0

      ALPH_EM = ALPH

      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)
      
      do i=-6,6
         XPQ1(i) = 0.d0
         XPQ2(i) = 0.d0
      end do

      stmax1 = stmax1 *scalfaf
      stmax2 = stmax2 *scalfaf
      
      
      If(Ires(1).ne.0) then
         if(stmax1.lt.sqrt(sss)) then
            call cauniglu(K(1,2),XG1SF,SFSCAL1,STMAX1,XPQ1)
         endif
      endif
      If(Ires(2).ne.0) then
         if(stmax2.lt.sqrt(sss)) then
            call cauniglu(K(2,2),XG2SF,SFSCAL2,STMAX2,XPQ2)
         endif
      endif
     
      
      IF(Ires(1).eq.0.AND.Ires(2).ne.0) THEN
C gamma glu -> q qbar
         Cab = 16.D0*PI*PI*ALPHA_S*ALPH_EM*QF2
         CabB = 0D0
         epsilon = yph**2/2d0/(1d0-yph+yph**2/2d0)
         if ((SumAab - epsilon*SumAL).le.0) then
         endif
         SumAab = SumAab - epsilon*SumAL
         WTG= dble(XPQ2(0))
         SUMA = (Cab*SumAab + CabB*SumAnab)*WTG
cc this is for gamma*q -> gq
c	   XPQC2 = dble(XPQ2(3)+xpq2(-3))
c	   Cab = 16.D0*PI*PI*ALPHA_S*ALPH_EM*QF2*8./3.
c	   CabB = 0D0
c         SUMC2 = (Cab*SumC2ab + CabB*SumC2nab)*XPQC2
c         SUMC =  SUMC2
c	   Rleg = 0.6
         K(NIA1,2) = 22
         K(NIA2,2) = 21
      ELSEIF(Ires(1).ne.0.and.Ires(2).ne.0) THEN
C glu glu -> qqbar
         Cab = 16D0*PI**2*ALPHA_S**2/2D0/NC
         CabB = 16D0*PI**2*ALPHA_S**2/2D0*NC/(NC**2-1D0)
         Cab  = Cab*NQ
         CabB = CabB*NQ
         SUMA = 0
         If(IRPAT.ne.0) then
           WTG= dble(XPQ1(0)*XPQ2(0))
           SUMA = (Cab*SumAab + CabB*SumAnab)*WTG
         Endif
C glu glu -> glu glu
         SUMB = 0
         IF(IRPBT.ne.0) then
           Cab = 0  
           CabB = 16D0*PI**2*ALPHA_S**2*9D0/8D0
c include factor 1/2 for identical final state particles
           WTG= dble(XPQ1(0)*XPQ2(0))*0.5

           SUMB = (Cab*SumB1ab + CabB*SumB1nab)*WTG + (Cab*SumB2ab + CabB*SumB2nab)*WTG
         Endif
C glu q -> glu q 
         SUMC = 0
         SUMC1 = 0
         SUMC2 = 0
         IF(IRPCT.ne.0) then        
           Cab = 16D0*PI**2*ALPHA_S**2*4D0/9D0
           CabB = 16D0*PI**2*ALPHA_S**2/2D0 
           XPQC1=XPQ1(1)+XPQ1(2)+XPQ1(3)+XPQ1(-1)+XPQ1(-2)+XPQ1(-3)
           XPQC2=XPQ2(1)+XPQ2(2)+XPQ2(3)+XPQ2(-1)+XPQ2(-2)+XPQ2(-3)
           WTQ1= dble(XPQ2(0)*XPQC1)
           WTQ2= dble(XPQ1(0)*XPQC2)
           If(Rleg.gt.0.5) then
c KIA2 is onshell          
c factor  2 because of Rleg selection, but both contributions must be added
   	        SUMC2=(Cab*SumC2ab+CabB*SumC2nab)*WTQ2*2
c              write(6,*) ' Rleg > 0.5 ',Sumc1,SumC1ab,SumC1nab,XPQ1(0),XPQC1
              else
c KIA1 is onshell           
c factor  2 because of Rleg selection, but both contributions must be added
   	        SUMC2=(Cab*SumC1ab+CabB*SumC1nab)*WTQ1*2
c              write(6,*) ' Rleg < 0.5 ',Sumc2,SumC1ab,SumC1nab,XPQ2(0),XPQC2
           endif
           SUMC = SUMC1 + SUMC2
         Endif
      ELSEIF(Ires(1).eq.0.AND.Ires(2).eq.0) THEN
C gamma gamma -> q qbar
         Cab = 2D0*(4*PI*ALPH_EM)**2*QF2*QF2*NC
         CabB = 0
         SUMA = (Cab*SumAab + CabB*SumAnab)         
         K(NIA1,2) = 22
         K(NIA2,2) = 22
         K(NF1,2) = -KPA
         K(NF2,2) =  KPA
      ENDIF
      
      SUM2 = SUMA + SUMB + SUMC
      rntest = dcasrn()

	IF((SUMA+SUMB)/SUM2.LT.rntest) THEN
         IRESPRO = 3
      ELSEIF(SUMA/SUM2.LT.rntest) THEN
         IRESPRO = 2
      ELSE
         IRESPRO = 1
      ENDIF

      If(Irespro.eq.3) then
         if(Rleg.gt.0.5) then
            NFLT = -NFLAV-1
            rntest = dcasrn()
            QSUM = -DBLE(XPQC2)*rntest
  200       NFLT = NFLT + 1
            IF(NFLT.EQ.0) GOTO 200
            QSUM = QSUM + DBLE(XPQ2(NFLT))
            IF(QSUM.LE.0.0D0) GOTO 200
            IF(NFLT.GT.NFLAV) write(6,*) ' SUMC1 meoffsh NFL > NFLAV ',NFLT, NFLAV,SUMC1
            K(NIA1,2) = 21
            K(NIA2,2) = NFLT
            K(NF1,2) = 21
            K(NF2,2) =  NFLT
         Else
            NFLT = -NFLAV-1
            rntest = dcasrn()
            QSUM = -DBLE(XPQC1)*rntest
  201       NFLT = NFLT + 1
            IF(NFLT.EQ.0) GOTO 201
            QSUM = QSUM + DBLE(XPQ1(NFLT))
            IF(QSUM.LE.0.0D0) GOTO 201
            IF(NFLT.GT.NFLAV) write(6,*) ' SUMC2 meoffsh NFL > NFLAV ',NFLT, NFLAV,SUMC2
            K(NIA2,2) = 21
            K(NIA1,2) = NFLT
            K(NF2,2) = 21
            K(NF1,2) =  NFLT
         
         endif
c get color configuartion for gq --> gq
         A1 = 0.5
         A2 = 0.5
         ASUM = A1 + A2
         ICOLORA = 0
         IF(A1/ASUM.lt.dcasrn()) THEN
            ICOLORA=2
         ELSE
            ICOLORA=1
         ENDIF
      elseif(irespro.eq.2) then
         A1 = 1.
         A2 = 1.
         A3 = 1.
         ASUM = A1+A2+A3
         RNTEST = dcasrn()
         ICOLORA = 0
         IF((A1+A2)/ASUM.lt.rntest) THEN
            ICOLORA=3
         ELSEIF(A1/ASUM.lt.rntest) THEN
            ICOLORA=2
         ELSE
            ICOLORA=1
         ENDIF
         K(NIA2,2) = 21
         K(NIA1,2) = 21
         K(NF2,2) = 21
         K(NF1,2) = 21
      elseif(irespro.eq.1) then
            K(NIA2,2) = K(KIA2,2)
            K(NIA1,2) = K(KIA1,2)
            K(NF2,2) = K(KF2,2)
            K(NF1,2) = K(KF1,2)
      else
         Write(6,*) ' now valid irespro found: irespro = ',irespro
      endif

C new cut for on-shell DGLAP
	If(Iccfm.eq.0) then
	   If(ires(1).eq.1.and.(sfscal1.gt.scale*scalfa)) SUM2 = 0 
	   If(ires(1).eq.1.and.(sfscal1.gt.STMAX1**2)) SUM2 = 0 
	   If(ires(2).eq.1.and.(sfscal2.gt.scale*scalfa)) SUM2 = 0 
	   If(ires(2).eq.1.and.(sfscal2.gt.STMAX2**2)) SUM2 = 0 
	Endif

      IF(SUM2.LT.0D0) THEN
         if(sum2.ne.sum2) then
            write(6,*) ' meoffsh SUM2,xsi_hard ',SUM2,
     &        xsi_hard(1),xsi_hard(2)
            Write(6,*) ' meoffsh Cab,Aab,CabB,Anab ',Cab,SumAab,CabB,
     +      SumAnab
            write(6,*) ' meoffsh k12,k22 ',k12,k22
            write(6,*) ' meoffsh th,uh ',th,uh
         endif
         SUM2 = 0D0
      ENDIF
c apply pt cuts here
      PT2 = PYP(KF1,9)
ccc      PT2 = max(PYP(KF1,9),PYP(KF2,9)) ! hannes new
      CALL CUTG(PT2,NACC)
      IF(NACC.EQ.0) SUM2 = 0
      PT2 = PYP(KF2,9)
      CALL CUTG(PT2,NACC)
      IF(NACC.EQ.0) SUM2 = 0
c      PT2 = min(PYP(KF1,9),PYP(KF2,9)) ! hannes new
c      If(PT2.le.1.) SUM2 = 0 ! hannes new
c  cut on true rapidity      
c      if(Dabs(PYP(KF1,17)).gt.4.7) Sum2 = 0
c      if(Dabs(PYP(KF2,17)).gt.4.7) Sum2 = 0
cc      
c cut to ensure collinear treatment
      IF(Ires(1).ne.0..and.Ires(2).ne.0.and.IPRO.EQ.10.and.(IRPB.eq.1.or.IRPC.eq.1)) then
         If(Rleg.gt.0.5) then
            if(sfscal2.ge.scale) then
               sum2 = 0
            endif
	      If(sfscal2.gt.stmax2**2) SUM2 = 0 
         Else
            if(sfscal1.ge.scale) then
               sum2 = 0
            endif
	      If(sfscal1.gt.stmax1**2) SUM2 = 0 
         Endif
      endif
 
      WT1 = SUM2

      if(sum2.ne.sum2) then

         write(6,*) ' meoffsh xpq(0),alpha_s,alph_em',xpq1(0),xpq2(0),
     +      alpha_s,alph_em
         write(6,*) ' meoffsh SUM2,xsi_hard ',SUM2,
     &      xsi_hard(1),xsi_hard(2)
         Write(6,*) ' meoffsh Cab,Aab,CabB,Anab ',Cab,SumAab,CabB,SumAnab
         write(6,*) ' meoffsh k12,k22 ',k12,k22
         write(6,*) ' meoffsh th,uh ',th,uh
         write(6,*) ' P(KIA1) ',(P(KIA1,I),I=1,5)
         write(6,*) ' P(KIA2) ',(P(KIA2,I),I=1,5)
         write(6,*) ' P(KF1) ',(P(KF1,I),I=1,5)
         write(6,*) ' P(KF2) ',(P(KF2,I),I=1,5)

      endif
      CALL PYEDIT(13)
      RETURN
   20 WT1 = 0.d0
      CALL PYEDIT(13)
      RETURN
      END
      
