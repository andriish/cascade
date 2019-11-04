      SUBROUTINE mevq(WT1)
      IMPLICIT None
C
c
c use the momenta ala CCH
c qg* -> V q 
c S.Marzani and R. Ball. High Energy Resummation of Drell-Yan Processes. 
c Nucl. Phys., B814:246,264, 2009.
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
	
      Integer I,J,KIA1,KIA2,KF1,KF2,Kin1,Kin2
      Double Precision k_bold_sq, kq_bold, e_q, quark, mySH, myTH
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
	Double Precision rntest,Rleg,xsecnorm
      Integer NFLT
      Integer ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON /CASCOL/ ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      double Precision Aas
      Integer IRPAT,IRPBT,IRPCT,IRPDT,IRPET,IRPFT,IRPGT

      Double Precision dot,dot1,dcasrn
      Double Precision sL,rt,ru
      Double Precision Del,Nc
      Double Precision rap1,rap2,rapp,ptstm1,ptstm2,zz1,zz2
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
      Integer Nacc,itemp
      Double Precision gDY,GF,Vq2,Aq2
      Double Precision Vqq,Vqqsum(6),Vqqrn,S_in(6),SUMin
      Integer KW1,KW2,KPZ
      Integer kk,Inflav,NFLI

      
      Data Nc/3.d0/
      Logical ltest
c      data ltest/.true./
      data ltest/.false./
      Logical first
      data first/.true./
c      SAVE STMAX1,STMAX2, 
      Save XPQ2, phio   
c      write(6,*) ' begin mevq'
c      call pylist(1)
      IQBAR = Iqqbar
c      IQbar = 2
      Xsecnorm = 1.
      IF(INTER.EQ.2) THEN
         write(6,*) ' INTER = 2 not implemented for meoffsh'
         STOP
      ENDIF
      
c      write(6,*) IPRO,IRPA,IRPB,IRPC
      NPFIN = N
c	call pylist(1)
      SH = AM(1)**2 + AM(2)**2 + 2.D0 * DOT1(NF1,NF2)
      M2 = AM(1)**2
      WT1 = 0.d0
c      write(6,*) ' meoffsh ',am(1)
c      call pylist(1)
c	write(6,*) ' meoffsh: ',IRPAT,IRPBT,IRPCT,IRPDT,IRPET,IRPFT,IRPGT
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
c      call pylist(1)
         P(KIA1,3) = dsqrt(((sh_old-K12-K22)**2 - 4d0*K12*K22)/4d0/sh_old)
         P(KIA1,4) = 0.5d0*(sh_old+(-k22+k12))/sqrt(sh_old)
         P(KIA2,3) =-dsqrt(((sh_old-K12-K22)**2 - 4d0*K12*K22)/4d0/sh_old)
         P(KIA2,4) = 0.5d0*(sh_old-(-k22+k12))/sqrt(sh_old) 

c      write(6,*) ' Kia1,2 = ',KIA1,KIA2,' KF1,2 = ',KF1,KF2
c      call pylist(1)
            
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
            write(6,*) ' mevq: momentum checks ',P(kf1,1)+P(kf2,1),
     &      P(kf1,2)+P(kf2,2), 
     &      P(kf1,3)+P(kf2,3), 
     &      P(kf1,4)+P(kf2,4)-P(KIA2,4)-P(KIA1,4)
            call pylist(1)
            stop
            
         endif
ctest
         PCM_FIN=SQRT(P(NF1,4)**2 - AM(1)**2)
         COS_FIN=P(NF1,3)/PCM_FIN
cend     
         call PYROBO(0,0,DTHETA,0.D0,0.d0,0.d0,0.d0)
         call PYROBO(0,0,0.D0,dphi,0.d0,0.d0,0.d0)

         CALL PYROBO(0,0,0.D0,0.D0,DSH(1)/DSH(4) ,DSH(2)/DSH(4), DSH(3)/
     +   DSH(4))
      Endif
c      call pylist(1)

ctest
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
c      	write (6,*) ' masses ',P(NF1,5),P(NF2,5),P(NF1,5)-P(NF2,5)
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
     
     

ctesting
      K22 = DOT1(KIA2,KIA2)
     
	IF(Ires(1).ne.0.and.Ires(2).ne.0.and.(IPRO.eq.501.or.IPRO.eq.511).
     &	and.iccfm.ne.0) then
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
         Rleg=dcasrn()
c         Rleg = 0.
         xsecnorm = 2.
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
           itemp = KIA1
c           KIA1 = KIA2
c           KIA2 = itemp
         endif
         
	else
c      write(6,*) ' should not enter ',IRES(1),IRPBT,IRPCT
	ENDIF
Cend testing
c NIA1 = gluon
C NIA2 = quark
C NF1 = vectorboson
C NF2 = quark
      KPZ = K(NF1,2)

   10 continue

      K12 = DOT1(KIA1,KIA1)
      K22 = DOT1(KIA2,KIA2)
      
c      k12 = - (P(KIA1,1)**2+ P(KIA1,2)**2)
c      k22 = - (P(KIA2,1)**2+ P(KIA2,2)**2)
c      k12 = - (P(NIA1,1)**2+ P(NIA1,2)**2)
c      k22 = - (P(NIA2,1)**2+ P(NIA2,2)**2)
c      write(6,*) ' nia1 ',(P(Kia1,I),i=1,4),(P(Nia1,I),i=1,4)
c      write(6,*) ' nia2 ',(P(Kia2,I),i=1,4),(P(Nia2,I),i=1,4)
      
      SH = K12 + 2.D0 * DOT1(KIA1,KIA2) + K22
c define t/u according to paper
      TH = DOT1(KF1,KF1) - 2.D0 * DOT1(KIA2,KF1) + DOT1(KIA2,KIA2)
      UH = DOT1(KF1,KF1) - 2.D0 * DOT1(KIA1,KF1) + DOT1(KIA1,KIA1)
     


cc      test=Sh_new+th+uh-2*AM(1)**2-K12-K22
c      test=Sh+th+uh-2*AM(1)**2-K12-K22
	
	
      if(abs(test).gt.1e-4) then
         write(6,*) ' sh+th+uh ,m**2 ',Sh_new+th+uh,2*AM(1)**2+K12+
     +   K22,AM(1)
         write(6,*) ' sh+th+uh-m**2 ',Sh_new+th+uh-2*AM(1)**2-K12-
     +   K22
      endif
C Catani npb 366 (1991) 135
c      IF(M2.gt.0.001) then
c        E1 = -K12/M2
c        E2 = -K22/M2
c        write(6,*) ' mass ',m2
c        else
        E1 = -K12*CUT
        E2 = -K22*CUT
c        write(6,*) ' checks k12,k22 ',k12,k22, rleg
c        write(6,*) ' no mass ',m2
c      Endif

c      write(6,*) ' K12 ',K12,' K22 ',K22
c      IF(E2.GT.CUT) write(6,*) ' meoffsh ',E1,E2,CUT,M2
      IF(E1.LT.CUT.AND.E2.LT.CUT) THEN
C all incoming partons are on mass shell
C parton 1 on shell, parton 2 on shell
cc	   write(6,*) ' in on-shell case ',k12,k22,e1,e2
         SumC = (-1.)/3.*(TH/SH + SH/TH + 2*M2*UH/(SH*TH) 
     &                  + UH/SH + SH/UH + 2*M2*TH/(SH*UH)  )
         
      ELSEIF(E1.LT.CUT.AND.E2.GT.CUT) THEN
C parton 1 on shell, parton 2 off shell
c NIA1 = quark
C NIA2 = gluon
C NF1 = vectorboson
C NF2 = quark
C insert here the ME for DY: SumC1ab
         k_bold_sq = P(KIA2,1)*P(KIA2,1) + P(KIA2,2)*P(KIA2,2) 
         kq_bold = P(KIA2,1)*P(KF1,1) + P(KIA2,2)*P(KF1,2)
         

         SumC = (-1.)/3.*(UH/SH + SH/UH + M2*k_bold_sq*(1./(SH*SH)+1./(UH*UH)) +
     + 4.*kq_bold/SH - (1.-kq_bold/k_bold_sq)*4.*M2*kq_bold/(UH*UH) + 
     + 2./(SH*UH)*((M2+k_bold_sq)*(k_bold_sq-2.*kq_bold)+2.*kq_bold*kq_bold+k_bold_sq*(SH-UH)))
         if (SumC.LT.0) then 
c           SumC = 0.
           write(6,*) ' mevq: ERROR ',SumC,dot1(kia1,kia1),dot1(kia2,kia2)
         endif
         
      ELSEIF(E1.GT.CUT.AND.E2.LT.CUT) THEN
C parton 1 off shell, parton 2 on shell
         SumC = 0
c NIA1 = gluon
C NIA2 = quark
C NF1 = vectorboson
C NF2 = quark
C insert here the ME for DY: SumC1ab
         k_bold_sq = P(KIA1,1)*P(KIA1,1) + P(KIA1,2)*P(KIA1,2) 
         kq_bold = P(KIA1,1)*P(KF1,1) + P(KIA1,2)*P(KF1,2)
         

         SumC = (-1.)/3.*(TH/SH + SH/TH + M2*k_bold_sq*(1./(SH*SH)+1./(TH*TH)) +
     + 4.*kq_bold/SH - (1.-kq_bold/k_bold_sq)*4.*M2*kq_bold/(TH*TH) + 
     + 2./(SH*TH)*((M2+k_bold_sq)*(k_bold_sq-2.*kq_bold)+2.*kq_bold*kq_bold+k_bold_sq*(SH-TH)))
c         SumC = (-1.)*e_q*e_q*(4.*3.14159)/3.*(TH/SH + SH/TH + 2*M2*UH/(SH*TH))
c         SumC = (-1.)*e_q*e_q*(4.*3.14159)/3.*(TH/SH + SH/TH)

         if (SumC.LT.0) then 
c           SumC = 0.
           write(6,*) ' mevq ERROR ',SumC,dot1(kia1,kia1),dot1(kia2,kia2)
         endif
c        write(6,*) ' 2nd sumc = ',sumc
c       write(6,*) kia1,kia2,kf1,kf2
c       call pylist(1)

      ELSEIF(E1.GT.CUT.AND.E2.GT.CUT) THEN
C parton 1 off shell, parton 2 off shell
         SumC = 0
      ENDIF


      IF(FIRST) THEN
	   XG2SF = 0.001
	   SFSCAL2 = 2.
	   STMAX2 = 5.
         call cauniglu(K(2,2),XG2SF,SFSCAL2,STMAX2,XPQ2)
         IQBAR = Iqqbar
         IF(IQBAR.EQ.1) THEN
            write(6,*) ' meoffsh - scale: qq-pair '
         ELSE
            write(6,*) ' meoffsh - scale: q or q_bar '
         ENDIF
	      write(6,*) ' meoffsh - scale for alphas ',Qscal
         FIRST=.FALSE.
      ENDIF

c	write(6,*) K12,k22,K(nia1,2),nia1,K(nia1,2),nia2
c	call pylist(1)
      QG1 = P(NIA1,1)**2+ P(NIA1,2)**2
      QG2 = P(NIA2,1)**2+ P(NIA2,2)**2
c      QG1 = P(KIA1,1)**2+ P(KIA1,2)**2
c      QG2 = P(KIA2,1)**2+ P(KIA2,2)**2
      SFSCAL1 = DABS(QG1)
      SFSCAL2 = DABS(QG2)
c      write(6,*) 'k12 ',k12,' k22 ',k22
c      write(6,*) 'qg1 ',qg1,' qg2 ',qg2
c new scale
      PT2Q = (P(KF1,1)**2 +P(KF1,2)**2+P(KF2,1)**2 +P(KF2,2)**2)/2.
c      PT2Q = (P(KF1,1)**2 +P(KF1,2)**2)
      xg_bar =  -qg2/(1.d0-xg)/SSS
      yph = dble(YY)
      yph_bar =  -qg1/(1.d0 - yph)/SSS

	if(sfscal1.le.sfcut) sfscal1=sfcut 
	if(sfscal2.le.sfcut) sfscal2=sfcut 
      xsi_hard(1) = (xg+yph_bar)/(yph+xg_bar)	
      xsi_hard(2) = (yph+xg_bar)/(xg+yph_bar)

c which of the quarks is closer to the proton:
      rap1 = rapp-PYP(NF1,17)
      rap2 = rapp-PYP(NF2,17)
      zz1 = dabs(qg1)/yph/sss/xg
      zz2 = dabs(qg2)/yph/sss/xg
c      write(6,*) ' old zz1,2 ',zz1,zz2
      zz1 = 1. - (P(NF1,4)+ P(NF1,3))/2./P(2,4)/yph
      zz2 = 1. - (P(NF2,4)+ P(NF2,3))/2./P(2,4)/xg
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
         write(6,*) ' meoffsh. no valid scale selected ',IQBAR
      endif
	
      CALL PYROBO(0,N,0.D0,0.D0,-CM(1)/CM(4),-CM(2)/CM(4),-CM(3)/CM(4)
     +)

c new scale
c check on scale:
      SCALE = PT2Q+M2
      IF(IQ2.EQ.1) THEN
         SCALE = M2
      ELSEIF(IQ2.EQ.2) THEN
         SCALE = SH_NEW
      ELSEIF(IQ2.EQ.3) THEN
        SCALE = M2 + PT2Q
      ELSEIF(IQ2.EQ.6) THEN
         SCALE = DABS(K12)+DABS(K22)
      ELSEIF(IQ2.EQ.7) THEN
         Scale = max(stmax1,stmax2)**2
      ELSE
         WRITE(6,*) ' NO VALID Q2 SCALE. STOP'
         STOP
      ENDIF
	
c	write(6,*) ' scales stmax = ',stmax,' sqrt(scale) = ',sqrt(scale)
      Dscale = scale
      alpha_s = ALPHAS(sqrt(dscale))
	QQg = Qg0

c         write(6,*) ' new alphas ',alpha_s

      ALPH_EM = ALPH
c	alph_em = 1./137.
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(M2)
      
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
         SUMA = 0
      ELSEIF(Ires(1).ne.0.and.Ires(2).ne.0) THEN
         SUMA = 0
c calculate it for gamma production 
c        16\pi^2*alpha_s*\alpha_em *ME^2         
c calculate it for W production (from Ellis Stirling Webber, eq 9.53 +corr)
c        4\pi*alpha_s*\sqrt{2}*GF*Mw^2*|Vqq|^2 *ME^2         
   	   GF = 1.166E-5
         if(IABS(KPZ).eq.23.and.IPRO.eq.501) then
C glu q ->  Z q 
            if(Rleg.gt.0.5) then
               XPQC1=XPQ2(1)+XPQ2(2)+XPQ2(3)+XPQ2(-1)+XPQ2(-2)+XPQ2(-3)
               WTQ1= dble(XPQ1(0)*XPQC1)
            else
               XPQC2=XPQ1(1)+XPQ1(2)+XPQ1(3)+XPQ1(-1)+XPQ1(-2)+XPQ1(-3)
               WTQ1= dble(XPQ2(0)*XPQC2)
            endif
c calculate it for Z production (from Ellis Stirling Webber, eq 9.53 +corr)
c        4\pi*alpha_s*\sqrt{2}*GF*Mz^2*(Vq^2 + Aq^2) *ME^2 
c        gDY =  coupling to vector boson
c        dDY for W : \pi * g_weak *|Vqq|^2 
c        dDY for Z0: \pi * g_weak *(Vq^2 + Aq^2) 
c        g_weak = sqrt(2)*GF *Mz^2*(Vq^2 + Aq^2)
c  vector and axial couplings for Z         
c                (u,c,t, 2/3)                  (d,s,b, -1/3)         
c         Vq2 = 2*(0.5 - 2.*2./3.*sin2w)**2 + 3*(-0.5 + 2.*1./3.*sin2w)**2
c         Aq2 = 2*0.5*0.5 + 3*0.5*0.5
c    only valence quarks
             Vq2 = (0.5 - 2.*2./3.*sin2w)**2 + (-0.5 + 2.*1./3.*sin2w)**2
             Aq2 = 0.5*0.5 + 0.5*0.5
c         write(6,*) Vq2,Aq2
             gDY = 4.*pi *sqrt(2.)* M2 * GF * (Vq2 + Aq2)     
   	       SUMA=SumC*alpha_s*gDY*WTQ1
         elseif(IABS(KPZ).eq.24.and.IPRO.eq.511) then
C glu q ->  W q 
            SUMA = 0.
            
c parton in p only up to charm (not clear what remnant is when a b or t is removed)  
c sum over incoming flavors 
            Nflav = 2        
            do Inflav=1,Nflav
               S_in(Inflav) = 0.
               vqqsum(inflav) = 0.
               do J=1,3 
                  KW1 = Inflav
                  KW2 = 2*J-1+MOD(Inflav,2)
                  if(iabs(KW2).le.6) Then
                  if(KPZ.gt.0) then
                    if(KW1.eq.2.or.KW1.eq.4.or.KW1.eq.6) then
                      KW2 = KW2
                    endif
                    if(KW1.eq.1.or.KW1.eq.3.or.KW1.eq.5) then
                      KW1 = -KW1
                      KW2 = -KW2
                    endif
                  elseif(KPZ.lt.0) then
                    if(KW1.eq.2.or.KW1.eq.4.or.KW1.eq.6) then
                      KW1 = -KW1
                      KW2 = -KW2
                    endif
                    if(KW1.eq.1.or.KW1.eq.3.or.KW1.eq.5) then
                      KW2 = KW2
                    endif
                  else
                  
                  endif
                  
                  
c                  write(6,*) 'flavor ',KW1,KW2,KPZ
                  
                  if(Rleg.gt.0.5) then
                     WTQ1= dble(XPQ1(0)*XPQ2(KW1))
                     else
                     WTQ1= dble(XPQ2(0)*XPQ1(KW1))
                  endif
                  
                  Vqq=VCKM((Iabs(KW1)+1)/2,(Iabs(KW2)+1)/2)
                  vqqsum(inflav) = vqq+vqqsum(inflav)
c                  write(6,*) ' vqq ',vqq,vqqsum(inflav),inflav,kw1,kw2
                  gDY = 4.*pi *sqrt(2.)* M2 * GF * Vqq
                  SUMin = SumC*alpha_s*gDY*WTQ1
                  S_in(Inflav) = S_in(Inflav)+SUMin
c                  write(6,*) ' Suma ',Suma,S_in(inflav),inflav
   	            SUMA=SUMA+SUMin
c                  write(6,*) ' xpq1 = ',(xpq1(kk),kk=0,3)
c                  write(6,*) ' xpq2 = ',(xpq2(kk),kk=0,3)
                  
c                  write(6,*) ' suma ',suma,vqq,KW1    
    
                  Endif
               End do 
            End do

         endif
c 	 write(6,*) 'SUMC ',e_q,alpha_s,alph_em

      ELSEIF(Ires(1).eq.0.AND.Ires(2).eq.0) THEN
C gamma gamma -> 
      ENDIF
      SUM2 = SUMA*xsecnorm
c      if(vqqsum.le.0.9)       write(6,*) 'vqqsum = ',vqqsum,K(NIA1,2),K(NIA2,2)
c      write(6,*) ' mevq ',xsecnorm,rleg
      if(sum2.le.0) goto 20
      IF(KPZ.eq.23) Then
c select flavor of incoming partons      
      rntest = dcasrn()     
      if(Rleg.gt.0.5) then
         NFLT = -NFLAV-1
         rntest = dcasrn()
         QSUM = -DBLE(XPQC1)*rntest
  200    NFLT = NFLT + 1
         IF(NFLT.EQ.0) GOTO 200
c               write(6,*) ' test ',QSUM,XPQC1,NFLT,XPQ2(NFLT)
         QSUM = QSUM + DBLE(XPQ2(NFLT))
         IF(QSUM.LE.0.0D0) GOTO 200
         IF(NFLT.GT.NFLAV) write(6,*) ' meoffsh NFL > NFLAV ',
     +      NFLT, NFLAV
         K(NIA1,2) = 21
         K(NIA2,2) = NFLT
         K(NF2,2) =  NFLT
      Else
         NFLT = -NFLAV-1
         rntest = dcasrn()
         QSUM = -DBLE(XPQC2)*rntest
  201    NFLT = NFLT + 1
         IF(NFLT.EQ.0) GOTO 201
cc               write(6,*) ' test ',QSUM,XPQC1,NFLT,XPQ2(NFLT)
         QSUM = QSUM + DBLE(XPQ1(NFLT))
         IF(QSUM.LE.0.0D0) GOTO 201
         IF(NFLT.GT.NFLAV) write(6,*) ' mevq NFL > NFLAV ',
     +      NFLT, NFLAV
         K(NIA2,2) = 21
         K(NIA1,2) = NFLT
         K(NF2,2) =  NFLT         
      endif

      elseif(IABS(KPZ).eq.24) then
         NFLI = 0
         rntest = dcasrn()
         Qsum = -SUMA*rntest
  301    NFLI = NFLI+1
         IF(NFLI.eq.0) goto 301
c         write(6,*) ' qsum , S_in, NFLT ',qsum,s_in(NFLT),NFLT
         Qsum = Qsum + S_in(NFLI)
         IF(NFLI.GT.NFLAV) write(6,*) ' mevq NFL > NFLAV ',NFLI, NFLAV
         IF(QSUM.LE.0.0D0) GOTO 301
c select outgoing flavor
ccc         do Inflav=1,Nflav
         Inflav = Iabs(NFLI)
         NFLT = 0            
         rntest = dcasrn()
         Qsum = -Vqqsum(inflav)*rntest
 302     NFLT = NFLT+1
         KW1 = Inflav
         KW2 = 2*NFLT-1+MOD(Inflav,2)
         if(iabs(KW2).le.6) Then
                  if(KPZ.gt.0) then
                    if(KW1.eq.2.or.KW1.eq.4.or.KW1.eq.6) then
                      KW2 = KW2
                    endif
                    if(KW1.eq.1.or.KW1.eq.3.or.KW1.eq.5) then
                      KW1 = -KW1
                      KW2 = -KW2
                    endif
                   elseif(KPZ.lt.0) then
                    if(KW1.eq.2.or.KW1.eq.4.or.KW1.eq.6) then
                      KW1 = -KW1
                      KW2 = -KW2
                    endif
                    if(KW1.eq.1.or.KW1.eq.3.or.KW1.eq.5) then
                      KW2 = KW2
                    endif
                 endif
         endif
               Vqq=VCKM((Iabs(KW1)+1)/2,(Iabs(KW2)+1)/2)
c               write(6,*) ' KW1,KW2 ',KW1,KW2
               IF(NFLT.eq.0) goto 302
c               write(6,*) ' Vqqsum, qsum , Vqq, NFLT ',Vqqsum(inflav),qsum,Vqq,NFLT,inflav,J
               Qsum = Qsum + Vqq
c               IF(NFLT.GT.NFLAV) write(6,*) ' mevq final NFL > NFLAV ',NFLT, NFLAV
c               write(6,*) ' qsum',qsum
               IF(QSUM.LE.0.0D0) GOTO 302
               
cc         End do
         if(Rleg.gt.0.5) then
            K(NIA1,2) = 21
            K(NIA2,2) = NFLI
            K(NF2,2) =  NFLT
            else
            K(NIA2,2) = 21
            K(NIA1,2) = NFLI
            K(NF2,2) =  NFLT        
          endif
      endif
c      write(6,*) ' flavors  ',K(NIA1,2),K(NIA2,2),K(NF2,2)
C new cut for on-shell DGLAP
	If(Iccfm.eq.0) then
	   If(ires(1).eq.1.and.(sfscal1.gt.scale*scalfa)) SUM2 = 0 
	   If(ires(1).eq.1.and.(sfscal1.gt.STMAX1**2)) SUM2 = 0 
	   If(ires(2).eq.1.and.(sfscal2.gt.scale*scalfa)) SUM2 = 0 
	   If(ires(2).eq.1.and.(sfscal2.gt.STMAX2**2)) SUM2 = 0 
	Endif

c      Write(6,*) ' meoffsh Cab,Aab,CabB,Anab ',Cab,Aab,CabB,Anab
c      Write(6,*) ' meoffsh Wtg ',Wtg,XG1SF,XG2SF
c      Write(6,*) ' meoffsh Cab,Aab,CabB,Anab ',Cab,Aab,CabB,Anab
      IF(SUM2.LT.0D0) THEN
c     Write(6,*) ' meoffsh Cab,Aab,CabB,Anab ',Cab,Aab,CabB,Anab
c     write(6,*) ' SUM2 ',SUM2
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
c      PT2 = PYP(KF1,9)
c      CALL CUTG(PT2,NACC)
c      write(6,*) K(KF1,2),K(KF2,2)
c      IF(NACC.EQ.0) SUM2 = 0
      PT2 = PYP(KF2,9)
      CALL CUTG(PT2,NACC)
      IF(NACC.EQ.0) SUM2 = 0
c      write(6,*) ' ptcut ',pt2cut(ipro)
c      call pylist(1)
c cut to ensure collinear treatment
      IF(Ires(1).ne.0) then
         If(Rleg.gt.0.5) then
            if(sfscal2.ge.scale) then
               sum2 = 0
            endif
         Else
            if(sfscal1.ge.scale) then
               sum2 = 0
            endif
         Endif
      endif
c apply additonal factor of 2 in case we select one on shell leg
c since we do not add the xsections but rather put it symmetric only
 
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
c      write(6,*) ' end of mevq '
c      call pylist(1)
      RETURN
   20 WT1 = 0.d0
      CALL PYEDIT(13)
       RETURN
      END
