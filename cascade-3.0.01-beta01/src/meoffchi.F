      SUBROUTINE meoffchi(WT1)
      IMPLICIT None
C
C   Gluon* Gluon* ----> chi
C
C         P1-----*
C                * *
C                *   *-------- chi
C                * *
C         P2-----*
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
#include "caspsipol.inc"
      Integer IGENFL,IMEPHI
      COMMON/GENWEI/IGENFL,IMEPHI
      Double Precision XPQ1(-6:6),XPQ2(-6:6),SFSCAL1,SFSCAL2
	Double Precision STMAX1,STMAX2
	common/myscale/STMAX1,STMAX2
      DOUBLE PRECISION DSH(4),DSHN(4)
      Double Precision DScale
      Double Precision SH_NEW,M2,xh,sig,PT2Q
	Double Precision K1K2,K12,K22
      Double Precision  scale
      Double Precision Sfcut
      DATA Sfcut /1.d-3/ 
c      DATA Sfcut /1.d-9/ 
	Double Precision WT1,SH,qg2,qg1,wtg,sum2
      Double Precision alph_em,pcm_new,alphas,alpha_s,alpha_em
      Double Precision yph,yph_bar,xsi_hard(2),gf
      Double Precision dot,dot1,dcasrn
      Double Precision rap1,rap2,rapp,ptstm,zz1,zz2
      Integer npfin
      Integer IQBAR
	Integer KIA1,KIA2,KF1,KF2,Kin1,Kin2
	Integer I
      Integer ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Double Precision E1,E2
	Double Precision CUT
      DATA CUT /1.d-6/ ! as smallx
c      DATA CUT /1.d-2/ ! 
	
      DOUBLE PRECISION SIN2W,XMW2
      COMMON/ELWEAK/SIN2W,XMW2
	Double Precision QQg,lambda
	Double Precision Qscal
	Integer Iqqbar
	Common/GLUDAT2/Qscal,Iqqbar



c commons for S.baranov
      Integer ITYPE,MODEL,IB1,IB2,NJ,NL
      COMMON/TYPE2/ITYPE,MODEL,IB1,IB2,NJ,NL
      double Precision ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/CONST2/ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      Double Precision g1,g2,pj 
      COMMON/MOMEN2/g1(4),g2(4),pj(4) 
      Double Precision XChi0,XChi1,XChi2
      Integer kchitest
            
      Double Precision rntest
      CHARACTER * 12 CNAM

      Data IQBAR/1/
c      Data IQBAR/2/

      Logical ltest
c      data ltest/.true./
      data ltest/.false./
      Logical first
      data first/.true./
c     call pylist(1)
      NPFIN = N
      SH = DOT1(NF1,NF1) + 2.D0 * DOT1(NF1,NF2) * DOT1(NF2,NF2)
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
      if(ltest) then
         SH_NEW = DOT(DSH,DSH)
         goto 10
      endif
      P(KIA1,3) = 0.5d0*(P(KIA1,3)-P(KIA1,4))
      P(KIA1,4) = DABS(P(KIA1,3))
      P(KIA2,3) = 0.5d0*(P(KIA2,3)+P(KIA2,4))
      P(KIA2,4) = DABS(P(KIA2,3))

      IF(ICCFM.EQ.0) Then
	   If(ires(1).eq.1) then
            P(KIA1,1) = 0D0
            P(KIA1,2) = 0D0
         endif
	   If(ires(2).eq.1) then
            P(KIA2,1) = 0D0
            P(KIA2,2) = 0D0
         endif
      ENDIF

      DO I=1,4
         DSHN(I)= P(KIA1,I) + P(KIA2,I)
      ENDDO
      SH_NEW = DOT(DSHN,DSHN)
c      IF(SH_NEW.LE.AM(1)**2) GOTO 20


   10 continue
      K12 = DOT1(KIA1,KIA1)
      K22 = DOT1(KIA2,KIA2)

      

c here is the calculation of baranov et al. 

c         write(6,*) ' checking frame',M2,KIA1,KIA2,KF1,KF2
c         call pylist(1)
      PT2Q = P(KF1,1)**2 +P(KF1,2)**2
      SCALE = PT2Q+M2
      IF(IQ2.EQ.1) THEN
         SCALE = M2
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
c multiplication with Qscal done in alphas.F      
      Dscale = scale
      alpha_s = ALPHAS(sqrt(dscale))
      M2 = AM(1)**2
ccc      M2 = 3.4*3.4
      ALS = alpha_s
      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)
      ALP = ALPH_EM
      XJ = sqrt(M2)
      XC = XJ/2
      XC2 = XC**2
      XJ2 = XJ**2
C      PSI = 0.04/XJ  ! for eta_c (without strong corrections)
C      PSI = 0.40/XJ  ! for eta_b (without strong corrections)
      if(IHFLA.eq.4) then
c         PSI = 0.006/XJ ! for chi_c (from some potential model)
	    PSI = 0.075/4./pi/xj ! Eichten/Quigg  1P state
          if(i23s.eq.2) PSI = 0.102/4./pi/xj ! Eichten/Quigg  2P state
      elseif(IHFLA.eq.5) then
c         PSI = 0.12/XJ  ! for chi_b (from some potential model)
	    PSI = 1.417/4./pi/xj ! Eichten/Quigg  1P state
          if(i23s.eq.2) PSI = 1.653/4./pi/xj ! Eichten/Quigg  2P state
          if(i23s.eq.3) PSI = 1.79/4./pi/xj  ! Eichten/Quigg  3P state
cc         write(6,*) ' mechi ',xj,xc
      else
         write(6,*) ' meoffchi: wrong flavor selected, program stopped '
         stop
      endif
ccc      ALS = 0.25
      if(ICCFM.EQ.1) then 
           MODEL=1
           else
           Model = 0
      endif
      IF(ires(1).eq.0.and.Ires(2).eq.0) then
c            itype = 30 ! 30 = lepton*lepton [gam+gam->Eta(singlet)]
            itype = 32 ! 32 = lepton*lepton [gam+gam->Chi(singlet)]
            ib1 = 1
            ib2 = 1
      elseif(ires(1).eq.0.and.Ires(2).eq.1) then
c            itype = 50 ! 50 = lepton*proton [gam+glu->Eta(octet)]
            itype = 52 ! 52 = lepton*proton [gam+glu->Chi(octet)]
            ib1 = 1
            ib2 = 8
      elseif(ires(1).eq.1.and.Ires(2).eq.1) then
c      	itype = 10 ! 10 = proton*proton [glu+glu->Eta(singlet)]
      	itype = 12 ! 12 = proton*proton [glu+glu->Chi(singlet)]
c            itype = 40 ! 40 = proton*proton [glu+glu->Eta(octet)]
c            itype = 41 ! 41 = proton*proton [glu+glu->Psi(octet)]
c            itype = 42 ! 42 = proton*proton [glu+glu->Chi(octet)]
            ib1 = 8
            ib2 = 8
      endif
      NJ  = 1           ! 0 for eta_c; 1 for chi_c, J/psi(octet)
      NL  = 1           ! 0 for eta_c; 1 for chi_c

c include a very small kt for the calc of the polarisation
c othersise errors produced         
      g1(1) = P(KIA1,3)
      g1(2) = P(KIA1,1)+1.e-10*P(KIA1,4)
      g1(3) = P(KIA1,2)+1.e-10*P(KIA1,4)
      g1(4) = P(KIA1,4)
      g2(1) = P(KIA2,3)
      g2(2) = P(KIA2,1)+1.e-10*P(KIA2,4)
      g2(3) = P(KIA2,2)+1.e-10*P(KIA2,4)
      g2(4) = P(KIA2,4)

      pj(1) = P(KF1,3)
      pj(2) = P(KF1,1)
      pj(3) = P(KF1,2)
      pj(4) = P(KF1,4)
      
c         v1 = 
c         v2 = 
ccccc this is for serguei calculation
      CALL METRIC2
      Call xsec2(XChi0,XChi1,XChi2)
ccccc this is using Kniehls ME
cccc	Call xsec3(Xchi0,Xchi1,Xchi2)      
c      write(6,*) ' XChi0,XChi1,XChi2 ', XChi0,XChi1,XChi2
      sum2 = Xchi0+Xchi1+Xchi2
      if(sum2.ne.sum2) then
           sum2=-9999.
           write(6,*) ' sum2 NN ...'
      endif
      rntest = dcasrn()
            
      
      if(IHFLA.eq.4) then
         if(rntest.ge.(xchi0+xchi1)/sum2) then
            kchitest = 445
            elseif(rntest.ge.(xchi0)/sum2) then
            kchitest = 20443
            else
            kchitest = 10441
         endif
      
      elseif(IHFLA.eq.5) then
         if(rntest.ge.(xchi0+xchi1)/sum2) then
            kchitest = 555
            elseif(rntest.ge.(xchi0)/sum2) then
            kchitest = 20553
            else
            kchitest = 10551
         endif
      endif
      if(i23s.eq.0) then
      elseif(i23s.eq.2) then
         kchitest = 100000 + kchitest
      elseif(i23s.eq.3.and.IHFLA.eq.5) then
         kchitest = 200000 + kchitest
      endif
      
      K(NF1,2) = kchitest
      
c         write(6,*) xtotal

      

c	write(6,*) K12,k22,K(nia1,2),nia1,K(nia1,2),nia2
c	call pylist(1)
      QG1 = P(NIA1,1)**2+ P(NIA1,2)**2
      QG2 = P(NIA2,1)**2+ P(NIA2,2)**2
      SFSCAL1 = DABS(QG1)
      SFSCAL2 = DABS(QG2)
c      write(6,*) 'k12 ',k12,' k22 ',k22
c      write(6,*) 'qg1 ',qg1,' qg2 ',qg2
c new scale
      PT2Q = P(KF1,1)**2 + P(KF1,2)**2  
      xg_bar =  -qg2/(1.d0-xg)/SSS
      yph = dble(YY)
      yph_bar =  -qg1/(1.d0 - yph)/SSS

	if(sfscal1.le.sfcut) sfscal1=sfcut 
	if(sfscal2.le.sfcut) sfscal2=sfcut 

      xsi_hard(1) = (xg+yph_bar)/(yph+xg_bar)	
      xsi_hard(2) = (yph+xg_bar)/(xg+yph_bar)
c which of the quarks is closer to the proton:
      rapp = PYP(NPFIN,17)
      zz1 = 1. - (P(NF1,4)+ P(NF1,3))/2./P(2,4)/yph
      zz2 = 1. - (P(NF1,4)+ P(NF1,3))/2./P(2,4)/xg
      ptstm=dsqrt(dabs(qg1)+P(NF1,5)**2)

      If(IQBAR.eq.1) then
         stmax1 = yph*sqrt(xsi_hard(1)*SSS)
         stmax2 = xg*sqrt(xsi_hard(2)*SSS)
      elseif(IQBAR.EQ.2) then
         stmax1 = ptstm/(1.-zz1)
         stmax2 = ptstm/(1.-zz2)
      else         
         write(6,*) ' meoffchi. no valid scale selected ',IQBAR
      endif
		
	if(stmax1.ne.stmax1) then
c	   write(6,*) ' meoffchi: stmax problem ',xg,xsi_hard,SSS
	   stmax1=1d20
	endif

      stmax1 = stmax1 *scalfaf
      stmax2 = stmax2 *scalfaf
      
cccc      stmax = scale ! ala Serguei
      
      IF(FIRST) THEN
         IF(IQBAR.EQ.1) THEN
            write(6,*) ' meoffchi - scale: qq-pair '
         ELSE
            write(6,*) ' meoffchi - scale: q or q_bar '
         ENDIF
         FIRST=.FALSE.
         CALL PYNAME(k(nf1,2),CNAM)
	   write(6,*) ' meoffchi - scale for alphas ',Qscal
         write(6,*) ' meoffchi: width used for ',cnam,' Meson: ',PSI*4*pi*Xj
      ENDIF
c            write(6,*) 'meoffchi xsi ',xsi_hard,' stmax ',stmax
c            write(6,*) ' new xg,xg_bar ',xg,xg_bar
c            write(6,*) ' new yph,yph_bar ',yph,yph_bar
      CALL PYROBO(0,N,0.D0,0.D0,-CM(1)/CM(4),-CM(2)/CM(4),-CM(3)/CM(4)
     +)

	
      IF(K(NIA1,2).EQ.21.AND.K(NIA2,2).EQ.21) THEN
         if(stmax1.ge.sqrt(sss)) then
            XPQ1(0) = 0.d0
         elseif(stmax2.ge.sqrt(sss)) then
            XPQ2(0) = 0.d0
         else
            call cauniglu(K(1,2),XG1SF,SFSCAL1,STMAX1,XPQ1)
            call cauniglu(K(2,2),XG2SF,SFSCAL2,STMAX1,XPQ2)
         endif
c	   write(6,*) ' meoffchi: stmax,xpq1,xpq2 ',
c     &	   STMAX,XPQ1(0),XPQ2(0)
         WTG= dble(XPQ1(0)*XPQ2(0))
      ELSE
	   write(6,*) ' no valid process for higgs.... ',
     &   K(NIA1,2),K(NIA2,2)
         WTG= 0.d0
      ENDIF


c	write(6,*) ' meoffchi :',sum2
C new cut for on-shell DGLAP
	If(Iccfm.eq.0) then
	   If(ires(1).eq.1.and.(Qg1.gt.scale*scalfa)) SUM2 = 0 
	   If(ires(1).eq.1.and.(Qg1.gt.STMAX1**2)) SUM2 = 0 
	   If(ires(2).eq.1.and.(Qg2.gt.scale*scalfa)) SUM2 = 0 
	   If(ires(2).eq.1.and.(Qg2.gt.STMAX2**2)) SUM2 = 0 
	Endif

c      Write(6,*) ' meoffchi Wtg ',Wtg,XG1SF,XG2SF
      IF(SUM2.LT.0D0) THEN
c     write(6,*) ' SUM2 ',SUM2
         if(sum2.ne.sum2) then
            write(6,*) ' meoffchi SUM2,xsi_hard ',SUM2,xsi_hard
            write(6,*) ' meoffchi k12,k22 ',k12,k22
         endif
         SUM2 = 0D0
      ENDIF
 
      SUM2 = SUM2*WTG
c      write(6,*) ' SUM2 ',SUM2

      WT1 = SUM2
            
      if(sum2.ne.sum2) then

         write(6,*) ' meoffchi xpq(0),alpha_s',xpq1(0),xpq2(0),alpha_s
         write(6,*) ' meoffchi SUM2,xsi_hard ',SUM2,xsi_hard
         write(6,*) ' meoffchi k12,k22 ',k12,k22
         write(6,*) ' P(KIA1) ',(P(KIA1,I),I=1,5)
         write(6,*) ' P(KIA2) ',(P(KIA2,I),I=1,5)
         write(6,*) ' P(KF1) ',(P(KF1,I),I=1,5)

      endif
      CALL PYEDIT(13)
      RETURN
   20 WT1 = 0.d0
      CALL PYEDIT(13)
      RETURN
      END
      
      SUBROUTINE XSEC3(XChi0,XChi1,XChi2)
C*******************************************************************C
C      Partonic differential cross section  using Kniehl ME                        C
C*******************************************************************C
C
      IMPLICIT None
      double Precision ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/CONST2/ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      Double Precision g1,g2,pj 
      COMMON/MOMEN2/g1(4),g2(4),pj(4) 
      Double Precision XChi0,XChi1,XChi2
      
      Double Precision t1,t2,pt2,phi,M2,M,sinp,cosp
      Double Precision F0,F1,F2
      Double Precision pi
      PI=4.D0*DATAN(1.D0)
     
      t1 = g1(2)**2 + g1(3)**2
      t2 = g2(2)**2 + g2(3)**2
      pt2 = pj(2)**2 + pj(3)**2
      phi = (pt2 - t1 - t2)/2.d0/sqrt(t1*t2) 
      if(phi.le.-1.) phi=-1.
      if(phi.ge.1.) phi=1.
      
      phi = acos(phi)
      sinp =sin(phi)
      cosp =cos(phi)
      M =XJ
      M2 = M*M
      
      F0 = 2./9.*M2*(M2+pt2)**2*
     &      ((3*M2+t1+t2)*cosp+2*sqrt(t1*t2))**2/  
     &      (M2+t1+t2)**4
      F1 = 2./9.*M2*(M2+pt2)**2*
     &     ((t1+t2)**2*sinp**2 +M2*(t1+t2-2*sqrt(t1*t2)*cosp))/
     &      (M2+t1+t2)**4
      F2 = 1./3.*M2*(M2+pt2)**2/(M2+t1+t2)**4*
     &      (3*M2*M2+3*M2*(t1+t2)+4*t1*t2+
     &       (t1+t2)**2*cosp*cosp+ 2*sqrt(t1*t2)*(3*M2+2*(t1+t2))*cosp) 
      
c      write(6,*) F0,F1,F2,M2
      XChi0 = 8./3.*pi*pi*als*als*psi/M2/M2 * F0
      XChi1 = 16./3.*pi*pi*als*als*psi/M2/M2 * F1
      XChi2 = 32./45. *pi*pi*als*als*psi/M2/M2 * F2
      
      if(xchi0.ne.xchi0) then
         write(6,*) 'xchi0 ',xchi0,F0,t1,t2,phi,m2
         write(6,*) 'test',(pt2 - t1 -t2)/2./sqrt(t1*t2),pt2
      endif 
      if(xchi1.ne.xchi1) then
         write(6,*) 'xchi1 ',xchi1,F1
      endif
      if(xchi2.ne.xchi2) then
         write(6,*) 'xchi2 ',xchi2,F2
      endif
c      write(6,*) XChi1,COLOR,FACTOR
      RETURN
      END

      SUBROUTINE XSEC2(XChi0,XChi1,XChi2)
C*******************************************************************C
C      Partonic differential cross section                          C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Parameter(GeVmb=0.4D0, GeVmub=4.D2, GeVnb=4.D5, GeVpb=4.D8)
      COMMON/TYPE2/ITYPE,MODEL,IB1,IB2,NJ,NL
      COMMON/CONST2/ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/LOOPJ2/AMPJ(2,-1:1,-1:1)
      COMMON/LOOPX2/AMPX0(2),AMPX1(2,-1:1),AMPX2(2,-2:2)
      SAVE
      PI=4.D0*DATAN(1.D0)
C...Gluon polarization vectors
      CALL GAUGEG2
C...Meson polarization vectors
      IF(NJ.EQ.1) CALL GAUGEJ2
C...Evaluating the Matrix Elements
      CALL FEYNJ2
      IF(NL.EQ.1 .AND. NJ.EQ.1) CALL CHIJ 
      XChi0=0.D0
      XChi1=0.D0
      XChi2=0.D0
C
      COLOR=0.D0    !Color factors
      IF(ITYPE.EQ.10) COLOR=  8./3.      ! glu+glu->Eta(1)
      IF(ITYPE.EQ.12) COLOR=  8./3.      ! glu+glu->Chi(1)
      IF(ITYPE.EQ.30) COLOR= 12.*QC2*QC2 ! gam+gam->Eta(1)
      IF(ITYPE.EQ.32) COLOR= 12.*QC2*QC2 ! gam+gam->Chi(1)
      IF(ITYPE.EQ.40) COLOR=  5./6.      ! glu+glu->Eta(8)
      IF(ITYPE.EQ.41) COLOR=  3./2.      ! glu+glu->Psi(8)   
      IF(ITYPE.EQ.42) COLOR=  5./6.      ! glu+glu->Chi(8)
      IF(ITYPE.EQ.50) COLOR=  2.*QC2     ! gam+glu->Eta(8)
      IF(ITYPE.EQ.52) COLOR=  2.*QC2     ! gam+glu->Chi(8)
      FACTOR=0.D0   !Coupling constants
      IF((IB1*IB2).EQ.64) FACTOR = (4*pi)**2 *PSI     *ALS**2 /64.
      IF((IB1*IB2).EQ.8)  FACTOR = (4*pi)**2 *PSI *ALP*ALS    /8.
      IF((IB1*IB2).EQ.1)  FACTOR = (4*pi)**2 *PSI *ALP**2
C                   !Extra compensation for matrix elements
      IF(NJ.EQ.0 .AND. NL.EQ.0) FACTOR =FACTOR*16.
      IF(NJ.EQ.1 .AND. NL.EQ.0) FACTOR =FACTOR*16.*XC2
      IF(NJ.EQ.1 .AND. NL.EQ.1) FACTOR =FACTOR*48./XJ2
C
C... Summing up
C 1P0 state
      IF(NJ.EQ.0 .AND. NL.EQ.0) XChi0 = AMPJ(1,0,0)*AMPJ(1,0,0)
C 3S1 state
      IF(NJ.EQ.1 .AND. NL.EQ.0) THEN
      DO 11 j=-1,1
  11  XChi1 =XChi1 +AMPJ(1,j,0)*AMPJ(1,j,0) +AMPJ(2,j,0)*AMPJ(2,j,0)
      ENDIF
C 3PJ states
      IF(NJ.EQ.1 .AND. NL.EQ.1) THEN
      XChi0 =       AMPX0(1)*AMPX0(1)     +AMPX0(2)*AMPX0(2)
      DO 13 k=-1,1
  13  XChi1 =XChi1 +AMPX1(1,k)*AMPX1(1,k) +AMPX1(2,k)*AMPX1(2,k)
      DO 15 k=-2,2
  15  XChi2 =XChi2 +AMPX2(1,k)*AMPX2(1,k) +AMPX2(2,k)*AMPX2(2,k)
      ENDIF
C
      XChi0 = XChi0*COLOR*FACTOR !*GeVnb
      XChi1 = XChi1*COLOR*FACTOR !*GeVnb
      XChi2 = XChi2*COLOR*FACTOR !*GeVnb
c      write(6,*) XChi1,COLOR,FACTOR
      RETURN
      END

      SUBROUTINE GAUGEJ2
C********************************************************************C
C       Polarization vector of outgoing meson                        C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION OX(4),OY(4),OZ(4)
      COMMON/CONST2/ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/MOMEN2/g1(4),g2(4),pj(4)
      COMMON/SPINJ2/EJ(2,4,-1:1)
      SAVE
      sq2=sqrt(2.)
      QJ=dsqrt(pj(3)**2 +pj(2)**2 +pj(1)**2)
      Qt=dsqrt(pj(3)**2 +pj(2)**2)
      CT=pj(1)/QJ
      ST=Qt/QJ
      CP=pj(2)/Qt
      SP=pj(3)/Qt
C... Longitudinal polarization OZ
      OZ(1)= pj(1)*pj(4)/QJ/XJ
      OZ(2)= pj(2)*pj(4)/QJ/XJ
      OZ(3)= pj(3)*pj(4)/QJ/XJ
      OZ(4)= QJ/XJ
C... Transverse polarization OX
      OX(1)=-ST
      OX(2)= CT*CP
      OX(3)= CT*SP
      OX(4)= 0.D0
C... Transverse polarization OY
      OY(1)= 0.D0
      OY(2)=-SP
      OY(3)= CP
      OY(4)= 0.D0
C... Alltogether
      DO 4 j=1,4
      EJ(1,j, 0)= OZ(j)
      EJ(2,j, 0)= 0.D0
      EJ(1,j,-1)= OX(j)/sq2
      EJ(2,j,-1)= OY(j)/sq2
      EJ(1,j, 1)= OX(j)/sq2
      EJ(2,j, 1)=-OY(j)/sq2
   4  CONTINUE
      RETURN
      END

      SUBROUTINE GAUGEG2
C*******************************************************************C
C      Photon and Gluon polarization vectors                        C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/MOMEN2/g1(4),g2(4),pj(4)
      COMMON/POLAR2/a1(4),a2(4)
      SAVE
      g1T =dsqrt(g1(2)**2+g1(3)**2)
      g2T =dsqrt(g2(2)**2+g2(3)**2)
      a1(1)=0.D0
      a1(4)=0.D0
      a2(1)=0.D0
      a2(4)=0.D0
      sq=sqrt(2.)
      DO 2 i=2,3
      a1(i)=g1(i)/g1T
      a2(i)=g2(i)/g2T
    2 Continue
C     
      RETURN !Comment out this line for gauge invariance tests
      DO 9 i=1,4
      a1(i)=g1(i) 
      a2(i)=g2(i)
    9 Continue
      RETURN
      END

      SUBROUTINE FEYNJ2
C*******************************************************************C
C      Matrix Elements for all Feynman diagrams                     C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION g3(4)
      COMMON/TYPE2/ITYPE,MODEL,IB1,IB2,NJ,NL
      COMMON/CONST2/ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/MOMEN2/g1(4),g2(4),pj(4)
      COMMON/POLAR2/a1(4),a2(4)
      COMMON/LOOPJ2/AMPJ(2,-1:1,-1:1)
      COMMON/GMUNU2/DF(4,4),DC(4),EP(4,4,4,4)
      COMMON/SPINJ2/EJ(2,4,-1:1)
      SAVE
      W1=DOTB(g1,g1)
      W2=DOTB(g2,g2)
      Z0=(W1+W2-XJ2)/2. ! = -DOTB(g1,g2)
C...Clean main array
      DO 1 j=-1,1
      DO 1 l=-1,1
      AMPJ(1,j,l)= 0.D0
      AMPJ(2,j,l)= 0.D0
   1  CONTINUE
      DO 2 i=1,4
   2  g3(i) =0.5*(g1(i)-g2(i))
C
C... Eta = 1S0 state
      IF(NJ.NE.0 .OR. NL.NE.0) GOTO 110
      DO 105 m1=1,4
      DO 105 m2=1,4
      DO 105 k1=1,4
      DO 105 k2=1,4
      !Diagram << 1 >> = << 2 >>
      TERM0 = EP(k1,k2,m1,m2)*g1(k1)*g2(k2)*DC(k1)*DC(k2)/2./Z0
      AMPJ(1,0,0) = AMPJ(1,0,0) + TERM0*a1(m1)*a2(m2)*DC(m1)*DC(m2)
 105  CONTINUE
C
C... Psi = 3S1 state
 110  IF(NJ.NE.1 .OR. NL.NE.0) GOTO 120
      DO 115 m1=1,4
      DO 115 m2=1,4
      DO 115 mj=1,4      
      !Diagram << 1 >> = - << 2 >> = Diagram << 3 >>/2 
      TERM0 = DF(m1,m2)*g3(mj) - DF(m1,mj)*pj(m2) + DF(m2,mj)*pj(m1)
      TERM0 = TERM0*(1./Z0 +2./XJ2)*a1(m1)*a2(m2)*DC(m1)*DC(m2)
      DO 115 j=-NJ,NJ
      AMPJ(1,j,0) = AMPJ(1,j,0) + TERM0*EJ(1,mj,j)*DC(mj)
      AMPJ(2,j,0) = AMPJ(2,j,0) + TERM0*EJ(2,mj,j)*DC(mj)
 115  CONTINUE
C
C... Chi = 3PJ state
 120  IF(NJ.NE.1 .OR. NL.NE.1) GOTO 199
      DO 125 m1=1,4
      DO 125 m2=1,4
      DO 125 mj=1,4
      DO 125 ml=1,4
      !Zeroth order in relative momentum, Diagram << 1 >> = - << 2 >>
      TERM0 = - DF(m1,m2)*g3(mj) + DF(m1,mj)*g1(m2) - DF(m2,mj)*g2(m1)
      !Linear terms in relative momentum, Diagram << 1 >> = << 2 >>
      TERM1 = - DF(m1,mj)*DF(m2,ml)*W1  + DF(m1,mj)*DF(m2,ml)*Z0
     .        + DF(m1,mj)*pj(m2)*g3(ml) - DF(m1,ml)*DF(m2,mj)*W2
     .        + DF(m1,ml)*DF(m2,mj)*Z0  - DF(m1,ml)*pj(m2)*g3(mj)
     .        - DF(m2,mj)*pj(m1)*g3(ml) + DF(m2,ml)*pj(m1)*g3(mj)
      Z1  = -2.*g3(ml) 
      SUM = TERM1/Z0 - 2*XC2*TERM0*Z1/(Z0)**2
      SUM = SUM*a1(m1)*a2(m2)*DC(m1)*DC(m2)
      DO 125 j=-NJ,NJ
      DO 125 l=-NL,NL
      AMPJ(1,j,l) = AMPJ(1,j,l)                           !real part
     .            + SUM*EJ(1,mj,j)*DC(mj)*EJ(1,ml,l)*DC(ml)
     .            - SUM*EJ(2,mj,j)*DC(mj)*EJ(2,ml,l)*DC(ml)
      AMPJ(2,j,l) = AMPJ(2,j,l)                           !imag part
     .            + SUM*EJ(1,mj,j)*DC(mj)*EJ(2,ml,l)*DC(ml)
     .            + SUM*EJ(2,mj,j)*DC(mj)*EJ(1,ml,l)*DC(ml)
 125  CONTINUE
 199  RETURN
      END


      SUBROUTINE CHIJ
C*******************************************************************C
C     Converts (S,L) amplitudes into (J,mj) amplitudes              C
C       using the Clebsch-Gordan coefficients                       C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/LOOPJ2/AMPJ(2,-1:1,-1:1)
      COMMON/LOOPX2/AMPX0(2),AMPX1(2,-1:1),AMPX2(2,-2:2)
      SAVE
      sq2=sqrt(2.)
      sq3=sqrt(3.)
      sq6=sqrt(6.)
      DO 1  i=1,2 !real and imag parts
      AMPX2(i, 2)=AMPJ(i, 1, 1)
      AMPX2(i, 1)=AMPJ(i, 1, 0)/sq2 + AMPJ(i, 0, 1)/sq2
      AMPX2(i, 0)=AMPJ(i, 1,-1)/sq6 + AMPJ(i,-1, 1)/sq6
     .           -AMPJ(i, 0, 0)*sq2/sq3       !sign changed here
      AMPX2(i,-1)=AMPJ(i,-1, 0)/sq2 + AMPJ(i, 0,-1)/sq2
      AMPX2(i,-2)=AMPJ(i,-1,-1)
      AMPX1(i, 1)=AMPJ(i, 1, 0)/sq2 - AMPJ(i, 0, 1)/sq2
      AMPX1(i, 0)=AMPJ(i, 1,-1)/sq2 - AMPJ(i,-1, 1)/sq2
      AMPX1(i,-1)=AMPJ(i,-1, 0)/sq2 - AMPJ(i, 0,-1)/sq2
      AMPX0(i)   =AMPJ(i, 1,-1)/sq3 + AMPJ(i,-1, 1)/sq3
     .           +AMPJ(i, 0, 0)/sq3           !sign changed here
    1 CONTINUE
      RETURN
      END

      SUBROUTINE METRIC2
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/GMUNU2/DF(4,4),DC(4),EP(4,4,4,4)
      DO 1 I=1,4
      DO 1 J=1,4
      DF(I,J)= 0.D0
      DO 1 M=1,4
      DO 1 N=1,4
      EP(I,J,M,N)= 0.D0
    1 CONTINUE !P(i) =(Beam, P_t, P_t, Energy)
      DO 2 I=1,3
      DF(I,I)=-1.D0
    2 DC(I)  =-1.D0
      DF(4,4)= 1.D0
      DC(4)  = 1.D0
      EP(1,2,3,4)= 1.D0
      EP(1,2,4,3)=-1.D0
      EP(1,3,2,4)=-1.D0
      EP(1,3,4,2)= 1.D0
      EP(1,4,2,3)= 1.D0
      EP(1,4,3,2)=-1.D0
      EP(2,1,3,4)=-1.D0
      EP(2,1,4,3)= 1.D0
      EP(2,3,1,4)= 1.D0
      EP(2,3,4,1)=-1.D0
      EP(2,4,1,3)=-1.D0
      EP(2,4,3,1)= 1.D0
      EP(3,1,2,4)= 1.D0
      EP(3,1,4,2)=-1.D0
      EP(3,2,1,4)=-1.D0
      EP(3,2,4,1)= 1.D0
      EP(3,4,1,2)= 1.D0
      EP(3,4,2,1)=-1.D0
      EP(4,1,2,3)=-1.D0
      EP(4,1,3,2)= 1.D0
      EP(4,2,1,3)= 1.D0
      EP(4,2,3,1)=-1.D0
      EP(4,3,1,2)=-1.D0
      EP(4,3,2,1)= 1.D0
      RETURN
      END


      DOUBLE PRECISION FUNCTION EPS2(q1,q2,q3,q4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION q1(4),q2(4),q3(4),q4(4)
      G1 = q1(2)*q2(3)*q3(4) +q1(3)*q2(4)*q3(2) +q1(4)*q2(2)*q3(3)
     .    -q1(4)*q2(3)*q3(2) -q1(2)*q2(4)*q3(3) -q1(3)*q2(2)*q3(4)
      G2 = q1(1)*q2(3)*q3(4) +q1(3)*q2(4)*q3(1) +q1(4)*q2(1)*q3(3)
     .    -q1(4)*q2(3)*q3(1) -q1(1)*q2(4)*q3(3) -q1(3)*q2(1)*q3(4)
      G3 = q1(1)*q2(2)*q3(4) +q1(2)*q2(4)*q3(1) +q1(4)*q2(1)*q3(2)
     .    -q1(4)*q2(2)*q3(1) -q1(1)*q2(4)*q3(2) -q1(2)*q2(1)*q3(4)
      G4 = q1(1)*q2(2)*q3(3) +q1(2)*q2(3)*q3(1) +q1(3)*q2(1)*q3(2)
     .    -q1(3)*q2(2)*q3(1) -q1(1)*q2(3)*q3(2) -q1(2)*q2(1)*q3(3)
      EPS= -G1*q4(1) +G2*q4(2) -G3*q4(3) +G4*q4(4)
      RETURN
      END


