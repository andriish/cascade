      SUBROUTINE meoffjpsi(WT1)
      IMPLICIT None
C
C    PHOTON GLUON ----> JPSI GLUON
C
C         P1(G)-----//////--------Q1(JPSI)
C                  ////////
C         P2(PH)-----//////--------Q2(G)
C
C  ( calculation of A. Lipatov, N.Zotov 
C    new ala S. Baranov
C
C    k_photon = y * e + kt_gamma
C    k_gluon = x * p + kt_gluon
C    where y = a1 + a2, x = b1 + b2 and
C    a1,a2,b1,b2 - standard Sudakov variables 

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
      DOUBLE PRECISION K12,K22
      DOUBLE PRECISION PT2Q
      Double Precision XPQ1(-6:6),XPQ2(-6:6),SFSCAL1,SFSCAL2 
      Double Precision STMAX1,STMAX2
      DOUBLE PRECISION DSHN(4),DSH(4),dshnf(4)
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
      Double Precision QF2,alpha_s,qg2,qg1,sh,t,u,wt1,sum2,wtg
      Double Precision SH_NEW,pcm_fin,cos_fin,cphi,cphio,test
      Double Precision sh_old,shf,check,dphi,dtheta
      Double Precision yph,yph_bar,xsi_hard(2),epsilon
      Double Precision alph_em,pcm_new,alphas,nq
      Double Precision eq,width,wf

      Double Precision m2,ValT,ValL,y,xq2,xkt2
c commons for S.baranov
      Integer ITYPE,MODEL,IB1,IB2
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2
      double Precision ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/CONST/ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      Double Precision g1,g2,pj,g3
      COMMON/MOMEN/g1(4),g2(4),pj(4),g3(4)
      Double Precision v1,v2
      COMMON/MUONS/v1(4),v2(4)
      Double Precision Xtotal
            
      
      Double Precision phi,pt2,pt2g

      Double Precision dot,dot1,dcasrn
      Double Precision Del,Nc
      Double Precision rap1,rap2,rapp,ptstm1,ptstm2,zz1,zz2
      Integer npfin
      Integer IQBAR
      Integer  ikincut,Ipgg,ns_sel
      Double Precision  QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
	
	Double Precision QQg,lambda
	Double Precision Qscal
	Integer Iqqbar
	Common/GLUDAT2/Qscal,Iqqbar
c      Data IQBAR/1/
c      Data IQBAR/2/

      Double Precision PT2T
      Integer Nacc
      
      Data Nc/3.d0/
      Logical ltest
c      data ltest/.true./
      data ltest/.false./
      Logical first
      data first/.true./
      SAVE STMAX1,STMAX2, XPQ2  
      CHARACTER * 12 CNAM
      Integer ksave(4000),nold
      Integer Izotov
      data Izotov/0/
      
c     call pylist(1)
      IQBAR = Iqqbar
      ipsiel1 = 0
      ipsiel2 = 0
c      IQbar = 2
      eq = DFLOAT(PYCHGE(ihfla))/3.D0
c      if(ihfla.eq.3) then
c        width = 1.26D-6 ! omega
c        elseif(ihfla.eq.4) then
c        width = 5.55D-6 ! j/psi
c        elseif(ihfla.eq.5) then
c        width = 1.34d-6 ! upsilon
c        else
c        write(6,*) ' meoffjpsi: wrong flavor selected, program stopped '
c        stop
c      endif
      if(K(nf1,2).eq.333) then
        width = 1.26D-6 ! phi
        elseif(K(nf1,2).eq.443) then
        width = 5.55D-6 ! j/psi
        elseif(K(nf1,2).eq.100443) then
        width = 2.38D-6 ! j/psi
c        elseif(K(nf1,2).eq.300443) then
c        width = 5.55D-6 ! j/psi
        elseif(K(nf1,2).eq.553) then
        width = 1.34d-6 ! upsilon
        elseif(K(nf1,2).eq.100553) then
        width = 0.612d-6 ! upsilon
        elseif(K(nf1,2).eq.200553) then
        width = 0.443d-6 ! upsilon
        else
        write(6,*) ' meoffjpsi: wrong flavor selected, program stopped '
        stop
      endif
      IF(INTER.EQ.2) THEN
         write(6,*) ' INTER = 2 not implemented for meoffjpsi'
         STOP
      ENDIF
      NPFIN = N
c	call pylist(1)
      Sh = AM(1)**2 + AM(2)**2 + 2.D0 * DOT1(NF1,NF2)
      WT1 = 0.d0
      IF(Sh.LE.0.) RETURN
c boost  cms
ccc      CALL PYROBO(0,N,0.D0,0.D0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
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
      SH_old = DOT(DSH,DSH)
      if(ltest) then
         SH_NEW = DOT(DSH,DSH)
         goto 10
      endif

      if(iccfm.eq.1) then
      
c      P(KIA1,3) = 0.5d0*(P(KIA1,3)-P(KIA1,4))
c      P(KIA1,4) = DABS(P(KIA1,3))
c      P(KIA2,3) = 0.5d0*(P(KIA2,3)+P(KIA2,4))
c      P(KIA2,4) = DABS(P(KIA2,3))


      DO I=1,4
         DSHN(I)= P(KIA1,I) + P(KIA2,I)
      ENDDO
      SH_NEW = DOT(DSHN,DSHN)
      IF(SH_NEW.LE.AM(1)**2) GOTO 20

c      call pylist(1)

      elseif(iccfm.eq.0) then

      CALL PYROBO(0,0,0.D0,0.D0,-DSH(1)/DSH(4), -DSH(2)/DSH(4), -DSH(3)
     +/DSH(4))
      DPHI = PYANGL(P(NIA1,1),P(NIA1,2))
      call PYROBO(0,0,0.D0,-dphi,0.d0,0.d0,0.d0)
      DTHETA = PYANGL(P(NIA1,3),P(NIA1,1))
      call PYROBO(0,0,-DTHETA,0.D0,0.d0,0.d0,0.d0)
      IF(ICCFM.EQ.0) Then
	   If(ires(1).eq.1) then
            P(KIA1,1) = 0D0
            P(KIA1,2) = 0D0
            P(KIA1,5) = 0
         endif
	   If(ires(2).eq.1) then
            P(KIA2,1) = 0D0
            P(KIA2,2) = 0D0
            P(KIA2,5) = 0
         endif
      ENDIF
      P(KIA1,3) = 0.5d0*sqrt(sh_old)
      P(KIA1,4) = DABS(P(KIA1,3))
      P(KIA2,3) = -0.5d0*sqrt(sh_old)
      P(KIA2,4) = DABS(P(KIA2,3))
c      call pylist(1)
            
      do i=1,4
      dshnf(i) = P(kf1,i)+P(kf2,i)
      enddo
      shf = dot(dshnf,dshnf)

      check = abs(P(kf1,1)+P(kf2,1)+
     &  P(kf1,2)+P(kf2,2) +
     &  P(kf1,3)+P(kf2,3) +
     &  P(kf1,4)+P(kf2,4)-P(KIA2,4)-P(KIA1,4))
      if(check.gt.1e-5) then 
      write(6,*) 'shat ',sqrt(sh_old),sqrt(sh_new),sqrt(shf),check
      write(6,*) ' meoffjpsi momentum checks ',P(kf1,1)+P(kf2,1),
     &  P(kf1,2)+P(kf2,2), 
     &  P(kf1,3)+P(kf2,3), 
     &  P(kf1,4)+P(kf2,4)-P(KIA2,4)-P(KIA1,4)
      endif

      
      call PYROBO(0,0,DTHETA,0.D0,0.d0,0.d0,0.d0)
      call PYROBO(0,0,0.D0,dphi,0.d0,0.d0,0.d0)
      CALL PYROBO(0,0,0.D0,0.D0,DSH(1)/DSH(4) ,DSH(2)/DSH(4), DSH(3)/
     +DSH(4))
c      call pylist(1)
cccc      write(6,*) ' here not enter',iccfm
      endif
   10 continue
      K12 = DOT1(KIA1,KIA1)
      K22 = DOT1(KIA2,KIA2)

      M2 = AM(1)**2
      T = AM(1)**2 - 2.D0 * DOT1(KIA1,KF1) + K12
      U = AM(1)**2 - 2.D0 * DOT1(KIA2,KF1) + K22
      Sh = K12 + 2.D0 * DOT1(KIA1,KIA2) + K22
      test=Sh+t+u-AM(1)**2-AM(2)**2-K12-K22

c new scale
      PT2Q = P(KF1,1)**2 +P(KF1,2)**2
      pt2 = pt2q
	
c	call pylist(1)
      if(abs(test).gt.1e-5) then
         write(6,*) ' sh+t+u ,m**2 ',Sh_new+t+u,AM(1)**2+AM(2)**2+K12+
     +   K22,AM(1)
         write(6,*) ' sh+t+u-m**2 ',Sh_new+t+u-AM(1)**2-AM(2)**2-K12-
     +   K22
      endif
c new scale
c check on scale:
      SCALE = PT2Q+M2
      IF(IQ2.EQ.1) THEN
         SCALE = M2
      ELSEIF(IQ2.EQ.2) THEN
         SCALE = SH_NEW
      ELSEIF(IQ2.EQ.3) THEN
c         SCALE = 4D0*M2 + PT2Q
         SCALE = M2 + PT2Q
      ELSEIF(IQ2.EQ.4) THEN
         SCALE = Q2
      ELSEIF(IQ2.EQ.5) THEN
         SCALE = Q2 + PT2Q + M2
      ELSEIF(IQ2.EQ.6) THEN
         SCALE = DABS(K12)+DABS(K22)
      ELSE
         WRITE(6,*) ' NO VALID Q2 SCALE. STOP'
         STOP
      ENDIF
	
c	write(6,*) ' scales stmax = ',stmax,' sqrt(scale) = ',sqrt(scale)
c	stmax = sqrt(scale)*scalfa
c multiplication with Qscal done in alphas.F      
      Dscale = scale
      alpha_s = ALPHAS(sqrt(dscale))

      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)
      
c multiply xsection with wavefct
      wf=width/(1d0-16d0/3d0*ALPHAS(sqrt(m2))/pi)*m2/16d0/pi/eq/eq/alph_em**2

      y = dble(YY)
      if(izotov.eq.1) then 
c here is the calculation of zotov et al.      
         phi = P(KIA2,1)*P(KF1,1)+P(KIA2,2)*P(KF1,2)
         phi= phi/sqrt(P(KIA2,1)**2+P(KIA2,2)**2)
         phi= phi/sqrt(P(KF1,1)**2+P(KF1,2)**2)
         phi=acos(phi)
         xq2 =  P(KIA1,1)**2 +P(KIA1,2)**2
         xkt2 = P(KIA2,1)**2 +P(KIA2,2)**2
         pt2  = P(KF1,1)**2 +P(KF1,2)**2
         pt2g = P(KF2,1)**2 +P(KF2,2)**2
c      if(pt2g.lt.1) goto 20
         call mepsi(t,u,m2,xq2,xkt2,pt2,phi,eq,alph_em,alpha_s,valT,valL)
c      write(6,*) valt,vall
         if(abs(k12).gt.0.) then 
            epsilon = 2*(1d0-yph)/(1d0-yph+yph**2/2d0)
            sum2 = valT - epsilon *valL
         else
            sum2 = valT
         endif
c multiply xsection with wavefct
c         wf=width/(1d0-16d0/3d0*alpha_s/pi)*m2/16d0/pi/eq/eq/alph_em**2
c         write(6,*) wf,eq,alph_em,m2
      
cccc         sum2 = sum2 * 0.0876 
         sum2 = sum2 * wf
      else
c here is the calculation of baranov et al. 
c The value of PSI is calculated from the experimentally known
c leptonic width; Gamma_{ll} is proportional to |Psi(0)|^2.
c When relying upon the Leading-Order formula for the process
c J/psi --> l+ l-  one gets |Psi(0)|^2 \simeq 0.038 GeV^3.
c When including the NLO strong corrections one gets
c approximately twice as large value. The difference between
c 0.08 and 0.0876 comes from the different choice of alpha_s
c in the strong corrections.

c         write(6,*) ' checking frame',M2,KIA1,KIA2,KF1,KF2
c gymnastics for polasisation of decay leptons
c 	   write(6,*) 'before jpsi decay '
c         call pylist(1)
         if(ipsipol.ge.1) then
            do i=1,n
              ksave(I) = k(i,1)
              k(i,1)= 21
            end do 
            k(nf1,1)= 1
            nold = n
            call pyexec
            do i=1,nold
            if(i.eq.nf1) then 
              else
              k(i,1)=ksave(I)
            endif
            enddo 
c 	   write(6,*) 'after jpsi decay '
c         call pylist(1)
         endif
         
         ALS = alpha_s
         ALP = ALPH_EM

cc         M2 = 3.1*3.1 ! check this.... this is only for J/psi
         XJ = sqrt(M2)
         XC = XJ/2
         XC2 = XC**2
         XJ2 = XJ**2
         QC2 = eq**2 
ccc         PSI = 0.08/XJ ! check this.... this is only for J/psi
ccc         PSI = 0.0876/XJ ! check this.... this is only for J/psi
         PSI = wf/XJ
c         write(6,*) wf
         if(ICCFM.EQ.1) then 
           MODEL=1
           else
           Model = 0
         endif
         IF(ires(1).eq.0.and.Ires(2).eq.0) then
            itype = 3 ! 3 = lepton*lepton [gam + gam -> Jpsi + gam]
c            itype = 8 ! 8 = lepton*lepton [gam + gam -> Jpsi8 + glu]
            ib1 = 1
            ib2 = 1
         elseif(ires(1).eq.0.and.Ires(2).eq.1) then
            itype = 2 ! 2 = lepton*proton [gam + glu -> Jpsi + glu]
c            itype = 6 ! 6 = lepton*proton [gam + glu -> Jpsi8 + gam]
c            itype = 7 ! 7 = lepton*proton [gam + glu -> Jpsi8 + glu]
            ib1 = 1
            ib2 = 8
         elseif(ires(1).eq.1.and.Ires(2).eq.1) then
c            itype = 0 ! 0 = proton*proton [glu + glu -> Jpsi + gam]
            itype = 1 ! 1 = proton*proton [glu + glu -> Jpsi + glu]
c            itype = 4 ! 4 = proton*proton [glu + glu -> Jpsi8 + gam]
c            itype = 5 ! 5 = proton*proton [glu + glu -> Jpsi8 + glu]
            ib1 = 8
            ib2 = 8
         endif
c include a very small kt for the calc of the polarisation
c othersise errors produced         
         g1(1) = P(KIA1,3)
         g1(2) = P(KIA1,1)+1.e-10*P(KIA1,4)
         g1(3) = P(KIA1,2)+1.e-10*P(KIA1,4)
         g1(4) = P(KIA1,4)
         g2(1) = P(KIA2,3)+1.e-10*P(KIA2,4)
         g2(2) = P(KIA2,1)+1.e-10*P(KIA2,4)
         g2(3) = P(KIA2,2)
         g2(4) = P(KIA2,4)

         pj(1) = P(KF1,3)
         pj(2) = P(KF1,1)
         pj(3) = P(KF1,2)
         pj(4) = P(KF1,4)
         g3(1) = P(KF2,3)
         g3(2) = P(KF2,1)
         g3(3) = P(KF2,2)
         g3(4) = P(KF2,4)
         If(ipsipol.ge.1) then 
           v1(1) = P(N-1,3)
           v1(2) = P(N-1,1)
           v1(3) = P(N-1,2)
           v1(4) = P(N-1,4)
           v2(1) = P(N,3)
           v2(2) = P(N,1)
           v2(3) = P(N,2)
           v2(4) = P(N,4)
         Endif
         CALL METRIC
         Call xsec(xtotal)
         sum2 = xtotal
         

         if(sum2.ne.sum2) then
           sum2=-9999.
           write(6,*) ' sum2 NN ...'
         endif
c         write(6,*) xtotal
      endif
      IF(FIRST) THEN
	   XG2SF = 0.001
	   SFSCAL2 = 2.
	   STMAX2 = 5.
         call cauniglu(K(2,2),XG2SF,SFSCAL2,STMAX2,XPQ2)
         IQBAR = Iqqbar
         IF(IQBAR.EQ.1) THEN
            write(6,*) ' meoffjpsi - scale: qq-pair '
         ELSE
            write(6,*) ' meoffjpsi - scale: q or q_bar '
         ENDIF
	      write(6,*) ' meoffjpsi - scale for alphas ',Qscal
         CALL PYNAME(k(nf1,2),CNAM)
         write(6,*) ' meoffjpsi: width used for ',cnam,' Meson: ',width
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
      xg_bar =  -qg2/(1.d0-xg)/SSS
      yph = dble(YY)
      yph_bar =  -qg1/(1.d0 - yph)/SSS

	if(sfscal1.le.sfcut) sfscal1=sfcut 
	if(sfscal2.le.sfcut) sfscal2=sfcut 
      xsi_hard(1) = (xg+yph_bar)/(yph+xg_bar)	
      xsi_hard(2) = (yph+xg_bar)/(xg+yph_bar)
      if(xsi_hard(2).le.0) then
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
         write(6,*) ' meoffjpsi. no valid scale selected ',IQBAR
      endif

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
      IF(K(NIA1,2).EQ.22.AND.K(NIA2,2).EQ.21) THEN
         WTG= dble(XPQ2(0))
      ELSEIF(K(NIA1,2).EQ.21.AND.K(NIA2,2).EQ.21) THEN
c	   Sum2 = 1
  	   wtg =dble(XPQ1(0)*XPQ2(0))  
      Else
         write(6,*) ' meoffjpsi: wrong subprocesses selected :',
     &      K(NIA1,2),K(NIA2,2)
         write(6,*) ' meoffjpsi: program stopped '
         stop
      ENDIF

C new cut for on-shell DGLAP
	If(Iccfm.eq.0) then
	   If(ires(1).eq.1.and.(Qg1.gt.(PT2Q+M2))) SUM2 = 0 
	   If(ires(1).eq.1.and.(Qg1.gt.scale*scalfa)) SUM2 = 0 
	   If(ires(1).eq.1.and.(Qg1.gt.STMAX1**2)) SUM2 = 0 
	   If(ires(2).eq.1.and.(Qg2.gt.(PT2Q+M2))) SUM2 = 0 
	   If(ires(2).eq.1.and.(Qg2.gt.scale*scalfa)) SUM2 = 0 
	   If(ires(2).eq.1.and.(Qg2.gt.STMAX2**2)) SUM2 = 0 
	Endif

      IF(SUM2.LT.0D0) THEN
c         write(6,*) ' SUM2 ',SUM2
         if(sum2.ne.sum2) then
            write(6,*) ' meoffjpsi SUM2,xsi_hard ',SUM2,
     &        xsi_hard(1),xsi_hard(2)
            write(6,*) ' meoffjpsi k12,k22 ',k12,k22
            write(6,*) ' meoffjpsi t,u ',t,u
         endif
         SUM2 = 0D0
      ENDIF
c apply pt cuts here
      PT2T = PYP(KF2,9)
      CALL CUTG(PT2T,NACC)
      IF(NACC.EQ.0) SUM2 = 0
 
      SUM2 = SUM2*WTG
   
      WT1 = SUM2
      
      

      if(sum2.ne.sum2) then

         write(6,*) ' meoffjpsi xpq(0),alpha_s,alph_em',xpq1(0),xpq2(0),
     +      alpha_s,alph_em
         write(6,*) ' meoffjpsi SUM2,xsi_hard ',SUM2,
     &      xsi_hard(1),xsi_hard(2)
         write(6,*) ' meoffjpsi k12,k22 ',k12,k22
         write(6,*) ' meoffjpsi t,u ',t,u
         write(6,*) ' P(KIA1) ',(P(KIA1,I),I=1,5)
         write(6,*) ' P(KIA2) ',(P(KIA2,I),I=1,5)
         write(6,*) ' P(KF1) ',(P(KF1,I),I=1,5)
         write(6,*) ' P(KF2) ',(P(KF2,I),I=1,5)

      endif
      CALL PYEDIT(14)
      if(ipsipol.eq.1) then 
        ipsiel1 = n-1
        ipsiel2 = n
      endif
      RETURN
   20 WT1 = 0.d0
      CALL PYEDIT(14)
      RETURN
      END
      
      subroutine mepsi(t,u,m2,q2,kt2,pt2,phi,eq,alpe,alps,valT,valL)
      implicit none
      
      Double Precision m2,eC,eq,NT,ValT,NL,ValL
      
      Double Precision kt2,kt4,kt6,kt8,kt10,pt2
      Double Precision q2,q4,q6,q8
      Double Precision m4,m6,m8,m10,m12
      Double Precision t,t2,u,u2,u3,t3,t4,u4
      Double Precision phi
      Double Precision alpe,alps,pi
      
      double precision sqr,pow
      
      eC = eq
      PI=4.D0*DATAN(1.D0)
      
      m4 = m2**2
      m6 = m4*m2
      m8 = m6*m2
      m10 = m8*m2
      m12 = m10*m2
      
      q4 = q2*q2
      q6 = q4*q2
      q8 = q6*q2
      
      kt4 = kt2*kt2
      kt6 = kt4*kt2
      kt8 = kt6*kt2
      kt10 = kt8*kt2
      
      t2 = t*t
      t3 = t2*t
      t4 = t3*t
      
      u2 = u*u
      u3 = u2*u
      u4 = u3*u
      
c      write(6,*) phi,m2,Q2,kt2,t,u
c      write(6,*) pi,alph_em,ec,alpha_s
c      write(6,*) ' sqrts ',pt2*kt2
c      write(6,*) 'us ',u,u2,u3,u4
c      write(6,*) 'ts ',t,t2,t3,t4
c      write(6,*) 'kts ',kt2,kt4,kt6,kt8,kt10
c      write(6,*) 'qs ',q2,q4,q6,q8


      NT = 64.0*4.0*Pi*alpe*(eC)**2*(4.0*Pi*alps)**2 
        
	ValT = NT*(kt2*(2*kt10*Q4 + kt8*Q2*(-m4 + 4*Q4 - t2 + 
     &  2*m2*(Q2 + t) - 2*Q2*(t - u) + u2) + 
     &kt6*(-m8 + 2*Q8 - 4*Q6*(t - u) - t*pow(t + u,3) - 
     & 2*Q4*t*(2*t + 3*u) + m6*(-3*Q2 + 4*t + 3*u) + 
     & Q2*(t3 + t2*u + t*u2 + u3) - m4*(2*Q4 - 7*Q2*(t + u) + 
     & 3*(2*t2 + 3*t*u + u2)) + m2*(4*Q6 + 
     &sqr(t + u)*(4*t + u) + 2*Q4*(t + 6*u) - 
     & Q2*(5*t2 + 8*t*u + u2))) - kt4*(3*m10 + 3*m8*(2*Q2 - 
     &4*t - 3*u) + 2*Q8*(t - u) + t*u*pow(t + u,3) - 
     &Q2*t*sqr(t + u)*(2*t + u) + Q6*(t2 + 4*t*u - u2) +
     &2*Q4*(-t3 + t2*u + 3*t*u2 + u3) + m6*(Q4 + 18*t2 + 
     & 26*t*u + 8*u2 - Q2*(16*t + 9*u)) + 
     &m2*(2*Q8 + 2*Q6*(t - 4*u) + Q4*(3*t2 + 2*t*u - 17*u2) + 
     & Q2*u*(3*t2 + 8*t*u - 3*u2) + 
     & sqr(t + u)*(3*t2 - u2)) + m4*(9*Q6 - 12*t3 - 
     & 2*Q4*(t - 7*u) - 24*t2*u - 13*t*u2 - u3 + 
     & Q2*(12*t2 + 11*t*u + 6*u2))) + sqr(-m2 + t + u)*
     & (-(Q4*t*u*(Q2 + t + u)) + m6*sqr(Q2 + t + u) + 
     & m2*(Q6*t + Q4*(t2 + u2) + sqr(t2 + t*u + u2) + 
     & Q2*u*(t2 + t*u + 2*u2)) + m4*(Q4*(t - 2*u) - 
     & Q2*(t2 + 3*t*u + 4*u2) - 2*(t3 + 2*t2*u + 2*t*u2 + 
     & u3))) + kt2*(-2*m12 - 2*Q8*t*u - 
     & Q4*u*pow(t + u,3) + 2*Q2*t*pow(t + u,4) + 
     & t2*pow(t + u,4) + m10*(8*t + 6*u) - 
     & Q6*(t3 + 5*t2*u + 5*t*u2 + u3) - m8*(Q4 + 11*t2 + 
     & 15*t*u + 6*u2 + Q2*(t + 3*u)) + 
     & m6*(3*Q6 + 4*t3 + 5*t2*u + t*u2 + 2*u3 + 
     & Q4*(3*t + 8*u) + Q2*(3*t2 + 8*t*u + 11*u2)) + 
     & m4*(Q6*(t - 7*u) - Q4*(3*t2 + 3*t*u + 14*u2) + t*(4*t3 + 
     & 15*t2*u + 22*t*u2 + 11*u3) - 
     & Q2*(t3 + 5*t2*u + 9*t*u2 + 13*u3)) + m2*(2*Q8*t +
     &  Q6*(5*t2 + 4*t*u + 5*u2) - 
     & t*sqr(t + u)*(4*t2 + 7*t*u + 5*u2) + Q4*(3*t3 + 
     & 6*t2*u + 3*t*u2 + 8*u3) - Q2*(3*t4 + 8*t3*u + 
     & 6*t2*u2 - 5*u4)))) + 2*kt2*(-2*m12 + 4*kt8*Q4 + 
     & m10*(7*t + 6*u) + kt6*Q2*(-m4 + 8*Q4 - t2 + 
     & 2*m2*(3*Q2 + t) - 2*Q2*(t - u) + u2) - 2*m8*(Q4 + 
     & 4*t2 + 7*t*u + 3*u2 + 2*Q2*(t + u)) + 
     & Q2*t*sqr(t + u)*(-Q4 + sqr(t + u)) + 2*m6*(t3 + 
     & 3*Q4*u + 3*t2*u + 2*t*u2 + u3 + 
     & Q2*(6*t2 + 8*t*u + 6*u2)) - m4*(Q6*t - 
     & 2*Q4*(t2 - 3*u2) - 2*t*(t3 + 3*t2*u + 5*t*u2 + 3*u3) + 
     & Q2*(11*t3 + 22*t2*u + 19*t*u2 + 12*u3)) + 
     & m2*(2*Q6*t*(t + u) - t*sqr(t + u)*(t2 + 2*t*u + 3*u2) + 
     & Q4*(-2*t2*u + 2*u3) + 2*Q2*(t4 + 3*t3*u + 
     & 3*t2*u2 + 3*t*u3 + 2*u4)) + kt4*(-m8 + 4*Q8 - 
     & 4*Q6*(t - u) + Q2*u*sqr(t + u) - t*pow(t + u,3) 
     & + m6*(-Q2 + 4*t + 3*u) - 2*Q4*(3*t2 + 
     & 4*t*u + u2) + m2*(12*Q6 + 20*Q4*u + sqr(t + u)*(4*t + u) 
     & + Q2*(-t2 - 8*t*u + u2)) + 
     & m4*(2*Q4 + Q2*(2*t + 7*u) - 3*(2*t2 + 3*t*u + u2))) 
     & + kt2*(-3*m10 - 2*Q8*(t - u) + 
     & Q4*(t - u)*sqr(t + u) + Q2*t*pow(t + u,3) - 
     & t*pow(t + u,4) + m8*(-5*Q2 + 11*t + 9*u) + 
     & Q6*(-t2 + u2) + m6*(Q4 - 14*t2 - 23*t*u - 
     & 7*u2 + 5*Q2*(3*t + 2*u)) - m4*(9*Q6 - 6*t3 - 
     & 15*t2*u - 8*t*u2 + u3 + Q4*(3*t + 17*u) + 
     & Q2*(14*t2 + 17*t*u + 11*u2)) + m2*(-2*Q8 + Q6*(-6*t + 
     & 8*u) + sqr(t + u)*(t2 + t*u + 2*u2) + Q4*(-3*t2 - 
     & 2*t*u + 17*u2) + Q2*(3*t3 + 4*t2*u - t*u2 + 
     & 6*u3))))*( - sqrt(pt2*kt2)*cos(phi) ) + 
     & 2*(4*kt8*Q4 + 8*kt6*Q4*(m2 + Q2) - 
     & m2*sqr(-m2 + t + u)*(m6 + m2*(t + u)*(2*Q2 + t + u) 
     & - Q2*(3*t2 + 2*t*u + u2) - 
     & m4*(Q2 + 2*(t + u))) - kt2*m2*(m2 - t - u)*(m6 - 
     & t3 - 2*Q4*(t - 3*u) - 3*t2*u - t*u2 + u3 - 
     & 3*m4*(t + u) + 4*Q2*(t2 + u2) + m2*(-6*Q4 + 3*t2 + 
     & 6*t*u + u2 - 4*Q2*(t + u))) + 
     & 2*kt4*Q2*(m6 + 2*m4*(2*Q2 - t + u) + 
     & m2*(8*Q4 + sqr(t - u) + 8*Q2*u) + 
     & 2*Q2*(Q4 - sqr(t + u))))*sqr( - sqrt(pt2*kt2)*cos(phi) ))
     & /(3.0*kt2*sqrt(m2)*sqr(m2 + 
     & Q2 - t)*sqr(kt2 + m2 - u)*sqr(-m2 + t + u)*
     & sqr(kt2 + Q2 + t + u))

c        write(6,*) '  new valuet = ',valt

c      if(valt.lt.0) then
c        write(6,*) ' valuet < 0 for t,u,m2,kt2,q2,pt2,phi'
c        write(6,*) ' valuet = ',valt
c        write(6,*) ' t,u,m2,kt2,q2,pt2,phi = ',t,u,m2,kt2,q2,pt2,phi
c      endif


      NL = - 32.0*4.0*Pi*alpe*(eC)**2*(4.0*Pi*alps)**2 

 	ValL = NL*Q2*(kt2*(4*kt10*Q2 + sqr(m2 - u)*
     &  sqr(-m2 + t + u)*sqr(Q2 + t + u) + 
     & kt8*(-3*m4 + 8*Q4 - 3*t2 - 4*Q2*(t - u) - 
     & 2*t*u + u2 + 2*m2*(2*Q2 + 3*t + u)) - 
     & 2*kt6*(3*m6 - 2*Q6 - t3 + m4*(3*Q2 - 7*t - 6*u) + 
     & 4*Q4*(t - u) + t*u2 + Q2*(5*t2 + 8*t*u + u2) - 
     & m2*(4*Q4 - 5*t2 - 6*t*u + u2 + 2*Q2*(2*t + 7*u))) - 
     & kt4*(7*m8 + 2*m6*(2*Q2 - 9*t) + 4*Q6*(t - u) - 
     & sqr(t + u)*(5*t2 + 2*t*u - 2*u2) + Q4*(3*t2 + 
     & 10*t*u - u2) + Q2*(-6*t3 + 2*t2*u + 14*t*u2 + 6*u3) + 
     & m4*(19*Q4 + 10*t2 - 12*t*u - 7*u2 + Q2*(-10*t + 26*u)) + 
     & 2*m2*(2*Q6 + 3*t3 + Q4*(t - 9*u) + 12*t2*u + 
     & 16*t*u2 - u3 + 2*Q2*(3*t2 + 2*t*u - 9*u2))) - 
     & 2*kt2*(2*m10 - m8*(Q2 + 5*t + 9*u) + m6*(-3*Q4 + t2 - 
     & 2*Q2*u + 16*t*u + 13*u2) + u*(2*Q6*t + Q2*u*
     & sqr(t + u) - t*pow(t + u,3) + Q4*(3*t2 + 4*t*u + u2)) + 
     & m4*(5*t3 - 20*t*u2 - 7*u3 + Q4*(-2*t + 7*u) + 
     & Q2*(3*t2 - 4*t*u + 8*u2)) - m2*(2*Q6*t + 3*t4 + 6*t3*u - 
     & 2*t2*u2 - 10*t*u3 - u4 + Q4*(3*t2 + 2*t*u + 5*u2) + 
     & Q2*(4*t3 + 4*t2*u - 2*t*u2 + 6*u3)))) + 
     & 8*kt2*(2*kt8*Q2 - kt6*(m4 - 4*Q4 + Q2*(t - u) + 
     & t*(t + u) - m2*(3*Q2 + 2*t + u)) - m2*(m2 - t - 
     & u)*(m6 + t3 + Q2*t*(t - u) + t2*u - 3*t*u2 - u3 - 
     & m4*(2*t + 3*u) + m2*(Q2*t + u*(5*t + 3*u))) + 
     & kt4*(-m6 + 2*m4*(t + 2*u) + 2*Q2*(Q4 - 2*t2 - 3*t*u - 
     & u2 + Q2*(-t + u)) + m2*(6*Q4 - t2 - 4*t*u + 
     & u2 + 2*Q2*(t + 6*u))) - kt2*(m8 + Q6*(t - u) + 
     & Q4*t*(t + u) - Q2*(t - u)*sqr(t + u) - t*pow(t + u,3) + 
     & m6*(Q2 - 3*t + 2*u) + m4*(5*Q4 - 2*Q2*(t - 3*u) + 
     & 2*(t2 - 3*t*u - 2*u2)) + m2*(Q6 + t3 + Q4*(2*t - 5*u) + 
     & 7*t2*u + 11*t*u2 + u3 + 4*Q2*(t2 + t*u - 2*u2))))*
     & ( - sqrt(pt2*kt2)*cos(phi) ) + 4*(4*kt8*Q2 - 
     & m2*sqr(-m2 + t + u)*(m4 - t2 + 2*t*u + u2 - 2*m2*(t + u)) - 
     & kt6*(m4 - 8*Q4 + sqr(t + u) - 
     & 2*m2*(4*Q2 + t + u)) + kt2*(m2 - t - u)*(2*m4*(Q2 + t - u) + 
     & Q4*(t + u) - pow(t + u,3) - 
     & m2*(Q4 + t2 - 10*t*u - 3*u2 + Q2*(-6*t + 2*u))) + 
     & 2*kt4*(m4*(3*Q2 + 4*u) + Q2*(2*Q4 - 3*sqr(t + u)) + 
     & 2*m2*(4*Q4 - 2*t*u + Q2*(t + 5*u))))*
     & sqr( - sqrt(pt2*kt2)*cos(phi) ))/
     & (3.0*kt2*sqrt(m2)*sqr(m2 + Q2 - 
     & t)*sqr(kt2 + m2 - u)*sqr(-m2 + t + u)*sqr(kt2 + Q2 + t + u))
      
      return
      end
      double precision function sqr(a)
      double precision a
      sqr= a * a
      return
      end
      double precision function pow(a,i)
      double precision a
      Integer i
      pow = a**i
      return
      end      

c baranovs code

      SUBROUTINE XSEC(XTOTAL)
C*******************************************************************C
C      Partonic matrix element squared                              C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      Parameter(pi=3.1415926)
      Parameter(GeVmb=0.4D0, GeVmub=4.D2, GeVnb=4.D5, GeVpb=4.D8)
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2
      COMMON/CONST/ALP,ALS, XC,XC2,XJ,XJ2, QC,PSI
      COMMON/LOOP8/AMP8(4,4,9)
      COMMON/LOOPJ/AMPJ(4,4)
      COMMON/SPINJ/SPJ(4,4)
      COMMON/GMUNU/DF(4,4),DC(4)
      SAVE
C...Momentum dot-products
      CALL SCALAR
C...Gluon polarization vectors
      CALL GAUGEG
C...Meson polarization vectors
      CALL GAUGEJ(1) ! g_mu_nu -(p_mu*p_nu)/XJ2
c     CALL GAUGEJ(2) ! Include decay into mu+mu- 
C...Evaluating the Matrix Elements
      CALL FEYNJ
      XTOTAL=0.D0
      COLOR=0.D0    !Color and charge coefficients 
      IF(ITYPE.EQ.0) COLOR= 8./3. *QC    ! glu*glu -> Jpsi*gam
      IF(ITYPE.EQ.1) COLOR=10./9.        ! glu*glu -> Jpsi*glu
      IF(ITYPE.EQ.2) COLOR= 8./3. *QC    ! gam*glu -> Jpsi*glu
      IF(ITYPE.EQ.3) COLOR=   12. *QC**3 ! gam*gam -> Jpsi*gam
      IF(ITYPE.EQ.4) COLOR= 5./6. *QC    ! glu*glu -> Jpsi8*gam
      IF(ITYPE.EQ.5) COLOR= 1.           ! glu*glu -> Jpsi8*glu
      IF(ITYPE.EQ.6) COLOR= 2.    *QC**2 ! gam*glu -> Jpsi8*gam
      IF(ITYPE.EQ.7) COLOR= 5./6. *QC    ! gam*glu -> Jpsi8*glu
      IF(ITYPE.EQ.8) COLOR= 2.    *QC**2 ! gam*gam -> Jpsi8*glu
      FACTOR=0.D0   !Coupling constants
      IF(ITYPE.EQ.0) FACTOR = (4*pi)**3 *PSI *ALP    *ALS**2 /4./64.
      IF(ITYPE.EQ.1) FACTOR = (4*pi)**3 *PSI         *ALS**3 /4./64.
      IF(ITYPE.EQ.2) FACTOR = (4*pi)**3 *PSI *ALP    *ALS**2 /4./8.
      IF(ITYPE.EQ.3) FACTOR = (4*pi)**3 *PSI *ALP**3         /4.
      IF(ITYPE.EQ.4) FACTOR = (4*pi)**3 *PSI *ALP    *ALS**2 /4./64.
      IF(ITYPE.EQ.5) FACTOR = (4*pi)**3 *PSI         *ALS**3 /4./64.
      IF(ITYPE.EQ.6) FACTOR = (4*pi)**3 *PSI *ALP**2 *ALS    /4./8.
      IF(ITYPE.EQ.7) FACTOR = (4*pi)**3 *PSI *ALP    *ALS**2 /4./8.
      IF(ITYPE.EQ.8) FACTOR = (4*pi)**3 *PSI *ALP**2 *ALS    /4.
      DO 12  l=1,4
      DO 11 mj=1,4
      DO 11 nj=1,4
      IF(ITYPE.EQ.5) GOTO 5
      XTOTAL = XTOTAL -AMPJ(mj,l)*AMPJ(nj,l)*SPJ(mj,nj)*DC(l)
      if(xtotal.ne.xtotal) then
         write(6,*) ' xsec: AMPJ(mj,l),AMPJ(nj,l),mj,nj',
     &      AMPJ(mj,l),AMPJ(nj,l),mj,nj
       call pylist(1)
       endif
      IF(ITYPE.NE.5) GOTO 10
   5  ANS0=
     .  AMP8(mj,l,1)*AMP8(nj,l,1)*72.
     . +AMP8(mj,l,2)*AMP8(nj,l,2)*200./9.
     . +AMP8(mj,l,3)*AMP8(nj,l,3)*72.
     . +AMP8(mj,l,4)*AMP8(nj,l,4)*200./9.
     . +AMP8(mj,l,5)*AMP8(nj,l,5)*72.
     . +AMP8(mj,l,6)*AMP8(nj,l,6)*200./9.
     . +AMP8(mj,l,7)*AMP8(nj,l,7)*64.
     . +AMP8(mj,l,8)*AMP8(nj,l,8)*64.
     . +AMP8(mj,l,9)*AMP8(nj,l,9)*64.
     . +AMP8(mj,l,7)*AMP8(nj,l,8)*16.
     . +AMP8(mj,l,7)*AMP8(nj,l,9)*16.
     . +AMP8(mj,l,8)*AMP8(nj,l,9)*16.
      ANS1=
     .  AMP8(mj,l,1)*AMP8(nj,l,3)*72.
     . -AMP8(mj,l,1)*AMP8(nj,l,4)*40.
     . +AMP8(mj,l,1)*AMP8(nj,l,5)*72.
     . +AMP8(mj,l,1)*AMP8(nj,l,6)*40.
     . +AMP8(mj,l,1)*AMP8(nj,l,7)*48.
     . -AMP8(mj,l,1)*AMP8(nj,l,9)*48.
     . -AMP8(mj,l,2)*AMP8(nj,l,3)*40.
     . -AMP8(mj,l,2)*AMP8(nj,l,4)*40./3.
     . +AMP8(mj,l,2)*AMP8(nj,l,5)*40.
     . -AMP8(mj,l,2)*AMP8(nj,l,6)*40./3.
     . +AMP8(mj,l,2)*AMP8(nj,l,7)*80./3.
     . +AMP8(mj,l,2)*AMP8(nj,l,9)*80./3.
      ANS2=
     . -AMP8(mj,l,3)*AMP8(nj,l,5)*72.
     . +AMP8(mj,l,3)*AMP8(nj,l,6)*40.
     . +AMP8(mj,l,3)*AMP8(nj,l,7)*48.
     . -AMP8(mj,l,3)*AMP8(nj,l,8)*48.
     . -AMP8(mj,l,4)*AMP8(nj,l,5)*40.
     . -AMP8(mj,l,4)*AMP8(nj,l,6)*40./3.
     . +AMP8(mj,l,4)*AMP8(nj,l,7)*80./3.
     . +AMP8(mj,l,4)*AMP8(nj,l,8)*80./3.
     . +AMP8(mj,l,5)*AMP8(nj,l,8)*48.
     . -AMP8(mj,l,5)*AMP8(nj,l,9)*48.
     . +AMP8(mj,l,6)*AMP8(nj,l,8)*80./3.
     . +AMP8(mj,l,6)*AMP8(nj,l,9)*80./3.
      XTOTAL = XTOTAL - (ANS0 +ANS1 +ANS2)*SPJ(mj,nj)*DC(l)
      if(xtotal.ne.xtotal) then
         write(6,*) ' xsec 5',XTOTAL,ANS0,ANS1,ANS2,SPJ(mj,nj),DC(l)
       endif
  10  CONTINUE
  11  CONTINUE
  12  CONTINUE
      XTOTAL = XTOTAL*COLOR*FACTOR*16.*XC2 !*GeVnb
      if(xtotal.ne.xtotal) then
         write(6,*)' contd ',COLOR,FACTOR,XC2
      endif
      RETURN
      END

      SUBROUTINE SCALAR
C*******************************************************************C
C      Dot-products of particle momenta                             C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      COMMON/MOMEN/g1(4),g2(4),pj(4),g3(4)
      COMMON/DOTPR/DG1P,DG2P,DG3P,DG12,DG13,DG23,G1W,G2W
      SAVE
      DG1P=DOTB(g1,pj)
      DG2P=DOTB(g2,pj)
      DG3P=DOTB(g3,pj)
      DG12=DOTB(g1,g2)
      DG13=DOTB(g1,g3)
      DG23=DOTB(g2,g3)
       G1W=DOTB(g1,g1)
       G2W=DOTB(g2,g2)
      RETURN
      END

      SUBROUTINE GAUGEJ(Mode)
C********************************************************************C
C       Spin density of the outgoing J/psi via its Leptonic decay    C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      COMMON/CONST/ALP,ALS, XC,XC2,XJ,XJ2, QC,PSI
      COMMON/MOMEN/g1(4),g2(4),pj(4),g3(4)
      COMMON/GMUNU/DF(4,4),DC(4)
      COMMON/MUONS/v1(4),v2(4)
      COMMON/SPINJ/SPJ(4,4)
      SAVE
      GOTO(20,10),Mode
   10 DO 14 mj=1,4
      DO 14 nj=1,4
      SPJ(mj,nj)=3*(v1(mj)*v2(nj)+v2(mj)*v1(nj)-XJ2*DF(mj,nj)/2.)/XJ2
   14 CONTINUE
      GOTO 99
   20 DO 24 mj=1,4
      DO 24 nj=1,4
      SPJ(mj,nj)= -(DF(mj,nj) - pj(mj)*pj(nj)/XJ2)
   24 CONTINUE
   99 RETURN
      END

      SUBROUTINE GAUGEG
C*******************************************************************C
C      Photon and Gluon polarization vectors                        C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      COMMON/MOMEN/g1(4),g2(4),pj(4),g3(4)
      COMMON/POLAR/            a1(4),a2(4)
      COMMON/KINEM/x1,x2,g1T,g2T, pjT,yj,zj
      COMMON/GMUNU/DF(4,4),DC(4)
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2
      SAVE
      sq=sqrt(2.)
      g1T =dsqrt(g1(2)**2+g1(3)**2)
      g2T =dsqrt(g2(2)**2+g2(3)**2)
      a1(1)=0.D0
      a1(4)=0.D0
      a2(1)=0.D0
      a2(4)=0.D0
      DO 2 i=2,3
      a1(i)=g1(i)*sq/g1T
      a2(i)=g2(i)*sq/g2T
    2 Continue
C
      GOTO 9 !Comment out this line for gauge invariance tests
      DO 4 i=1,4
      a1(i)=g1(i)
      a2(i)=g2(i)
    4 Continue
    9 RETURN
      END

      SUBROUTINE FEYNJ
C*******************************************************************C
C      Matrix Elements for all Feynman diagrams                     C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2
      COMMON/CONST/ALP,ALS, XC,XC2,XJ,XJ2, QC,PSI
      COMMON/MOMEN/g1(4),g2(4),pj(4),g3(4)
      COMMON/POLAR/            a1(4),a2(4)
      COMMON/DOTPR/DG1P,DG2P,DG3P,DG12,DG13,DG23,W1,W2
      COMMON/LOOP8/AMP8(4,4,9)
      COMMON/LOOPJ/AMPJ(4,4)
      COMMON/GMUNU/DF(4,4),DC(4)
      SAVE
      W3=0.D0 !real final gluon
      SHAT=(W1+W2+2*DG12)
      THAT=(W2+W3-2*DG23)
      UHAT=(W1+W3-2*DG13)
C...Fill main array
      DO 1 mj=1,4
      DO 1 m3=1,4
      AMPJ(mj,m3)=0.D0
      DO 1 ic=1,9
      AMP8(mj,m3,ic)=0.D0
   1  CONTINUE
      DO 10 m1=1,4
      DO 10 m2=1,4
      DO 10 m3=1,4
      DO 10 mj=1,4
C...Diagram  =<< 1 >> = << 6 >>
      TERM01=  - DF(m1,m2)*DF(m3,mj)*DG13 - DF(m1,m2)*g1(m3)*g1(mj) +
     + DF(m1,m2)*g1(m3)*g3(mj) - DF(m1,m2)*g1(mj)*g2(m3) + DF(m1,m3)*
     + DF(m2,mj)*DG13 - DF(m1,m3)*g1(m2)*g3(mj) - DF(m1,m3)*g1(mj)*
     + g3(m2) - DF(m1,mj)*DF(m2,m3)*DG13 + DF(m1,mj)*g1(m2)*g1(m3) +
     + DF(m1,mj)*g1(m2)*g2(m3) + DF(m1,mj)*g1(m3)*g3(m2) + DF(m2,m3)*
     + g1(mj)*g3(m1) + DF(m2,m3)*g2(m1)*g3(mj) - DF(m2,m3)*g3(m1)*
     + g3(mj) - DF(m2,mj)*g1(m3)*g2(m1) - DF(m2,mj)*g2(m1)*g2(m3) +
     + DF(m2,mj)*g2(m3)*g3(m1) + DF(m3,mj)*g1(m2)*g3(m1) - DF(m3,mj)*
     + g2(m1)*g3(m2) + DF(m3,mj)*g3(m1)*g3(m2)
      TERM01=TERM01/((W1-DG1P)*(W3+DG3P))*a1(m1)*a2(m2)*DC(m1)*DC(m2)
C...Diagram  << 2 >> = << 5 >>
      TERM02= - DF(m1,m2)*DF(m3,mj)*DG12 + DF(m1,m2)*g1(m3)*g2(mj) +
     + DF(m1,m2)*g1(mj)*g2(m3) + DF(m1,m3)*DF(m2,mj)*DG12 - DF(m1,m3)*
     + g1(m2)*g1(mj) - DF(m1,m3)*g1(m2)*g2(mj) + DF(m1,m3)*g1(mj)*
     + g3(m2) + DF(m1,mj)*DF(m2,m3)*DG12 + DF(m1,mj)*g1(m2)*g1(m3) -
     + DF(m1,mj)*g1(m2)*g2(m3) - DF(m1,mj)*g1(m3)*g3(m2) - DF(m2,m3)*
     + g1(mj)*g2(m1) - DF(m2,m3)*g2(m1)*g2(mj) + DF(m2,m3)*g2(mj)*
     + g3(m1) - DF(m2,mj)*g1(m3)*g2(m1) + DF(m2,mj)*g2(m1)*g2(m3) -
     + DF(m2,mj)*g2(m3)*g3(m1) + DF(m3,mj)*g1(m2)*g3(m1) + DF(m3,mj)*
     + g2(m1)*g3(m2) - DF(m3,mj)*g3(m1)*g3(m2)
      TERM02=TERM02/((W1-DG1P)*(W2-DG2P))*a1(m1)*a2(m2)*DC(m1)*DC(m2)
C...Diagram  << 3 >> = << 4 >>
      TERM03= - DF(m1,m2)*DF(m3,mj)*DG23 - DF(m1,m2)*g1(m3)*g2(mj) -
     + DF(m1,m2)*g2(m3)*g2(mj) + DF(m1,m2)*g2(m3)*g3(mj) - DF(m1,m3)*
     + DF(m2,mj)*DG23 + DF(m1,m3)*g1(m2)*g3(mj) + DF(m1,m3)*g2(mj)*
     + g3(m2) - DF(m1,m3)*g3(m2)*g3(mj) + DF(m1,mj)*DF(m2,m3)*DG23 -
     + DF(m1,mj)*g1(m2)*g1(m3) - DF(m1,mj)*g1(m2)*g2(m3) + DF(m1,mj)*
     + g1(m3)*g3(m2) - DF(m2,m3)*g2(m1)*g3(mj) - DF(m2,m3)*g2(mj)*
     + g3(m1) + DF(m2,mj)*g1(m3)*g2(m1) + DF(m2,mj)*g2(m1)*g2(m3) +
     + DF(m2,mj)*g2(m3)*g3(m1) - DF(m3,mj)*g1(m2)*g3(m1) + DF(m3,mj)*
     + g2(m1)*g3(m2) + DF(m3,mj)*g3(m1)*g3(m2)
      TERM03=TERM03/((W2-DG2P)*(W3+DG3P))*a1(m1)*a2(m2)*DC(m1)*DC(m2)
C
      IF(ITYPE.NE.5) GOTO 5 !Color-octet section starts
C...Diagram << 21 >> = - << 22 >>
      TERM21= - DF(m1,m2)*DF(m3,mj)*DG13 + DF(m1,m2)*DF(m3,mj)*DG23 -
     + DF(m1,m2)*g1(m3)*g1(mj) + DF(m1,m2)*g1(m3)*g2(mj) + DF(m1,m2)*
     + g1(m3)*g3(mj) - DF(m1,m2)*g1(mj)*g2(m3) + DF(m1,m2)*g2(m3)*
     + g2(mj) - DF(m1,m2)*g2(m3)*g3(mj) - 2*DF(m1,m3)*g1(m2)*g3(mj) + 2
     + *DF(m1,mj)*g1(m2)*g1(m3) + 2*DF(m1,mj)*g1(m2)*g2(m3) + 2*
     + DF(m2,m3)*g2(m1)*g3(mj) - 2*DF(m2,mj)*g1(m3)*g2(m1) - 2*
     + DF(m2,mj)*g2(m1)*g2(m3) + 2*DF(m3,mj)*g1(m2)*g3(m1) - 2*
     + DF(m3,mj)*g2(m1)*g3(m2)
      TERM21=TERM21/((W3+DG3P)*SHAT)*a1(m1)*a2(m2)*DC(m1)*DC(m2)
C...Diagram << 23 >> = - << 24 >>
      TERM23= - 2*DF(m1,m2)*g1(mj)*g2(m3) - 2*DF(m1,m3)*g1(mj)*g3(m2) -
     + DF(m1,mj)*DF(m2,m3)*DG12 - DF(m1,mj)*DF(m2,m3)*DG13 + 2*
     + DF(m1,mj)*g1(m2)*g2(m3) + 2*DF(m1,mj)*g1(m3)*g3(m2) + DF(m2,m3)*
     + g1(mj)*g2(m1) + DF(m2,m3)*g1(mj)*g3(m1) + DF(m2,m3)*g2(m1)*
     + g2(mj) + DF(m2,m3)*g2(m1)*g3(mj) - DF(m2,m3)*g2(mj)*g3(m1) -
     + DF(m2,m3)*g3(m1)*g3(mj) - 2*DF(m2,mj)*g2(m1)*g2(m3) + 2*
     + DF(m2,mj)*g2(m3)*g3(m1) - 2*DF(m3,mj)*g2(m1)*g3(m2) + 2*
     + DF(m3,mj)*g3(m1)*g3(m2)
      TERM23=TERM23/((W1-DG1P)*THAT)*a1(m1)*a2(m2)*DC(m1)*DC(m2)
C...Diagram << 25 >> = - << 26 >>
      TERM25= - 2*DF(m1,m2)*g1(m3)*g2(mj) - DF(m1,m3)*DF(m2,mj)*DG12 -
     + DF(m1,m3)*DF(m2,mj)*DG23 + DF(m1,m3)*g1(m2)*g1(mj) + DF(m1,m3)*
     + g1(m2)*g2(mj) + DF(m1,m3)*g1(m2)*g3(mj) - DF(m1,m3)*g1(mj)*
     + g3(m2) + DF(m1,m3)*g2(mj)*g3(m2) - DF(m1,m3)*g3(m2)*g3(mj) - 2*
     + DF(m1,mj)*g1(m2)*g1(m3) + 2*DF(m1,mj)*g1(m3)*g3(m2) - 2*
     + DF(m2,m3)*g2(mj)*g3(m1) + 2*DF(m2,mj)*g1(m3)*g2(m1) + 2*
     + DF(m2,mj)*g2(m3)*g3(m1) - 2*DF(m3,mj)*g1(m2)*g3(m1) + 2*
     + DF(m3,mj)*g3(m1)*g3(m2)
      TERM25=TERM25/((W2-DG2P)*UHAT)*a1(m1)*a2(m2)*DC(m1)*DC(m2)
C...Diagram << 41 + 51 >>
      TERM41 =
     +  - DF(m1,m2)*DF(m3,mj)*W1 + DF(m1,m2)*DF(m3,mj)*W2 + 2*DF(m1,m2)
     + *DF(m3,mj)*DG13 - 2*DF(m1,m2)*DF(m3,mj)*DG23 + DF(m1,m2)*g1(m3)*
     + g1(mj) - 3*DF(m1,m2)*g1(m3)*g2(mj) - DF(m1,m2)*g1(m3)*g3(mj) + 3
     + *DF(m1,m2)*g1(mj)*g2(m3) - DF(m1,m2)*g2(m3)*g2(mj) + DF(m1,m2)*
     + g2(m3)*g3(mj) + 2*DF(m1,m3)*g1(m2)*g1(mj) + 2*DF(m1,m3)*g1(m2)*
     + g2(mj) + 2*DF(m1,m3)*g1(m2)*g3(mj) - 4*DF(m1,mj)*g1(m2)*g1(m3)
     +  - 4*DF(m1,mj)*g1(m2)*g2(m3) - 2*DF(m2,m3)*g1(mj)*g2(m1) - 2*
     + DF(m2,m3)*g2(m1)*g2(mj) - 2*DF(m2,m3)*g2(m1)*g3(mj) + 4*
     + DF(m2,mj)*g1(m3)*g2(m1) + 4*DF(m2,mj)*g2(m1)*g2(m3) - 4*
     + DF(m3,mj)*g1(m2)*g3(m1) + 4*DF(m3,mj)*g2(m1)*g3(m2)
      TERM41=(TERM41/SHAT +DF(m1,mj)*DF(m2,m3)-DF(m1,m3)*DF(m2,mj))/XJ2
     .      *a1(m1)*a2(m2)*DC(m1)*DC(m2)
C...Diagram << 42 + 52 >>
      TERM42 =
     + 2*DF(m1,m2)*g1(mj)*g2(m3) - 2*DF(m1,m2)*g2(m3)*g2(mj) + 2*
     + DF(m1,m2)*g2(m3)*g3(mj) + 2*DF(m1,m3)*g1(mj)*g3(m2) - 2*
     + DF(m1,m3)*g2(mj)*g3(m2) + 2*DF(m1,m3)*g3(m2)*g3(mj) + DF(m1,mj)*
     + DF(m2,m3)*W2 - DF(m1,mj)*DF(m2,m3)*W3 + 2*DF(m1,mj)*DF(m2,m3)*
     + DG12 + 2*DF(m1,mj)*DF(m2,m3)*DG13 - 4*DF(m1,mj)*g1(m2)*g2(m3) -
     + 4*DF(m1,mj)*g1(m3)*g3(m2) - DF(m2,m3)*g1(mj)*g2(m1) - DF(m2,m3)*
     + g1(mj)*g3(m1) - DF(m2,m3)*g2(m1)*g2(mj) - 3*DF(m2,m3)*g2(m1)*
     + g3(mj) + 3*DF(m2,m3)*g2(mj)*g3(m1) + DF(m2,m3)*g3(m1)*g3(mj) + 4
     + *DF(m2,mj)*g2(m1)*g2(m3) - 4*DF(m2,mj)*g2(m3)*g3(m1) + 4*
     + DF(m3,mj)*g2(m1)*g3(m2) - 4*DF(m3,mj)*g3(m1)*g3(m2)
      TERM42=(TERM42/THAT +DF(m1,m2)*DF(mj,m3)-DF(m1,m3)*DF(m2,mj))/XJ2
     .      *a1(m1)*a2(m2)*DC(m1)*DC(m2)
C...Diagram << 43 + 53 >>
      TERM43 =
     +  - 2*DF(m1,m2)*g1(m3)*g1(mj) + 2*DF(m1,m2)*g1(m3)*g2(mj) + 2*
     + DF(m1,m2)*g1(m3)*g3(mj) + DF(m1,m3)*DF(m2,mj)*W1 - DF(m1,m3)*
     + DF(m2,mj)*W3 + 2*DF(m1,m3)*DF(m2,mj)*DG12 + 2*DF(m1,m3)*
     + DF(m2,mj)*DG23 - DF(m1,m3)*g1(m2)*g1(mj) - DF(m1,m3)*g1(m2)*
     + g2(mj) - 3*DF(m1,m3)*g1(m2)*g3(mj) + 3*DF(m1,m3)*g1(mj)*g3(m2)
     +  - DF(m1,m3)*g2(mj)*g3(m2) + DF(m1,m3)*g3(m2)*g3(mj) + 4*
     + DF(m1,mj)*g1(m2)*g1(m3) - 4*DF(m1,mj)*g1(m3)*g3(m2) - 2*
     + DF(m2,m3)*g1(mj)*g3(m1) + 2*DF(m2,m3)*g2(mj)*g3(m1) + 2*
     + DF(m2,m3)*g3(m1)*g3(mj) - 4*DF(m2,mj)*g1(m3)*g2(m1) - 4*
     + DF(m2,mj)*g2(m3)*g3(m1) + 4*DF(m3,mj)*g1(m2)*g3(m1) - 4*
     + DF(m3,mj)*g3(m1)*g3(m2)
      TERM43=(TERM43/UHAT +DF(m1,m2)*DF(mj,m3)-DF(m1,mj)*DF(m2,m3))/XJ2
     .      *a1(m1)*a2(m2)*DC(m1)*DC(m2)
    5 CONTINUE !End of Color-octet section
C
      AMPJ(mj,m3)=AMPJ(mj,m3)+(TERM01+TERM02+TERM03)*DC(mj)
      IF(ITYPE.NE.5) GOTO 8
      AMP8(mj,m3,1)=AMP8(mj,m3,1)+
     + ( TERM01/8. +TERM23/4. -TERM42/4.)*DC(mj)
      AMP8(mj,m3,2)=AMP8(mj,m3,2) +TERM01*DC(mj)/8.
      AMP8(mj,m3,3)=AMP8(mj,m3,3) +
     + (-TERM02/8. +TERM25/4. -TERM43/4.)*DC(mj)
      AMP8(mj,m3,4)=AMP8(mj,m3,4) +TERM02*DC(mj)/8.
      AMP8(mj,m3,5)=AMP8(mj,m3,5) +
     + (-TERM03/8. +TERM21/4. -TERM41/4.)*DC(mj)
      AMP8(mj,m3,6)=AMP8(mj,m3,6) +TERM03*DC(mj)/8.
      AMP8(mj,m3,7)=AMP8(mj,m3,7) +TERM03*DC(mj)/12.
      AMP8(mj,m3,8)=AMP8(mj,m3,8) +TERM01*DC(mj)/12.
      AMP8(mj,m3,9)=AMP8(mj,m3,9) +TERM02*DC(mj)/12.
    8 CONTINUE 
   10 CONTINUE
      RETURN
      END

      SUBROUTINE METRIC
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      COMMON/GMUNU/DF(4,4),DC(4)
      SAVE
      DO 1 I=1,4
      DO 1 J=1,4
      DF(I,J)= 0.D0
    1 CONTINUE !P(i) =(Beam, P_t, P_t, Energy)
      DO 2 I=1,3
      DF(I,I)=-1.D0
    2 DC(I)  =-1.D0
      DF(4,4)= 1.D0
      DC(4)  = 1.D0
      RETURN
      END

      DOUBLE PRECISION FUNCTION DOTB(xx,yy)
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      DIMENSION xx(4),yy(4)   !P(i) = (Beam, P_t, P_t, Energy)
      DOTB= -xx(1)*yy(1) -xx(2)*yy(2) -xx(3)*yy(3) +xx(4)*yy(4)
      RETURN
      END


