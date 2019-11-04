      SUBROUTINE meoffjj(WT1)
      IMPLICIT None
C
C    gluon gluon ----> jpsi jpsi
C
C     A LA SERGEY BARANOV  NEW(31.7.92)
C     new A LA SERGEY BARANOV  NEW(oct.92)
C  ( calculation of A. Lipatov, N.Zotov 
C    new ala S. Baranov
C
C    k_gluon1 = y * e + kt_gamma
C    k_gluon2 = x * p + kt_gluon
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

      Double Precision m2,y,xq2,xkt2
      Double Precision als,alp,psi
c commons for S.baranov
      double Precision XC,XC2,XB,XB2,XM,XM2,CHC,CHB
      COMMON/QUARK1/XC,XC2,XB,XB2,XM,XM2,CHC,CHB
      Double Precision g1,a1,g2,a2,pj1,pj2
      COMMON/MOMEN1/g1(4),a1(4),g2(4),a2(4),pj1(4),pj2(4)
      Integer IB1,IB2
      COMMON/BEAMS1/IB1,IB2
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
      
c      call pylist(1)
      IQBAR = Iqqbar
      eq = DFLOAT(PYCHGE(ihfla))/3.D0
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
        write(6,*) ' meoffjj: wrong flavor selected, program stopped '
        stop
      endif
      IF(INTER.EQ.2) THEN
         write(6,*) ' INTER = 2 not implemented for meoffjj'
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
      if(abs(test).gt.1e-3) then
         write(6,*) ' sh+t+u ,m**2 ',Sh_new+t+u,AM(1)**2+AM(2)**2+K12+
     +   K22,AM(1)
         write(6,*) ' sh+t+u-m**2 ',Sh_new+t+u-AM(1)**2-AM(2)**2-K12-
     +   K22
      endif
c new scale
c check on scale:
      SCALE = PT2Q+4*M2
      IF(IQ2.EQ.1) THEN
         SCALE = 4*M2
      ELSEIF(IQ2.EQ.2) THEN
         SCALE = SH_NEW
      ELSEIF(IQ2.EQ.3) THEN
c         SCALE = 4D0*M2 + PT2Q
         SCALE = 4*M2 + PT2Q
      ELSEIF(IQ2.EQ.4) THEN
         SCALE = Q2
      ELSEIF(IQ2.EQ.5) THEN
         SCALE = Q2 + PT2Q + 4*M2
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
         

cc         M2 = 3.1*3.1 ! check this.... this is only for J/psi
         XM = sqrt(M2)
         XM2 = XM**2
         XC = XM/2
         XC2 = XC**2
         XB = XM/2
         XB2 = XB**2

         CHC = eq**2 
         CHB = eq**2 
         PSI = 0.08/XM ! check this.... this is only for J/psi
ccc         PSI = 0.0876/XJ ! check this.... this is only for J/psi
ccc         PSI = wf/XM

c include a very small kt for the calc of the polarisation
c othersise errors produced         

         IB1 = 8
         IB2 = 8
         
         g1(1) = P(KIA1,3)
         g1(2) = P(KIA1,1)+1.e-10*P(KIA1,4)
         g1(3) = P(KIA1,2)+1.e-10*P(KIA1,4)
         g1(4) = P(KIA1,4)
         g2(1) = P(KIA2,3)+1.e-10*P(KIA2,4)
         g2(2) = P(KIA2,1)+1.e-10*P(KIA2,4)
         g2(3) = P(KIA2,2)
         g2(4) = P(KIA2,4)

         pj1(1) = P(KF1,3)
         pj1(2) = P(KF1,1)
         pj1(3) = P(KF1,2)
         pj1(4) = P(KF1,4)
         pj2(1) = P(KF2,3)
         pj2(2) = P(KF2,1)
         pj2(3) = P(KF2,2)
         pj2(4) = P(KF2,4)

c polarisation vectors  (0,  kt_x, kt_y, 0)/|kt|
         
         A1(1) = 0
         A1(2) = g1(2)/sqrt(g1(2)**2+g1(3)**2)
         A1(3) = g1(3)/sqrt(g1(2)**2+g1(3)**2)
         A1(4) = 0
         A2(1) = 0
         A2(2) = g2(2)/sqrt(g2(2)**2+g2(3)**2)
         A2(3) = g2(3)/sqrt(g2(2)**2+g2(3)**2)
         A2(4) = 0


         CALL METRIC1
         Xtotal = 0.
         Call xsecmm(xtotal)
c         write(6,*) ' meoffjj ',xtotal
c        multiply xtotal with (4*pi)**4 * alpha_s**4 /9./8. * (PSI/XM)**2 
c        with PSI = 0.8/(4.*pi) for J/psi 
         sum2 = xtotal* (4.*pi*alpha_s)**4 /9./8.*psi*psi/(4.*pi)/(4.*pi)
          if(sum2.lt.0) then
c           write(6,*) ' sum2 < 0', sum2
c           sum2=0
         endif
c check only phase space:        
ccc         sum2 = 1 
         
         
         if(sum2.ne.sum2) then
           sum2=-9999.
           write(6,*) ' sum2 NN ...'
         endif
c         write(6,*) xtotal

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
c         stmax1 = yph*sqrt(xsi_hard(1)*SSS)
c         stmax2 = xg*sqrt(xsi_hard(2)*SSS)
         stmax1 = (yph+xg_bar)*sqrt(xsi_hard(1)*SSS)
         stmax2 = (xg+yph_bar)*sqrt(xsi_hard(2)*SSS)
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
      
      
c check with rapidity cut
	if( abs(PYP(KF1,17)).gt.5.or.abs(PYP(KF2,17)).gt.5) sum2=0    

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
      RETURN
   20 WT1 = 0.d0
      CALL PYEDIT(14)
      RETURN
      END
      
