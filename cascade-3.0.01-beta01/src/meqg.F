      SUBROUTINE meqg(WT1)
      IMPLICIT None
C
C   q g* ----> q
C
C         P1-----*
C                  *
C                    *-------- q
C                  *
C         P2-----*
C
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
      Double Precision XPQ1(-6:6),XPQ2(-6:6),SFSCAL1,SFSCAL2
	Double Precision STMAX1,STMAX2
	common/myscale/STMAX1,STMAX2
      DOUBLE PRECISION DSH(4),DSHN(4)
      Double Precision DScale
      Double Precision SH_NEW,M2,xh,sig,PT2Q
	Double Precision K1K2,K12,K22
      Double Precision scale
      Double Precision Sfcut
      DATA Sfcut /1.d-3/ 
c      DATA Sfcut /1.d-9/ 
	Double Precision WT1,SH,qg2,qg1,wtg,sum2
      Double Precision alph_em,pcm_new,alphas,alpha_s,alpha_em
      Double Precision yph,yph_bar,xsi_hard(2),gf
      Double Precision dot,dot1,dcasrn
      Double Precision rap1,rap2,rapp,ptstm,zz1,zz2,XPQC2
      Integer npfin
      Integer IQBAR
	Integer KIA1,KIA2,KF1,KF2,Kin1,Kin2
	Integer I
      Integer ikincut,Ipgg,ns_sel
      Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Double Precision E1,E2
	Double Precision CUT
      Double Precision QFT
      DATA CUT /1.d-6/ ! as smallx
c      DATA CUT /1.d-2/ ! 
	
      DOUBLE PRECISION SIN2W,XMW2
      COMMON/ELWEAK/SIN2W,XMW2
	Double Precision QQg,lambda
      Integer NFL
      Data IQBAR/1/
c      Data IQBAR/2/

      Logical ltest
c      data ltest/.true./
      data ltest/.false./
      Logical first
      data first/.true./
      ALPH_EM = ALPH
c     call pylist(1)
      NPFIN = N
c	call pylist(1)
      SH = DOT1(NF1,NF1) + 2.D0 * DOT1(NF1,NF2) * DOT1(NF2,NF2)
      WT1 = 0.d0
      IF(SH.LE.0.) RETURN
c boost  cms
c      write(6,*) ' before boost'
c      call pylist(1)
      CALL PYROBO(0,N,0.D0,0.D0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
c      call pylist(1)
      KIA1=N+1
      KIA2=N+2
      KF1=N+3
      KF2=N+4
      Kin1=N+5
      Kin2=N+6
      N=N+6
            
c select incoming flavor
      IF(Iabs(KE).eq.11) then
        K(NIA2,2) = K(NF1,2)       
      else
        K(NIA1,2) = K(NF1,2)
      endif
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


c      K1K2 = DOT1(KIA1,KIA2)
c      K12 = DOT1(KIA1,KIA1)
c      K22 = DOT1(KIA2,KIA2)
      k12 = P(KIA1,1)**2 + P(KIA1,2)**2
      k22 = P(KIA2,1)**2 + P(KIA2,2)**2
	K1K2 = P(KIA1,1)*P(KIA2,1) + P(KIA1,2)*P(KIA2,2)
      M2 = AM(1)**2
      E1 = K12/M2
      E2 = K22/M2
	sig=0
c	write(6,*) k12,e1,k22,e2
	  
 	sig = 0.5


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
         write(6,*) ' meqg. no valid scale selected ',IQBAR
      endif

      stmax1 = stmax1 *scalfaf
      stmax2 = stmax2 *scalfaf
		
	if(stmax1.ne.stmax1) then
c	   write(6,*) ' meqg: stmax problem ',xg,xsi_hard,SSS
	   stmax1=1d20
	endif
      IF(FIRST) THEN
         IF(IQBAR.EQ.1) THEN
            write(6,*) ' meqg - scale: qq-pair '
         ELSE
            write(6,*) ' meqg - scale: q or q_bar '
         ENDIF
         FIRST=.FALSE.
      ENDIF
c            write(6,*) 'meqg xsi ',xsi_hard,' stmax ',stmax
c            write(6,*) ' new xg,xg_bar ',xg,xg_bar
c            write(6,*) ' new yph,yph_bar ',yph,yph_bar
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
c this is for testing
      ELSEIF(IQ2.EQ.7) THEN
         Scale = max(stmax1,stmax2)**2
      ELSE
         WRITE(6,*) ' NO VALID Q2 SCALE. STOP'
         STOP
      ENDIF
	
      Dscale = scale
      alpha_s = ALPHAS(sqrt(dscale))
c	QQg = Qg0
c	Lambda = 0.2	  
c	alpha_s=12.*3.1415927/(33.-8)
c     &	   /Log(dmax1(Dscale,QQg**2)/Lambda**2)

c      write(6,*) ' new alphas ',alpha_s

      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)
     
      IF(K(NIA1,2).eq.22.and.K(NIA2,2).NE.21) then
         if(stmax2.ge.sqrt(sss)) then
            XPQ2(0) = 0.d0
         else
            call cauniglu(K(2,2),XG2SF,SFSCAL2,STMAX2,XPQ2)
         endif
c	   write(6,*) ' meqg: stmax,xpq1,xpq2 ',
c     &	   STMAX,XPQ1(0),XPQ2(0)
         xpqc2 = 0
         NFL=3
         do i=-NFL,NFL
         XPQC2= xpqc2 + XPQ2(I)*DFLOAT(PYCHGE(I))**2/9.D0
c         write(6,*) ' xpqc2 ',i,xpqc2
         end do
c select the flavour         
         QFT = - dcasrn()*XPQC2
         KPA=-NFL-1
   90    KPA=KPA+1
         QFT = QFT + DBLE(PYCHGE(KPA))**2/9.D0*DBLE(XPQ2(KPA))
         IF(QFT.LT.0.0D0) GOTO 90
         if(abs(kpa).gt.3) then
           write(6,*) ' fatal error in meqg: KPA = ',KPA
         endif
         K(NIA2,2) = KPA
         K(NF1,2) = KPA
         WTG= dble(XPQC2)
         Q2 = -Dot1(nia1,nia1)
c         Sig = 2.D0 * ALPH_EM*ALPH_EM*PI/YY/Q2**2
c         Sig = sig*(1.D0+(1.D0-YY)**2)
c use photon flux in p_semih
         Sig = 4D0* ALPH_EM*PI*PI/Q2
c         write(6,*) ' yy,q2 ',yy,q2
c         write(6,*) sig,wtg,kpa
      else IF(K(NIA1,2).NE.21.AND.K(NIA2,2).EQ.21) THEN
         if(stmax1.ge.sqrt(sss)) then
            XPQ1(0) = 0.d0
         elseif(stmax2.ge.sqrt(sss)) then
            XPQ2(0) = 0.d0
         else
            call cauniglu(K(1,2),XG1SF,SFSCAL1,STMAX1,XPQ1)
            call cauniglu(K(2,2),XG2SF,SFSCAL2,STMAX2,XPQ2)
         endif
c	   write(6,*) ' meqg: stmax,xpq1,xpq2 ',
c     &	   STMAX,XPQ1(0),XPQ2(0)
         XPQC2=XPQ1(1)+XPQ1(2)+XPQ1(3)+XPQ1(-1)+XPQ1(-2)+XPQ1(-3)
         WTG= dble(XPQ2(0)*XPQC2)
      ELSE
	   write(6,*) ' no valid process for ipro = 4 ',
     &   K(NIA1,2),K(NIA2,2)
         WTG= 0.d0
      ENDIF
      
      
      
      SUM2 = Sig 
c	write(6,*) ' meqg :',sig,wtg
C new cut for on-shell DGLAP
	If(Iccfm.eq.0) then
	   If(ires(1).eq.1.and.(Qg1.gt.scale*scalfa)) SUM2 = 0 
	   If(ires(1).eq.1.and.(Qg1.gt.STMAX1**2)) SUM2 = 0 
	   If(ires(2).eq.1.and.(Qg2.gt.scale*scalfa)) SUM2 = 0 
	   If(ires(2).eq.1.and.(Qg2.gt.STMAX2**2)) SUM2 = 0 
	Endif

c      Write(6,*) ' meqg Wtg ',Wtg,XG1SF,XG2SF
      IF(SUM2.LT.0D0) THEN
c     write(6,*) ' SUM2 ',SUM2
         if(sum2.ne.sum2) then
            write(6,*) ' meqg SUM2,xsi_hard ',SUM2,xsi_hard
            write(6,*) ' meqg k12,k22 ',k12,k22
         endif
         SUM2 = 0D0
      ENDIF
 
      SUM2 = SUM2*WTG
c      write(6,*) ' SUM2 ',SUM2

      WT1 = SUM2
      if(sum2.ne.sum2) then

         write(6,*) ' meqg xpq(0),alpha_s',xpq1(1),xpq2(0),alpha_s
         write(6,*) ' meqg SUM2,xsi_hard ',SUM2,xsi_hard
         write(6,*) ' meqg k12,k22 ',k12,k22
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
