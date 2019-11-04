       Subroutine casbran_new(ifla,iflb,x,kx,ky,
     &   tmax,xsi_max,scalup,xp,qx,qy,z,tnew,s,Iflag)
      Implicit None

C Feb 2011: hack for KMR included to allow full parton shower, 
C           not only one emission
C April 2014: introduce Q0z to generate softer parton Q0z=0.2 
C
C              |
C       XB     |
C    (IFLB)    |
C              +-------
C              |
C    XA= XB/Z  |
C     (IFLA)   |
C
      Double precision xg,x0,Qs
      common /cbran/ xg,x0,Qs
	Integer kf
	Common /cbran1/ kf
#include "caspara5.inc"
#include "casshwr.inc"
      Integer IGLU
	Common/CAGLUON/Iglu
	Integer  ikincut,Ipgg,ns_sel
	double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
	Integer Isel
      double precision x,kx,ky,tmax,xp,qx,qy,z,t,xsi,xsi_max,s
      double precision q0,t0,t0p,Rsud,pi,tnew
      double precision Capb,rn,Nf,Ca
      Double Precision kt2,kt2p,qt2p,qt2max,qt2
      Double Precision delta_s,tmin

      Double Precision alp_max,alp_kt,alp_qt,zmin,zmax,zm_true
      Double Precision g0,g1,gtot
      Double Precision phi,yi,z0,Psimple,Ptilde,delta_ns,wtz,wt
      Double Precision a0,kt2_alp,qt2_alp
	Double Precision kt2p_alp,alp_ktp
      Double Precision Bfac,Bfactor,gtot_sud,tm_bran,kt2d
      Double Precision sx,sxp,skt2p,pp,a_max,a_min,pp_a0
      Double Precision skt2,pp_old
	Double Precision xfa_old(-6:6),xfa(-6:6),xfb(-6:6)
	Double Precision xfa_a0(-6:6)
      Double Precision pp_prime,qt2d
      Double Precision dx,dxp,dkt2,dkt2p,dpp_old,dpp_a0,dpp_prime,dpp
      Double Precision dcasrn
      Double Precision ktcut
      Double Precision scalup
      Double Precision ntmcount,ntrycount
      Data ktcut/0.25/
      Integer ntry,ntrym,Iflag
      Double Precision phis
      common/check/ phis,ntry
      Integer Iprint
      Integer Lerr
      COMMON/ERR/LERR(1000)
	Double Precision Splitt,Splitt0
	External Splitt,Splitt0
	Double Precision asb
	External asb
      Double Precision zm_low

c      Character*5 xxxx
c      Data xxxx/'bran'/
      logical first
      Data first/.true./
      Data Iprint/0/
c     Data first/.false./
      Integer icheck 
	Data icheck/1/
      External dcasrn
	Double Precision kmr,sz
	External kmr
	
	Double Precision  BB,ktqt,Qcdlam,scal,qtd,ktd
      Common /splitting/ BB,ktqt,Qcdlam,scal,qtd,ktd

      
      Double Precision qt2temp,q0temp,Q0z,zmax_start
      
      Integer ifla,iflb

      double precision qthard
      common/myhard/qthard
            
      Double Precision phip,phipp
      
      double precision Splitt_int(-6:6), b0
      double precision pdf_weight(-6:6), weight, wrn, wfac, alphas_ratio, wtzo, q
      Integer nflmax,i,iflnew, nflav
      Data nflmax/5/
      
      Double Precision Psplit
      External Psplit
      
	Isel = Iglu
c	Isel = 6

      
      IF(First) then
         Lerr(60) = 0
         First=.FALSE.
      Endif
      do i=-6,6
        Splitt_int(i)=0.
        pdf_weight(i)=0.
      end do
      weight = 0.
      
      
      ifla = iflb
      kt2p = 0.
      qt2p = 0.
      kt2 = 0.
      qx = 0.
      qy = 0.
	
      Q0 = Qg0
      Qs=max(Qs,Q0)
      q0temp = 0.3
      
      pi= 4.D0*datan(1.D0)

      ntry = 0
      ntrym = 15000
      ntrym = 500000
c      ntrym = 1000000
      Iflag = 0
      Bfac = 1.d0
      t0 = 0
      
      zmax_start=max((1.-sqrt(qthard)/SQRT(qt2max)),x)
      zmax_start=-9999.
c      write(6,*) ' casbran: new call '
   10 continue
      if(tmax.le.t0) then
c       write(6,*) ' tmax < 0 ',tmax
      endif
c       write(6,*) ' casbran start: tmax  ',tmax

      tm_bran = tmax
c pp_a0 = max angle for initial distribution
      pp_a0=Q0

c change that to 4 flavors
      Nf = 4.D0
      Ca = 3.D0
      Capb=12.d0*Ca/(11.d0*Ca-2.d0*nf)
      kt2=kx**2+ky**2
      skt2=kt2
      sx=x
      kt2_alp=max(q0**2,kt2)
      alp_kt=asb(dsqrt(kt2_alp))
      alp_max= Capb/log(q0*q0/QCDLAM/QCDLAM)
      alp_max=asb(q0)
c      write(6,*) ' casbran_new: alp_max ',alp_max,Capb/log(q0*q0/QCDLAM/QCDLAM),q0,QCDLAM

      phis = atan(ky/kx)
      pp_old=x*sqrt(xsi_max*s)
      if(pp_old.gt.sqrt(s)) then
         pp_old=sqrt(s)
      endif
	icheck = 0
c      write(6,*)' casbran ifla,iflb ',ifla,iflb
ccccccccccccccccccccccccccccccccccccccccccccccc
c new stopping condition
      if(pp_old.le.pp_a0) then
c      write(6,*) ' stopping cond reached at pp = ',pp_old,' x = ',sx,
c     &   ' kt =',sqrt(kt2)
         goto 30
      endif
c      write(6,*) ' casbran: after 2 goto 30 '
c      write(6,*) ' amax , amin = ',a_max,a_min,bfactor
   20 Ntry = Ntry + 1
      qt2max = QCDLAM**2*EXP(tm_bran)
      if(tm_bran.le.-100.) goto 30
      if(qt2max.le.0.0) then
         write(6,*) ' casbran: qt2max ',qt2max,' tmax ',tm_bran
         goto 30
      endif
      If(Iccfm.eq.0.or.ICCFM.ge.2) then
         if(qt2max.le.Q0**2) goto 30
      endif

c      zmax = 1.D0 - x*Q0/SQRT(qt2max)
c      Q0z =Q0  ! test 
      Q0z =0.01 ! default for PB set 1 and set 2
      zmax = 1.D0 - Q0z/SQRT(qt2max)
c zmax ala Martin Ryskin, Watt eq 13 in arXiv 0909.5529
      if(iglu.ge.410000.and.iglu.le.410010)      zmax = 1. - sqrt(kt2)/(pp+sqrt(kt2))
      zmax = max(x,zmax)
      if(zmax_start.gt.0) then
        zmax = min(zmax,(1.D0 - sqrt(qthard)/SQRT(qt2max)))
        zmax_start=-9999.
      endif

      zmin = x

      If(Ntry.gt.ntrym) then
         LST(21)=56
         ntrycount = ntrycount + 1
         if(ntrycount.le.10) then
           write(6,*) ' casbran ntym reached: at x ',x,' tmax = ',tm_bran,t0
         endif
         goto 30
      endif      
      
c new z generation
      pp=sqrt(QCDLAM**2*EXP(Tm_bran))
      CALL cauniglu(KF,x,kt2,pp,xfb)
c      write(6,*) KF,x,kt2,pp,xfb(0)
      if(xfb(iflb).le.1E-10) then
c         write(6,*) ' xfb small ',xfb(iflb),iflb,x,kt2,pp,TM_bran,QCDlam
         goto 30
      endif 

c IFLA -> IFLB + X
      if(iflb.eq.0) then 
        do i=-nflmax,nflmax
c q -> g q
c          Splitt_int(i) =  16.*(1.-sqrt(x/zmax))/3./sqrt(x)
c g -> g g
          if(i.eq.0) Splitt_int(i) = 6.*log((1.-x)*zmax/x/(1.-zmax))
c        write(6,*) ' integral ',Splitt_int(i),i,x,zmax
        end do
      else
c g-> q qbar
c        Splitt_int(0) = 0.5 * (zmax - x) 
c q -> q g
        Splitt_int(iflb) = 8.* log((1.-x)/(1.-zmax))/3. 
c        write(6,*) ' integral ',Splitt_int(i),i
      endif   
c calculate weight from pdfs 
110   weight = 0 .
      ntry = ntry + 1
      if(ntry.gt.ntrym) goto 30
      do i=-nflmax,nflmax
         pdf_weight(i) = xfb(i)/max(1E-10,xfb(iflb))
c         write(6,*) ' pdf ',pdf_weight(i),xfb(i),Splitt_int(i),i
         weight = weight + Splitt_int(i)*pdf_weight(i)
         if(weight.lt.0) then
c            write(6,*) ' casbran_new: problem with AP weights ',i,Splitt_int(i),pdf_weight(i),weight,xfb(i),xfb(iflb),iflb
         endif
      end do
ccc      weight = max(0.0001d0,weight)
      weight = max(1D-15,weight)
c select flavor
      
      wrn = dcasrn() * weight
      iflnew = -nflmax-1 
150   iflnew = iflnew + 1
c      write(6,*) ' check ',iflnew,wrn,weight,Splitt_int(iflnew)*pdf_weight(iflnew), Splitt_int(iflnew), pdf_weight(iflnew)
      wrn = wrn - Splitt_int(iflnew)*pdf_weight(iflnew)
      if(iflnew.lt.nflmax.and.wrn.gt.0) goto 150

      ifla = iflnew 
c
      nflav = 3
      q=99999.
c generate z and weights
      if(iflb.eq.0.and.ifla.eq.0) then
c g-> g + g
         z = 1./(1.+ ((1.-x)/x)*(((1-zmax)/zmax)*(x/(1.-x)))**dcasrn())            
         wtzo = (1.-z*(1.-z))**2
c         wtz = (6.*(1.-z*(1.-z))**2/(z*(1.-z)))/(6./(z*(1.-z)))
         wtz = Psplit(1,nflav,z,q)/(6./(z*(1.-z)))
c         write(6,*) ' g->gg ',wtz/wtzo,z
c         wtz=1.
c q-> g + q
      elseif(iflb.eq.0) THEN
         z = x/(1.-dcasrn()*(1.-sqrt(x/zmax))**2)
         wtzo=0.5*(1.+(1.-z)**2)*sqrt(z)
c         wtz = (4./3.*(1.+(1.-z)**2)/z)/(8./3./z/sqrt(z))
         wtz = Psplit(2,nflav,z,q)/(8./3./z/sqrt(z))
c         write(6,*) ' q->gq ',wtz/wtzo,z
c g-> q + q_bar
      elseif(ifla.eq.0) then
         z=x*(1.+dcasrn()*(zmax-x))
         wtzo=1.-2.*z*(1.-z)
c         wtz= (0.5*(z**2+(1.-z)**2))/0.5
         wtz= Psplit(4,nflav,z,q)/0.5
c         write(6,*) ' g->qqbar ',wtz/wtzo,z
c q-> q + g
      else
         z = 1.-(1.-x)*((1.-zmax)/(1.-x))**dcasrn()
         wtzo=0.5*(1.+z**2)
c         wtz = (4./3.*(1.+z**2)/(1-z))/(4./3.*2./(1.-z))
         wtz = Psplit(3,nflav,z,q)/(4./3.*2./(1.-z))
c         write(6,*) ' g->qg ',wtz/wtzo,z
      endif 
      if (wtz/wtzo.le.0 ) then
         write(6,*) ' wtz < 0 ',wtz/wtzo,z
      endif
      if(z.le.0.or.z.ge.1.) then
         write(6,*) ' casbran_new FATAL, z out of range ',z,zmax,zmin
      endif
      RN = dcasrn()
      RSUD=LOG(RN)
      T = Rsud/(weight*alp_max*Bfac) +tm_bran ! new t-generation
cc      T = Rsud/(Splitt_int(0)*alp_max) +tm_bran ! new t-generation
      
cccccc here ends the new part      
      
      
      
c translate t into Qt2p scale which is qt'**2
      Qt2p=QCDLAM**2*EXP(T)
      tnew = t

      if(qt2p.gt.s) then
c generate new branching
         tm_bran = t
         Goto 20
      Endif
      
      
      
      
c check on actual zmax cut
c zmax = 1 - Q0/qt_prime
      zm_true = 1.D0 - Q0z/SQRT(Qt2p)
      If(zm_true.le.0.0d0) goto 30
      if(z.gt.zm_true) Then
c generate new branching
         tm_bran = t
         Goto 20
      Endif

c old kt2 calc
c now we have z,qt2p
c here generate phi angle
      phi=2.D0*pi*dcasrn()
      qx=sqrt(Qt2p)*cos(phi)*(1.d0-z)
      qy=sqrt(Qt2p)*sin(phi)*(1.d0-z)
c calculate kt**2 primed
      kt2p = (qx+kx)**2+(qy+ky)**2
c end of phi angle generation
      qt2 = qx**2 + qy**2
      
      
cc this is the new calc
ccc      pt2 = -(z*(1.-z)*Q2p - (1.-z)*Q2)
c      pt2 = (1.-z)**2*Q2p
cc      write(6,*) ' pt2 = ',pt2,' Q2p = ',Q2p,' Q2 = ',Q2,' z = ',z
c      px = sqrt(pt2)*cos(phi)
c      py = sqrt(pt2)*sin(phi)
c      kxp = (- px + kx)/z
c      kyp = (- py + ky)/z
c      kt2p = kxp**2 + kyp**2
c      qx = kxp - kx
c      qy = kyp - ky
c end new calc

      yi = (1.d0-z)* x/z
      xsi =qt2/s/yi**2
ccc      write(6,*) ' xsi ',xsi,' xsi_max',xsi_max
      if(xsi.ge.xsi_max) then
c         write(6,*) ' xsi ',xsi,' > xsi_max',xsi_max
         goto 20
      endif
cc
      if(sqrt(kt2p).gt.scalup) then
c         write(6,*) ' casbran: scalup ',scalup,sqrt(kt2p)
         tm_bran = t
         goto 20
      endif
cc

      xp=x/z
      sxp=xp


      skt2p = kt2p
cnew....
      if(sqrt(skt2p).le.Q0) then
         tm_bran = t
         goto 20
      endif
	
      ikincut = 0  ! check this
c      write(6,*) ' ikincut = ', ikincut
      qt2temp = qt2
      If(ikincut.eq.1) then
c LDC/ Catani kin constraint
         if(kt2.le.z*qt2temp) then
           tm_bran = t
           goto 20
         endif
      Elseif(ikincut.eq.2) then
c martin kwiecinski
         if(kt2.le.z*qt2temp/(1-z)) then
            tm_bran = t
            goto 20
         endif
      Elseif(ikincut.eq.3) then
c BFKL kin constraint
         if(kt2.le.z*kt2p) then
            tm_bran = t
            goto 20
         endif
      Endif
cnew...
      pp_prime = sqrt(qt2p)/z
      if(Iccfm.eq.0.or.ICCFM.ge.2) pp_prime = sqrt(qt2p)
      pp_prime=min(pp_prime,sqrt(s))
      if(pp_prime.gt.sqrt(s+0.01*s)) write(6,*) ' pp" > sqrt(s) ',pp_prime,sqrt(s)
      pp = sqrt(qt2p)         
      if(pp.gt.sqrt(s+0.01*s)) write(6,*) ' pp > sqrt(s) ',pp,sqrt(s)
      dxp=sxp
      dx=sx
      dkt2p=skt2p
      dkt2=skt2
      dpp_prime=pp_prime
      dpp=pp
c	write(6,*) ' pp = ',dpp,' pp" = ',dpp_prime

      CALL cauniglu(KF,dxp,dkt2p,dpp_prime,xfa)
      CALL cauniglu(KF,dx,dkt2,dpp,xfb)

      if(xfb(IFLB).le.0.0) then
c        write(6,*) ' xfb(0) < 0 ',xfb(0)
         tm_bran = t
         goto 20
      endif

      qt2_alp = max(q0temp**2,qt2)
      alp_qt=asb(dsqrt(qt2_alp))
      z0 = 1.d0
      qt2d = qt2p
      kt2d = kt2
	
      kt2p_alp=max(q0temp**2,kt2p)
      alp_ktp = asb(dsqrt(kt2p_alp))
	
      ktd = dsqrt(max(kt2d,ktcut**2))
      qtd = dsqrt(qt2d)
      ktqt = sqrt(max(kt2d,ktcut**2)/qt2d)
c 
c      write(6,*) ' Casbran: isel = ',Isel    
cend old wtz
      alp_qt=asb(max(Q0z,dsqrt(qt2p)))
      alp_kt = asb(max(Q0z,dsqrt(kt2p)))
      alphas_ratio = alp_qt/alp_max
cccc      write(6,*) ' casbran_new ', alphas_ratio,alp_qt,alp_max, Bfac
      wtz = wtz*alphas_ratio
      wt = wtz /bfac *(xfa(IFLA)/xfb(IFLB))

      if(iprint.eq.1) then
         write(6,*) ' wt ',wt,' wtz ',wtz,' xfb/a_max ',
     &	   xfb(IFLB)/a_max
         write(6,*) ' x  ',sx,' kt ',skt2,' p ',pp_old
         write(6,*) ' x"  ',sxp,' kt" ',skt2p,' p" ',pp
      endif
      if(wt.gt.1.) then
c         write(6,*) ' wt > 1 ',wt
         Bfac = Bfac * 2.d0
         Goto 10
      endif
      if(wt.gt.1.) then
         write(6,*) ' wt > 1 ',wt
         write(6,*) ' xfa(0) ',xfa(IFLA),' xfb(0) ',xfb(IFLB)
         write(6,*) ' a_max ',a_max,' a_min',a_min,' Bfactor ',Bfactor
         write(6,*) ' wt ',wt,' wtz ',wtz
     +      ,' xfa/xfb ',xfa(IFLA)/xfb(IFLB)
         write(6,*) ' x  ',sx,' kt ',skt2,' p ',pp
         write(6,*) ' x"  ',sxp,' kt" ',skt2p,' p" ',pp_prime
         write(6,*) ' delta_ns = ',delta_ns,' alp kt/qt',alp_kt,alp_qt
         write(6,*) ' qt2 = ',qt2,qcdlam**2
         write(6,*) ' z ',z,' zmin ',zmin,' zmax_true ',zm_true
         Lerr(60) = Lerr(60) + 1
      endif
      rn = dcasrn()
      if(wt.le.rn) then
         tm_bran = t
         goto 20
      endif



      Iflag = 1
c      write(6,*) ' casbran kt2= ',kt2,' kt2p = ',kt2p,iflag
      Return
   30 continue
c      write(6,*) ' casbran kt2= ',kt2,' kt2p = ',kt2p,iflag
      qt2p=0.d0
      qx=0.d0
      qy=0.d0
      Iflag=0
      icheck = 1
ccccccccccccccccccccccccccccccccccccccccccccccc
      return
      End


      Double Precision Function Psplit(Ip,nflav,z,q)
      Implicit None
include "casparam.inc"
c definition of QCDnum splitting fcts      
      Double Precision FF1SFUNC,dqcP0FFS,dqcP0FFR
      Double Precision GF1SFUNC,FG1SFUNC,dqcP0GFA,dqcP0FGA
      Double Precision GG1SFUNC,dqcP0GGA,dqcP0GGR,dqcP0GGS
      Double Precision PP1SFUNC,PM1SFUNC
      Double Precision p2gga,p2ggb,p2ggc,p2gqa,p2nssa,p2qga,p2psa,
     &                 p2nsb,p2nsma,p2nspa,p2nspc,p2nsmc
      External FF1SFUNC,dqcP0FFS,dqcP0FFR
      External GF1SFUNC,FG1SFUNC,dqcP0FGA
      External GG1SFUNC,dqcP0GGA,dqcP0GGR,dqcP0GGS
      External PP1SFUNC,PM1SFUNC

      
      Double Precision z,q,as1,asb
	External asb
      Integer nflav,Ip,Iordas
      
      Logical First/.True./
      
      Iordas = IordasPS
c      nflav = 3
      as1 = asb(q*(1.-z))	 
      if (first) then 
         call QCDnum_init
         First=.False.
      Endif
c use slitting fcts from QCDnum     
         if(Ip.eq.1) then 
c g -> gg         
            Psplit= dqcP0GGS(z,nflav)*dqcP0GGR(z,nflav) +dqcP0GGA(z,nflav)
c            write(6,*) 'suda ',Psplit,Splittgg(z)*2d0*Ca/z, Psplit/(Splittgg(z)*2d0*Ca/z),z,q
            if(Iordas.ge.2) Psplit = Psplit + as1*GG1SFUNC(z,nflav)
c            if(Iordas.ge.3) Psplit = Psplit + as1*as1*(P2GGA(z,nflav)+P2GGB(z,nflav))/8.d0
         elseif(Ip.eq.2) then 
c q -> gq   (in Petronzio notation:  FG)
	      Psplit= dqcP0GFA(z,nflav)
            if(Iordas.ge.2) Psplit = Psplit + as1*FG1SFUNC(z,nflav) 
c            if(Iordas.ge.3) Psplit = Psplit + as1*as1*P2GQA(z,nflav)/8.d0 
         elseif(Ip.eq.3) then 
c q_i -> q_i g 
            Psplit= dqcP0FFR(z,3)*dqcP0FFS(z,nflav)           
c PQQ (q_i -> q_i + g)
            if(Iordas.ge.2) 
     &       Psplit = Psplit + as1*(PP1SFUNC(z,nflav)+PM1SFUNC(z,nflav))/2.d0 
     &                + as1*(FF1SFUNC(z,nflav)-PP1SFUNC(z,nflav))/2./nflav
c            if(Iordas.ge.3) Psplit = Psplit +
c     &                as1*as1*((P2NSPA(z,nflav)+P2NSMA(z,nflav)+2*P2NSB(z,nflav) 
c     &                         + P2NSPC(z,nflav) )/2./8.+
c     &                         (P2PSA(z,nflav)+P2NSSA(z,nflav))/2./8.d0/nflav)
         elseif(Ip.eq.4) then 
c g -> qq  (in Petronzio notation: GF )
	      Psplit= dqcP0FGA(z,nflav)
	      Psplit= 2.d0*nflav*Psplit
            if(Iordas.ge.2) Psplit = Psplit + as1*GF1SFUNC(z,nflav) 
c            if(Iordas.ge.3) Psplit = Psplit + as1*as1*P2QGA(z,nflav) 
         elseif(Ip.eq.5) then 
c PQQB (q_i -> qbar_i + g)
            Psplit= 0 
c FF1SFUNC (eq11) *0.5 since P^s_{qq} = P^s_{qqbar}                     
           if(Iordas.ge.2) Psplit = Psplit + as1*(PP1SFUNC(z,1)-PM1SFUNC(z,1))/2.d0 
     &                            + as1*(FF1SFUNC(z,nflav)-PP1SFUNC(z,nflav))/2./nflav
c            if(Iordas.ge.3) Psplit = Psplit +
c     &                            as1*as1*((P2NSPA(z,nflav)-P2NSMA(z,nflav)
c     &                            )/2./8.+
c     &                            (P2PSA(z,nflav)-P2NSSA(z,nflav))/2./8.d0/nflav)
         elseif(Ip.eq.6) then 
c PQQ (q_i - > q_j + g (i ne j)
            Psplit= 0.           
            if(Iordas.ge.2) Psplit = Psplit + (nflav-1)*as1*(FF1SFUNC(z,nflav)-PP1SFUNC(z,nflav))/2./nflav 
c            if(Iordas.ge.3) Psplit = Psplit + (nflav-1)*as1*as1*(P2PSA(z,nflav)+P2NSSA(z,nflav))/2./8.d0/nflav
         elseif(Ip.eq.7) then 
c PQQB (q_i - > qbar_j + g (i ne j)
            Psplit= 0.           
            if(Iordas.ge.2) Psplit = Psplit + (nflav-1)*as1*(FF1SFUNC(z,nflav)-PP1SFUNC(z,nflav))/2./nflav 
c            if(Iordas.ge.3) Psplit = Psplit + (nflav-1)*as1*as1*(P2PSA(z,nflav)-P2NSSA(z,nflav))/2./8.d0/nflav
         else
            write(6,*)' Psplit: Ip = ',Ip,' not implemented '
         endif

      Return
      End
c
C the coding of the splitting functions is taken from QCDnum-17-01-11
C Botje, M. (2011).  QCDNUM: fast QCD evolution and convolution, 
C Comput.Phys.Commun., 182(), 490-532, arXiv 1005.1481
c thanks to M. Botje.
c

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0GGS(X,NF)
C     ========================================

      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
      
      IDUM     = NF        !avoid compiler warning      
      dqcP0GGS = 1.D0 / ( 1.D0 - X )
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0GGR(X,NF)
C     ========================================
      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
      


      IDUM     = NF        !avoid compiler warning       
      dqcP0GGR = 6.D0 * X 
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0GGA(X,NF)
C     ========================================

      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x,onemx
      

      IDUM     = NF        !avoid compiler warning       
      ONEMX    = 1.D0 - X
      dqcP0GGA = 6.D0 * ( ONEMX/X + X*ONEMX )
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0GFA(X,NF)
C     ========================================
      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
      
 
      IDUM     = NF        !avoid compiler warning 
      dqcP0GFA = 4.D0 * ( 1.D0 + (1.D0-X)*(1.D0-X) ) / ( 3.D0*X )
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0FGA(X,NF)
C     ========================================
      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
      
 
      IDUM     = NF        !avoid compiler warning 
      dqcP0FGA = 0.5D0 * ( X*X + (1.D0-X)*(1.D0-X) )
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0FFS(X,NF)
C     ========================================

      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
      
 
      IDUM     = NF        !avoid compiler warning 
      dqcP0FFS = 1.D0 / (1.D0-X)
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0FFR(X,NF)
C     ========================================
      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
      

      IDUM     = NF        !avoid compiler warning 
      dqcP0FFR = 4.D0 * ( 1.D0 + X*X ) / 3.D0
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION XF1TFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (10)
      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
      Double Precision FF1TFUNC,FG1TFUNC
      

      XF1TFUNC = X * ( FF1TFUNC(X,NF) + FG1TFUNC(X,NF) )

      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION XG1TFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (10)

      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
      Double Precision GG1TFUNC,GF1TFUNC    

      XG1TFUNC = X * ( GG1TFUNC(X,NF) + GF1TFUNC(X,NF) )

      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION PP1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio NPB175(1980)27

      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
#include "qcdnum.inc"
      
 
      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CPFFX  = (1.+CX2) / C1MX
      CPFFMX = (1.+CX2) / C1PX
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)
 
C--   P_F(x) in FP eq. (4.52)
      AAA = - CPFFX*2.*CLX*CL1MX - (2.*X+3./C1MX)*CLX - .5*C1PX*CLX2
     +      - 5.*C1MX
C--   0.5*P_G(x) in FP eq. (4.53)     
      BBB =   CPFFX*(.5*CLX2+C11S6*CLX+CPIA) + C1PX*CLX + C20S3*C1MX
C--   P_N(x) in FP eq. (4.54)   
      CCC = - CPFFX*C2S3*(C5S3+CLX) - C4S3*C1MX
C--   P_QQ(x) in FP eq. (4.50) 
      PQQ  = C16S9*AAA + 4.*BBB + C2S3*NF*CCC
C--   0.5*P_A(x) in FP eq. (4.55) and P_QQBAR(x) eq. (4.51)      
      PQQB = - C4S9 * ( CPFFMX*CS2X + C1PX*CLX + 2.*C1MX )
 
      PP1SFUNC = PQQ + PQQB
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION PM1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio NPB175(1980)27
      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
#include "qcdnum.inc"
      
 
      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CPFFX  = (1.+CX2) / C1MX
      CPFFMX = (1.+CX2) / C1PX
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

C--   P_F(x) in FP eq. (4.52) 
      AAA = - CPFFX*2.*CLX*CL1MX - (2.*X+3./C1MX)*CLX - .5*C1PX*CLX2
     +      - 5.*C1MX
C--   0.5*P_G(x) in FP eq. (4.53)     
      BBB =   CPFFX*(.5*CLX2+C11S6*CLX+CPIA) + C1PX*CLX + C20S3*C1MX
C--   P_N(x) in FP eq. (4.54)      
      CCC = - CPFFX*C2S3*(C5S3+CLX) - C4S3*C1MX
C--   P_QQ(x) in FP eq. (4.50) 
      PQQ  = C16S9*AAA + 4.*BBB + C2S3*NF*CCC
C--   0.5*P_A(x) in FP eq. (4.55) and P_QQBAR(x) eq. (4.51)      
      PQQB = - C4S9 * ( CPFFMX*CS2X + C1PX*CLX + 2.*C1MX )
 
      PM1SFUNC = PQQ - PQQB
 
      RETURN
      END
      
C     ========================================
      DOUBLE PRECISION FUNCTION PP1TFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio NPB175(1980)27

      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
#include "qcdnum.inc"
      
 
      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CPFFX  = (1.+CX2) / C1MX
      CPFFMX = (1.+CX2) / C1PX
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)
 
C--   P_F(x) in FP eq. (4.52)
      AAA = - CPFFX*2.*CLX*CL1MX - (2.*X+3./C1MX)*CLX - .5*C1PX*CLX2
     +      - 5.*C1MX
C--   0.5*P_G(x) in FP eq. (4.53)     
      BBB =   CPFFX*(.5*CLX2+C11S6*CLX+CPIA) + C1PX*CLX + C20S3*C1MX
C--   P_N(x) in FP eq. (4.54)   
      CCC = - CPFFX*C2S3*(C5S3+CLX) - C4S3*C1MX
C--   DEL(x) in FP eq. (6.40)      
      DDD = 4.D0*CPFFX*CLX*CL1MX + (6.D0/C1MX-5.D0-X)*CLX
     +      + (C1PX-2.D0*CPFFX)*CLX2
C--   P_QQ(x) in FP eq. (4.50) 
      PQQ  = C16S9*AAA + 4.*BBB + C2S3*NF*CCC
C--   0.5*P_A(x) in FP eq. (4.55) and P_QQBAR(x) eq. (4.51)      
      PQQB = - C4S9 * ( CPFFMX*CS2X + C1PX*CLX + 2.*C1MX )
 
      PP1TFUNC = PQQ + PQQB + C16S9*DDD
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION PM1TFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio NPB175(1980)27

      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
#include "qcdnum.inc"
      
      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CPFFX  = (1.+CX2) / C1MX
      CPFFMX = (1.+CX2) / C1PX
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

C--   P_F(x) in FP eq. (4.52) 
      AAA = - CPFFX*2.*CLX*CL1MX - (2.*X+3./C1MX)*CLX - .5*C1PX*CLX2
     +      - 5.*C1MX
C--   0.5*P_G(x) in FP eq. (4.53)     
      BBB =   CPFFX*(.5*CLX2+C11S6*CLX+CPIA) + C1PX*CLX + C20S3*C1MX
C--   P_N(x) in FP eq. (4.54)      
      CCC = - CPFFX*C2S3*(C5S3+CLX) - C4S3*C1MX
C--   DEL(x) in FP eq. (6.40)      
      DDD = 4.D0*CPFFX*CLX*CL1MX + (6.D0/C1MX-5.D0-X)*CLX
     +      + (C1PX-2.D0*CPFFX)*CLX2      
C--   P_QQ(x) in FP eq. (4.50) 
      PQQ  = C16S9*AAA + 4.*BBB + C2S3*NF*CCC
C--   0.5*P_A(x) in FP eq. (4.55) and P_QQBAR(x) eq. (4.51)      
      PQQB = - C4S9 * ( CPFFMX*CS2X + C1PX*CLX + 2.*C1MX )
 
      PM1TFUNC = PQQ - PQQB + C16S9*DDD
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION GG1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (11)
      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
#include "qcdnum.inc"
      

      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CL1MX2 = CL1MX**2
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

      CPGG   = 1./C1MX + 1./X -2. + X - CX2
      CMPGG  = 1./C1PX - 1./X -2. - X - CX2

      AAA   = -16.+ 8.*X+ C20S3*CX2 + C4S3/X + (-6.-10.*X)*CLX +
     +        (-2.)*C1PX*CLX2
      BBB   = 2.* C1MX +  26./9.*(CX2-1./X) - C4S3*C1PX*CLX -
     +        20./9.*CPGG
      CCC   = 27./2.*C1MX + 67./9.*(CX2-1./X)+(-25./3.+11./3.*x-
     +        44./3.*CX2)*CLX+4.*C1PX*CLX2+(67./9.-4.*CLX*CL1MX +
     +        CLX2-CPI2S3)*CPGG + 2.*CMPGG*CS2X

      GG1SFUNC = C2S3*NF*AAA + 1.5*NF*BBB + 9.* CCC

      RETURN
      END
      

C     ========================================
      DOUBLE PRECISION FUNCTION GG1TFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (12)
      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
#include "qcdnum.inc"
      
      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CL1MX2 = CL1MX**2
      CPGFX  = CX2 + C1MX**2
      CPGFMX = CX2 + C1PX**2
      CS1X   = -DMB_DILOG(1.D0-X)
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

      CPGG  = 1./C1MX + 1./X -2. + X - CX2
      CMPGG = 1./C1PX - 1./X -2. - X - CX2

      AAA     = -4.+12.*x-164./9.*CX2+92./9./X+(10.+14.*X+C16S3*CX2+
     +          C16S3/X)*CLX + 2.*C1PX*CLX2
      BBB     = 2.-2.*X+26./9.*(CX2-1./X)-C4S3*C1PX*CLX-
     +          (20./9.+8./3.*CLX)*CPGG
      CCC     = 27./2.*(C1MX)+67./9.*(CX2-1./X)+(11./3.-25./3.*X-
     +          44./3./X)*CLX -4.*(C1PX) * CLX2 + (4.*CLX*CL1MX -
     +          3.*CLX2+22./3.*CLX-CPI2S3+67./9.)*CPGG+
     +          2.*CMPGG*CS2X
      GG1TFUNC = 2./3.*NF*AAA+3./2.*NF*BBB+9.*CCC
      
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION GF1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (11)

*     Warning: 'Furmanski-Petronzio' notation where
*                          
*              |F'| = |FF GF| |F| 
*              |G'| = |FG GG| |G|
*                  
*     In ../src/qcdwfun.f we will swap to math notation
*
*              |F'| = |FF FG| |F| 
*              |G'| = |GF GG| |G|                   

      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
#include "qcdnum.inc"
      
      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CL1MX2 = CL1MX**2
      CPGFX  = CX2 + C1MX**2
      CPGFMX = CX2 + C1PX**2
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

      AAA =   4. - 9.*X + (4.*X-1.)*CLX + (2.*X-1.)*CLX2
     +      + 4.*CL1MX
     +      + (2.*CLX-2.*CLX*CL1MX+CLX2-2.*CL1MX+CL1MX2+CPIE)
     +      * 2. * CPGFX
      DDD =   C182S9 + C14S9*X + C40S9/X + (C136S3*X-C38S3)*CLX
     +      - 4.*CL1MX - (2.+8.*X)*CLX2 + 2.*CS2X*CPGFMX
     +      + (C44S3*CLX-CLX2-2.*CL1MX2+4.*CL1MX+CPIF) * CPGFX

      GF1SFUNC = C2S3*NF*AAA + 1.5*NF*DDD

      RETURN
      END
      
C     ========================================
      DOUBLE PRECISION FUNCTION GF1TFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (12)

*     Warning: 'Furmanski-Petronzio' notation where
*                          
*              |F'| = |FF GF| |F| 
*              |G'| = |FG GG| |G|
*                  
*     In ../src/qcdwfun.f we will swap to math notation
*
*              |F'| = |FF FG| |F| 
*              |G'| = |GF GG| |G|

      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
#include "qcdnum.inc"
      
      
      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CL1MX2 = CL1MX**2
      CPGFX  = CX2 + C1MX**2
      CPGFMX = CX2 + C1PX**2
      CS1X   = -DMB_DILOG(1.D0-X)
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

      CPGG   = 1./C1MX + 1./X -2. + X - CX2
      CMPGG  = 1./C1PX - 1./X -2. - X - CX2

      AAA     = -8./3.-(16./9.+8./3.*CLX+8./3.*CL1MX)*CPGFX
      BBB     = -2.+3.*X+(-7.+8.*X)*CLX-4.*CL1MX + (1.-2.*X)*CLX2
     +          +(-4.*CLX*CL1MX-2.*CLX2-2.*CL1MX+2.*CLX-2.*CL1MX2
     +          +16.*CS1X+ 2.*PI*PI - 10.)*CPGFX
      CCC     = -152./9.+166./9.*X-40./9./X+ (-C4S3-76./3.*X)*CLX+
     +          4.*CL1MX + (2.+8.*X)*CLX2+ (8.*CLX*CL1MX-CLX2-
     +          C4S3*CLX+10./3.*CL1MX+2.*CL1MX2-16.*CS1X-7.*CPI2S3+
     +          178./9.)*CPGFX+2.*CPGFMX*CS2X
      GF1TFUNC = (0.5*NF)**2*AAA+2./3.*NF*BBB+3./2.*NF*CCC
      
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION FG1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (11)

*     Warning: 'Furmanski-Petronzio' notation where
*                          
*              |F'| = |FF GF| |F| 
*              |G'| = |FG GG| |G|
*                  
*     In ../src/qcdwfun.f we will swap to math notation
*
*              |F'| = |FF FG| |F| 
*              |G'| = |GF GG| |G|

      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
#include "qcdnum.inc"      

      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CL1MX2 = CL1MX**2
      CPFGX  = (1.+C1MX**2) / X
      CPFGMX = - (1.+C1PX**2) / X
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

      AAA   = -5./2.- 7./2.*X+(2.+7./2.*X)*CLX+(-1.+0.5*X)*CLX2 
     +        -2.*X*CL1MX+ (-3.*CL1MX-CL1MX2)*CPFGX
      BBB   = 28./9.+65./18.*X+44./9.*CX2+(-12.-5.*X-8./3.*CX2)*CLX+
     +        (4.+X)*CLX2+2.*X*CL1MX+ (-2.*CLX*CL1MX+0.5*CLX2+
     +        11./3.*CL1MX+CL1MX2-0.5*CPI2S3+0.5)*CPFGX+CPFGMX*CS2X
      CCC   = -C4S3*X- (20./9.+C4S3*CL1MX)*CPFGX

      FG1SFUNC = C16S9*AAA+4.*BBB+2./3.*NF*CCC

      RETURN
      END
      
C     ========================================
      DOUBLE PRECISION FUNCTION FG1TFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (12)

*     Warning: 'Furmanski-Petronzio' notation where
*                          
*              |F'| = |FF GF| |F| 
*              |G'| = |FG GG| |G|
*                  
*     In ../src/qcdwfun.f we will swap to math notation
*
*              |F'| = |FF FG| |F| 
*              |G'| = |GF GG| |G|

      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
      
#include "qcdnum.inc"      

      
      idum   = nf      !avoid compiler warning
      
      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CL1MX2 = CL1MX**2
      CPFFX  = (1.+CX2) / C1MX
      CPFFMX = (1.+CX2) / C1PX
      CPFGX  = (1.+C1MX**2) / X
      CPFGMX = - (1.+C1PX**2) / X
      CS1X   = -DMB_DILOG(1.D0-X)
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

      AAA =   -.5 + 4.5*X + (-8.+.5*X)*CLX + 2.*X*CL1MX
     +        + (1.-.5*X)*CLX2
     +        + (CL1MX2+4.*CLX*CL1MX-8.*CS1X-CPIB) * CPFGX
      BBB =     C62S9 - C35S18*X - C44S9*CX2
     +        + (2.+12.*X+C8S3*CX2) * CLX
     +        - 2.*X*CL1MX - (4.+X)*CLX2 + CPFGMX*CS2X
     +        + ( - 2.*CLX*CL1MX - 3.*CLX - 1.5*CLX2
     +        - CL1MX2 + 8.*CS1X + CPIC ) * CPFGX
      FG1TFUNC = C16S9*AAA + 4.*BBB

      
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION FF1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (11)

      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
#include "qcdnum.inc"      
      

      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CPFFX  = (1.+CX2) / C1MX
      CPFFMX = (1.+CX2) / C1PX
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

      AAA = - CPFFX*CLX*(1.5+2.*CL1MX) + 2.*CPFFMX*CS2X
     +      - 1. + X + (.5-1.5*X)*CLX - .5*C1PX*CLX2
      BBB =   CPFFX*(C11S6*CLX+.5*CLX2+CPIA) - CPFFMX*CS2X
     +      + C14S3*C1MX
      CCC = - CPFFX*(C10S9+C2S3*CLX) + C40S9/X - 2.*C1PX*CLX2
     +      - C16S3 + C40S3*X + (10.*X+C16S3*CX2+2.)*CLX
     +      - C112S9*CX2

      FF1SFUNC = C16S9*AAA + 4.*BBB + C2S3*NF*CCC

      RETURN
      END
      
C     ========================================
      DOUBLE PRECISION FUNCTION FF1TFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (12)

      IMPLICIT None
 
      Integer NF,Idum
      Double Precision x
#include "qcdnum.inc"      
      
      
      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CL1MX2 = CL1MX**2
      CPFFX  = (1.+CX2) / C1MX
      CPFFMX = (1.+CX2) / C1PX
      CPFGX  = (1.+C1MX**2) / X
      CPFGMX = - (1.+C1PX**2) / X
      CS1X   = -DMB_DILOG(1.D0-X)
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

      AAA =     CPFFX*(1.5*CLX-2.*CLX2+2.*CLX*CL1MX) + 2.*CPFFMX*CS2X
     +        - 1. + X + (-1.5+.5*X)*CLX + .5*C1PX*CLX2
      BBB =     CPFFX*(C11S6*CLX+.5*CLX2+CPIA) - CPFFMX*CS2X
     +        + C14S3*C1MX
      CCC =   - CPFFX*(C2S3*CLX+C10S9)
     +        - C52S3 + C28S3*X + C112S9*CX2 - C40S9/X
     +        - (10.+18.*X+C16S3*CX2)*CLX + 2.*C1PX*CLX2
      FF1TFUNC = C16S9*AAA + 4.*BBB + C2S3*NF*CCC
      
      RETURN
      END            
        
C     ======================================
      double precision function dmb_dilog(x)
C     ======================================

C--   Cernlib routine DDILOG C332
C--
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../mathlib/gen/c/dilog64.F

      implicit double precision (a-h,o-z)

      DIMENSION C(0:19)

      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (PI3 = PI**2/3, PI6 = PI**2/6, PI12 = PI**2/12)

      DATA C( 0) / 0.42996 69356 08136 97D0/
      DATA C( 1) / 0.40975 98753 30771 05D0/
      DATA C( 2) /-0.01858 84366 50145 92D0/
      DATA C( 3) / 0.00145 75108 40622 68D0/
      DATA C( 4) /-0.00014 30418 44423 40D0/
      DATA C( 5) / 0.00001 58841 55418 80D0/
      DATA C( 6) /-0.00000 19078 49593 87D0/
      DATA C( 7) / 0.00000 02419 51808 54D0/
      DATA C( 8) /-0.00000 00319 33412 74D0/
      DATA C( 9) / 0.00000 00043 45450 63D0/
      DATA C(10) /-0.00000 00006 05784 80D0/
      DATA C(11) / 0.00000 00000 86120 98D0/
      DATA C(12) /-0.00000 00000 12443 32D0/
      DATA C(13) / 0.00000 00000 01822 56D0/
      DATA C(14) /-0.00000 00000 00270 07D0/
      DATA C(15) / 0.00000 00000 00040 42D0/
      DATA C(16) /-0.00000 00000 00006 10D0/
      DATA C(17) / 0.00000 00000 00000 93D0/
      DATA C(18) /-0.00000 00000 00000 14D0/
      DATA C(19) /+0.00000 00000 00000 02D0/

      IF(X .EQ. 1) THEN
       H=PI6
      ELSEIF(X .EQ. -1) THEN
       H=-PI12
      ELSE
       T=-X
       IF(T .LE. -2) THEN
        Y=-1/(1+T)
        S=1
        A=-PI3+HF*(LOG(-T)**2-LOG(1+1/T)**2)
       ELSEIF(T .LT. -1) THEN
        Y=-1-T
        S=-1
        A=LOG(-T)
        A=-PI6+A*(A+LOG(1+1/T))
       ELSE IF(T .LE. -HF) THEN
        Y=-(1+T)/T
        S=1
        A=LOG(-T)
        A=-PI6+A*(-HF*A+LOG(1+T))
       ELSE IF(T .LT. 0) THEN
        Y=-T/(1+T)
        S=-1
        A=HF*LOG(1+T)**2
       ELSE IF(T .LE. 1) THEN
        Y=T
        S=1
        A=0
       ELSE
        Y=1/T
        S=-1
        A=PI6+HF*LOG(T)**2
       ENDIF
       H=Y+Y-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 19,0,-1
       B0=C(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=-(S*(B0-H*B2)+A)
      ENDIF
      DMB_DILOG=H
      RETURN
      END
      
      Subroutine QCDnum_init
      Implicit None
#include "qcdnum.inc"      
      pi     = 3.14159265359
c      c1s3   = 1./3.
      c2s3   = 2./3.
      c4s3   = 4./3.
      c5s3   = 5./3.
      c8s3   = 8./3.
      c14s3  = 14./3.
      c16s3  = 16./3.
      c20s3  = 20./3.
      c28s3  = 28./3.
      c38s3  = 38./3.
      c40s3  = 40./3.
      c44s3  = 44./3.
      c52s3  = 52./3.
      c136s3 = 136./3.
      c11s6  = 11./6.
c      c2s9   = 2./9.
      c4s9   = 4./9.
      c10s9  = 10./9.
      c14s9  = 14./9.
      c16s9  = 16./9.
      c40s9  = 40./9.
      c44s9  = 44./9.
      c62s9  = 62./9.
      c112s9 = 112./9.
      c182s9 = 182./9.
c      c11s12 = 11./12.
      c35s18 = 35./18.
c      c11s3  = 11./3.
c      c22s3  = 22./3.
c      c61s12 = 61./12.
c      c215s1 = 215./12.
c      c29s12 = 29./12.
      cpi2s3 = pi**2/3.
      cpia   = 67./18. - cpi2s3/2.
      cpib   = 4.*cpi2s3
      cpic   = 17./18. + 3.5*cpi2s3
c      cpid   = 367./36. - cpi2s3
      cpie   = 5. - cpi2s3
      cpif   = cpi2s3 - 218./9.
c      cca    = 3.
c      ccf    = (cca*cca-1.)/(2.*cca)
c      ctf    = 0.5
c      catf   = cca*ctf
c      cftf   = ccf*ctf

      Write(6,*) ' QCDnum constants initialised '
c      write(6,*) C20S3, C2S3, C4S3, CPI2S3, C16S9, C11S6, C4S9, C5S3,
c     &                 CPIA, C10S9, C112S9, C16S3, C14S3, C40S3, C40S9, 
c     &                 C136S3, C14S9, C182S9, C38S3, C44S3, C35S18, C44S9, 
c     &                 C62S9, C8S3, C28S3, C52S3,
c     &                 CPIB, CPIC, CPIE, CPIF, PI
      Return
      End
      Double precision function asb(xkt)
      Double Precision TQtp,xkt,Qlam0

	Double Precision xtest
	integer nf 
      
      double precision alphasPDFM,alphasPDF
      external alphasPDFM,alphasPDF
      Double Precision mass(6)
      data mass/0.,0.,0.,1.47,4.5,173./
      
      Integer Ialps
c      data Ialps/1/
      data Ialps/2/
      Integer nfmin,nfmax
      Double Precision ALAM2,Q2EFF,Q2THR,pi,B0,B1,alphas,algq
	Double Precision as
	Double Precision Cf,Ca
      Double precision alfs,asfunc
      External asfunc
      Integer nft,ierrt
	Logical asbfirst
	Data asbfirst/.true./
	Double Precision mz,QCDlam
      Common/casalphas/mz,qcdlam
	Data mz/91.188d0/
      Data QCDlam/0.323 /
c      Ialps = Iordas
	Ca = 3d0
	Cf = 4d0/3d0
c	Qlam0 = Qg	
      Qlam0 = 1.0  ! def
cc     Qlam0 = 1.1 ! variation up	
cc     Qlam0 = 0.9 ! variation up
cc      QCDlam = 0.323! give 0.118 at Mz for 2 loops
      TQtp=Log((dmax1(xkt,Qlam0))**2/QCDLam**2)
      xkt = dmax1(xkt,Qlam0)
c new approach ala shirkov
      xtest=xkt/QCDLam
      if(xtest.lt.0.9999.or.xtest.gt.1.0001) then
	else 
c	  write(6,*) ' asbmy: xtest = ',xtest,xkt,QCDLam
	Endif
c old version
c asb = Ca *alphas/pi
      asb = Capb/TQtp
	pi=4.d0*datan(1.d0)
      q2eff = xkt**2     
      nfmin = 3
      nfmax = 5
      nf = 4
      alam2 = qcdlam**2
c alphas from PYTHIA      
C...Evaluate first or second order alpha_strong.
  100 IF(NF.GT.MAX(2,nfmin)) THEN
        Q2THR=mass(NF)**2
        IF(Q2EFF.LT.Q2THR) THEN
          NF=NF-1
          if(nf.eq.3) then
             ALAM2=ALAM2*(Q2THR/ALAM2)**(2./27.)
           if(Ialps.eq.2) ALAM2=ALAM2*(log(Q2THR/ALAM2))**(107./2025)
          endif
          GOTO 100
        ENDIF
      ENDIF
  110 IF(NF.LT.MIN(8,nfmax)) THEN
        Q2THR=mass(NF+1)**2
        IF(Q2EFF.GT.Q2THR) THEN
          NF=NF+1
          if(nf.eq.5) then 
             ALAM2=ALAM2*(ALAM2/Q2THR)**(2./23.)
             if(Ialps.eq.2) ALAM2=ALAM2*(log(Q2THR/ALAM2)**(-2.*963./13225.))
          endif
          if(nf.eq.6) then 
             ALAM2=ALAM2*(ALAM2/Q2THR)**(2./21.)
             if(Ialps.eq.2) ALAM2=ALAM2*(log(Q2THR/ALAM2)**(-2.*321./3381.)) 
          endif
         GOTO 110
        ENDIF
      ENDIF
	if(nf.gt.6.or.nf.le.0) write(6,*) nf,q2

      B0=(33D0-2D0*NF)/6D0
      ALGQ=LOG(MAX(1.0001D0,Q2EFF/ALAM2))
      IF(Ialps.EQ.1) THEN
        alphas=MIN(10.0,2.*pi/(B0*ALGQ))
      ELSE
        B1=(153D0-19D0*NF)/6D0
        alphas=MIN(10.0,2.*pi/(B0*ALGQ)*(1.-B1*LOG(ALGQ)/
     &  (B0**2*ALGQ)))
      ENDIF
      asb = alphas/pi/2d0  
      if(Ipdf.eq.-1) then 
c         alfs = ASFUNC( xkt*xkt, nft, ierrt )
ccc         asbmy = Ca*alfs/pi ! hannes new
ccc         alfs = 0.1 
         asb = alfs/pi/2d0
         elseif(Ipdf.gt.0) then 
         alfs = alphasPDF(xkt)
cc         alfs = ASFUNC( xkt*xkt, nft, ierrt )
ccc         alfs = 0.1 
         asb = alfs/pi/2d0
       endif
c      write(6,*) ' asb = ',asb   
	return
	end
