       Subroutine casbran(ifla,iflb,x,kx,ky,
     &   tmax,xsi_max,xp,qx,qy,z,tnew,s,Iflag)
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
	Double Precision asbmy
	External asbmy
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

      Double Precision Nonsud
      External Nonsud
      
      Double Precision qt2temp,q0temp,Q0z,zmax_start
      
      Integer ifla,iflb

      double precision qthard
      common/myhard/qthard
            
      Double Precision phip,phipp
	Isel = Iglu
c	Isel = 6

      
      IF(First) then
         Lerr(60) = 0
         First=.FALSE.
      Endif
      
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
      alp_kt=asbmy(dsqrt(kt2_alp))
      alp_max= Capb/log(q0*q0/QCDLAM/QCDLAM)

      phis = atan(ky/kx)
      pp_old=x*sqrt(xsi_max*s)
      if(pp_old.gt.sqrt(s)) then
         pp_old=sqrt(s)
      endif
c      write(6,*)' casbran ifla,iflb ',ifla,iflb
ccccccccccccccccccccccccccccccccccccccccccccccc
cc here do the hacking for KMR
      If(Isel.eq.6) Then 
cccc this condition allows only for one emission      
         If(skt2.le.Q0**2.or.icheck.eq.0) Goto 30
	Endif
	icheck = 0
ccccccccccccccccccccccccccccccccccccccccccccccc
c new stopping condition
      if(pp_old.le.pp_a0) then
c      write(6,*) ' stopping cond reached at pp = ',pp_old,' x = ',sx,
c     &   ' kt =',sqrt(kt2)
         goto 30
      endif
c      write(6,*) ' casbran: after 2 goto 30 '
ccccccccccc	
      if(kf.eq.22)  goto 30
ccccccccccc
	if(kf.eq.2212) then 
         CALL CCFMMAX(sx,a_max)
         CALL CCFMMIN(sx,skt2,pp_old,a_min)
	   elseif(kf.eq.22) then
         CALL CCFMMAXG(sx,a_max)
         CALL CCFMMING(sx,skt2,pp_old,a_min)
	   endif

c new stopping condition
	dx=sx
	dkt2=skt2
	dpp_old = pp_old
	dpp_a0 = pp_a0
      
c      write(6,*) ' casbran ',skt2,sx
      CALL cauniglu(KF,dx,dkt2,dpp_old,xfa_old)
      CALL cauniglu(KF,dx,dkt2,dpp_a0,xfa_a0)

      T0 = 2.D0*LOG(Q0/QCDLAM)

      T0p = T0 + 2.d0*log(x)

      T0 = MAX(T0,0.0001d0)
c bug fix:E. Rodrigues
      T0p = MAX(T0p,0.0001d0)
* E. Rodrigues - debugging

      tmin=t0
      tmin=t0p
      t = tmax
* E. Rodrigues - debugging
c      write(*,*) 'tmin,tmax = ',tmin,tmax
      delta_s = log(t0p)*(t-tmin)-t*(log(t)-1.d0)+tmin*(log(tmin)-1.d0)
      delta_s=exp(0.5d0*Capb*delta_s)
      if(delta_s.gt.1.) then
         write(6,*) ' delta_s ',delta_s,t0,tmin,t
      endif
      a0=xfa_a0(IFLA)
      a0 = a0*delta_s
      if(a0.gt.xfa_old(IFLA)) then
c         write(6,*) ' casbran: a0 check ',a0,xfa_old(IFLA),ifla
         goto 30
      endif


      if(a_max.le.0.0) then
         write(6,*) ' stopping amax < 0 ',a_max
	   LST(21) = 59
         goto 30
      endif
      Bfactor = a_max/a_min
c      write(6,*) ' amax , amin = ',a_max,a_min,bfactor
   20 Ntry = Ntry + 1
      qt2max = QCDLAM**2*EXP(tm_bran)
      if(tm_bran.le.-100.) goto 30
      if(qt2max.le.0.0) then
         write(6,*) ' casbran: qt2max ',qt2max,' tmax ',tm_bran
         goto 30
      endif
      If(Iccfm.eq.0) then
         if(qt2max.le.Q0**2) goto 30
      endif

      zmax = 1.D0 - x*Q0/SQRT(qt2max)
ccc      Q0z =0.2 ! use this as new default
      Q0z =Q0  ! only for testing
      zmax = 1.D0 - Q0z/SQRT(qt2max)
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

      T0 = 2.D0*LOG(Q0/QCDLAM)

      T0 = MAX(T0,0.0001d0)

      tm_bran = min(tmax,tm_bran)


c      write(6,*) ' zmin,zmax ',zmin,zmax
      IF(zmax.le.zmin) then
         Goto 30
      endif
c do z with approximation P_simple = asb_max/(1-z) + asb(kt)/z
      g0 = alp_kt*log(zmax/zmin)*Bfac
      IF(IFLA.NE.0.and.IFLB.ne.0) g0 = 0.
      g1 = alp_max*log((1.d0-zmin)/(1.d0-zmax))*Bfac
      gtot = g0 + g1
      gtot_sud = Bfactor*gtot

c now do the sudakov very simple
      RN = dcasrn()
      RSUD=LOG(RN)
      T0p = T0 + 2.d0*log(x/(1.d0-x))
      if(tm_bran.lt.T0p) then
         ntmcount=ntmcount+1
         if(ntmcount.le.10) then
         write(6,*) ' casbran: tm_bran =',tm_bran,' t0p = ',t0p
         endif
         goto 30
      endif
      T =RSUD/gtot_sud + tm_bran

c


c generate z
c new
      if(dcasrn()*gtot.le.g1) then
c generate 1/(1-z) term
         z = 1.d0-(1.d0-zmin)*(((1.d0-zmax)/(1.d0-zmin))**dcasrn())
      else
c generate 1/z term
         z = zmin*(zmax/zmin)**dcasrn()
      endif
      if(z.lt.0.0.or.z.gt.1.0) then
         write(6,*) 'problem with new z ',z
      endif
      if(z.lt.zmin) then
         write(6,*) ' z < zmin z=',z
         write(6,*) ' zmin = ',zmin,' z_max = ',zmax,' x ',x
         goto 30
      endif
      if(z.gt.zmax) then
         write(6,*) ' z > zmax z=',z
         write(6,*) ' zmin = ',zmin,' z_max = ',zmax,' x ',x
         goto 30
      endif
c translate t into Qt2p scale which is qt'**2
      if(Iccfm.eq.1) Then
         Qt2p=QCDLAM**2*EXP(T) / z**2 
         tnew = dlog(qt2p/qcdlam**2)
         Qt2p=QCDLAM**2*EXP(T)
      else
         Qt2p=QCDLAM**2*EXP(T)
         tnew = t
      endif

      if(qt2p.gt.s) then
c generate new branching
         tm_bran = t
         Goto 20
      Endif
ccccccccccccccccccccccccccccccccccccccccccccccc
c      write(6,*) ' before Isel=6 '

      If(Isel.eq.6) Then 
cc here do the hacking for KMR
         Qt2p = kt2/(1.-z)
      Endif
ccccccccccccccccccccccccccccccccccccccccccccccc
c check on actual zmax cut
c zmax = 1 - Q0/qt_prime
      zm_true = 1.D0 - Q0z/SQRT(Qt2p)
      If(zm_true.le.0.0d0) goto 30
      if(z.gt.zm_true) Then
c generate new branching
         tm_bran = t
         Goto 20
      Endif
c now we have z,qt2p
c here generate phi angle
      phi=2.D0*pi*dcasrn()
      qx=sqrt(Qt2p)*cos(phi)*(1.d0-z)
      qy=sqrt(Qt2p)*sin(phi)*(1.d0-z)

c end of phi angle generation
ccccccccccccccccccccccccccccccccccccccccccccccc
cc here do the hacking for KMR
      If(Isel.eq.6) Then 
cccc check this: this condition removes all branchings      
c	   qx = -kx
c	   qy = -ky
      Endif
ccccccccccccccccccccccccccccccccccccccccccccccc
      qt2 = qx**2 + qy**2
      
      
      yi = (1.d0-z)* x/z
      xsi =qt2/s/yi**2
ccc      write(6,*) ' xsi ',xsi,' xsi_max',xsi_max
      if(xsi.ge.xsi_max) then
c         write(6,*) ' xsi ',xsi,' > xsi_max',xsi_max
         goto 20
      endif
c calculate kt**2 primed
      kt2p = (qx+kx)**2+(qy+ky)**2
      xp=x/z
      sxp=xp


      skt2p = kt2p
cnew....
      if(sqrt(skt2p).le.Q0) then
         tm_bran = t
         goto 20
      endif
	
c      ikincut = 0  ! check this
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
      if(Iccfm.eq.0) pp_prime = sqrt(qt2p)
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
ccccccccccccccccccccccccccccccccccccccccccccccc
cc here do the hacking for KMR
      If(Isel.eq.6) Then
	  dkt2p=max(dkt2p,0.001)
      Endif
ccccccccccccccccccccccccccccccccccccccccccccccc
      CALL cauniglu(KF,dxp,dkt2p,dpp_prime,xfa)
      CALL cauniglu(KF,dx,dkt2,dpp,xfb)

      if(xfb(IFLB).le.0.0) then
c        write(6,*) ' xfb(0) < 0 ',xfb(0)
         tm_bran = t
         goto 20
      endif

      qt2_alp = max(q0temp**2,qt2)
      alp_qt=asbmy(dsqrt(qt2_alp))
      z0 = 1.d0
      qt2d = qt2p
      kt2d = kt2
	
      kt2p_alp=max(q0temp**2,kt2p)
      alp_ktp = asbmy(dsqrt(kt2p_alp))
	
      ktd = dsqrt(max(kt2d,ktcut**2))
      qtd = dsqrt(qt2d)
      ktqt = sqrt(max(kt2d,ktcut**2)/qt2d)
c 
c      write(6,*) ' Casbran: isel = ',Isel    
         IF(IFLA.eq.0.and.IFLB.eq.0) then
            if(Iccfm.eq.1) Then
	         Delta_ns = Nonsud(z)
            else
               Delta_ns = 1.d0
            endif
c here calculate simple splitting function
            Psimple = (alp_max/(1.d0-z) + alp_kt/z) * Bfac
c end of Psimple	
c now calculate full splitting function
            If(Ipgg.eq.0.or.Ipgg.eq.1) then
	         Ptilde = alp_qt*Splitt(1.d0-z) + alp_ktp*Delta_ns*Splitt0(z) 
	      elseif(Ipgg.eq.2.or.ipgg.eq.3) then
	         Ptilde = alp_qt*Splitt(1.d0-z) + alp_qt*Delta_ns*Splitt0(z)
	      else
	         write(6,*) ' wrong splitting fct selected ',Ipgg
	         stop
	      endif
         elseif(IFLA.ne.0.and.IFLB.ne.0) then
            Psimple = (alp_max/(1.d0-z)) 
	      Ptilde = alp_qt*Splitt(1.d0-z)
         else
            write(6,*) ' casbran: IFLA,IFLB =',IFLA,
     &        ' ',IFLB, 'not implemented'
         endif
      
      
      wtz = Ptilde/Psimple
      if(a_max.le.0.0) then
         write(6,*) ' amax < 0 ',a_max
         tm_bran = t
         goto 30
      endif
ccccccccccccccccccccccccccccccccccccccccccccccc
cc here do the hacking for KMR
      If(Isel.eq.6) Then
         sz = z
         wtz = alp_kt*kmr(skt2,sx,sz)
      Endif
ccccccccccccccccccccccccccccccccccccccccccccccc

      wt =  wtz / Bfactor*(xfa(IFLA)/xfb(IFLB))
c

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
ccccccccccccccccccccccccccccccccccccccccccccccc
cc here do the hacking for KMR
      If(Isel.eq.6) Then
cc	   icheck = 1
	Endif
      icheck = 1
ccccccccccccccccccccccccccccccccccccccccccccccc
      return
      End
      function kmr(kt2,sx,sz)
	Implicit None
	Double Precision kmr,kt2,sx,sz
c	write(6,*) ' kmr: kt2,x,z ',kt2,sx,sz
c	kmr = 1.
	kmr = (1.-sz)/sz + sz*(1.-sz) + sz/(1.-sz)
	return
	end
      Function Nonsud(zz)
	Implicit None
#include "casshwr.inc"
	
	Double Precision Nonsud,z,zz
	Double Precision Nonsud1,z1,z0,z0c,z01
	Double Precision Q0	

	Double Precision t0,B,k,q
	Double Precision tq,tk,tc,zc,zt
	Integer ncall

	Double Precision  BB,ktqt,Qcdlamp,scal,qtd,ktd
      Common /splitting/ BB,ktqt,Qcdlamp,scal,qtd,ktd
      Integer IGLU
	Common/CAGLUON/Iglu
	Integer  ikincut,Ipgg,ns_sel
	Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel

	Double Precision t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
	Double Precision t11,t12,t13,t14,t15,t16,t17,t18,t19
	Double Precision t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
	Double Precision t30,t31,t32,t33,t34,t35,t36,t37,t38,t39
	Double Precision t40,t42,t41,t43,t45,t46,t47,t48
	Double Precision t50,t51,t52,t53,t55,t56,t58,t59
	Double Precision t62,t63,t65,t66,t68
	Double Precision t70,t71,t74,t75,t76
	Double Precision t80,t82,t85,t92,t94,t95
	Double Precision t102,t107,t108
	Double Precision a0,b0,c0,d0
	Double Precision one
	Double Precision Ei
	External Ei
	Double Precision asbmy
	External asbmy
c ns_sel = -1: standard non sudakov old
c ns_sel = 0: standard non sudakov with cut
c ns_sel = 1: full splitting
c ns_sel = 2: alphas(q) sing. terms
c ns_sel = 3: alphas(q) full splitting
      Integer Insq
c     Insq = 0 non-sudakov with alphasq with hannes treatment
c	Insq = 1 non-sudakov with alphasq with gavin treatment
      Data insq/0/
	data ncall/0/
	Logical first
	Data first/.true./
      Double Precision Nf,Ca,Capbp 
      Nf = 4.D0
      Ca = 3.D0
      Capbp=12.d0*Ca/(11.d0*Ca-2.d0*nf)	
	if(first) then 
	   write(6,*) ' non-sudakov selection : ns_sel = ',ns_sel
	   write(6,*) ' Insq = ',Insq
	   first=.false.
	endif
	if(Iccfm.eq.0) then
	  Nonsud = 1.
	  else 
	
	z1=zz

	if(ktqt.gt.1.) then
	  z0= 1.
	Elseif(z1.lt.ktqt.and.ktqt.le.1.) then
	  z0=ktqt
	elseif(ktqt.le.z1) then
	  z0=z1
	else
	  write(6,*) ' Nonsud : this should never happen '
	Endif
	if(z1.gt.z0) then
	   write(6,*) ' Nonsud z1,z0',z1,z0
	endif
c	  if(QQt.lt.Qtem) write(6,*) ' QQt,Qtem ',QQt,Qtem 
      k=ktd
      q=qtd
      Q0=0.9
      z=z1
      tk = log(k**2/QCDlamp**2)
	tq = log(q**2/QCDlamp**2)
	tc = log(Q0**2/QCDlamp**2)
	zc = Q0/q
	zt = max(z1,zc)
c	zt = z1
c now gavins version
	B=BB
	if(tk.lt.tc) then 
c	  write(6,*) ' Nonsud0: tk<tc '
	endif
	if(ns_sel.eq.-1) then
         Nonsud1 =  log(z0/z1)*log(k**2/z0/z1/q**2)
         Nonsud1 =  -asbmy(scal*k)*Nonsud1
c	   write(6,*) ' nonsudakov old : ',Nonsud1,ns_sel
c         write(6,*) ' nonsudakov old z0,z1,k,q:',z0,z1,k,q
c         write(6,*) ' nonsudakov old logs :',log(z0/z1),
c     &         log(k**2/z0/z1/q**2)
	elseif(ns_sel.eq.0) then
c	   write(6,*) ' nonsudakov old with cut : ',ns_sel
c the simple non - sudakov
	   if(tk.lt.tc) then
	     Nonsud1 = 0.
	   else
           Nonsud1 =  log(z0/zt)*log(k**2/z0/zt/q**2)
	   if(zc.gt.z1) then      
              Nonsud1 = Nonsud1+(log(zc/z1))*log(k**2/Q0**2)
	   endif
         Nonsud1 =  -asbmy(scal*k)*Nonsud1
	   endif  
c	write(6,*) z0,z1,qqt,qtem,asb
	
	elseif(ns_sel.eq.1) then
c	   write(6,*) ' nonsudakov full splitt with cut : ',ns_sel
	   if(tk.lt.tc) then
	     t0 = 0.
	   else           
c
c
c here comes the new part
           t0 = -log(k/z0/q)**2+log(k/z0/q)*z0**2+z0**2/2-
     &          B*log(k/z0/q)*z0**2-
     &          z0**2*B/2-2*z0*log(k/z0/q)-
     &          2*z0-2.D0/3.D0*z0**3*log(k/z0/q)-2.D0/9.D0*z0**3+
     &          2.D0/3.D0*B*z0**3*log(k/z0/q)+2.D0/9.D0*z0**3*B+
     &          log(k/zt/q)**2-log(k/zt/q)*zt**2-
     &          zt**2/2+B*log(k/zt/q)*zt**2+zt**2*B/2+2*zt*log(k/zt/q)+
     &          2*zt+2.D0/3.D0*zt**3*log(k/zt/q)+2.D0/9.D0*zt**3-
     &          2.D0/3.D0*B*zt**3*log(k/zt/q)-2.D0/9.D0*zt**3*B
           
            If(Insq.eq.0) then
	        if(zc.gt.z1) then      
                t1 = -log(k/zc/q)**2+log(k/zc/q)*zc**2+zc**2/2-
     &          B*log(k/zc/q)*zc**2-
     &          zc**2*B/2-2*zc*log(k/zc/q)-
     &          2*zc-2.D0/3.D0*zc**3*log(k/zc/q)-2.D0/9.D0*zc**3+
     &          2.D0/3.D0*B*zc**3*log(k/zc/q)+2.D0/9.D0*zc**3*B+
     &          log(k/z1/q)**2-log(k/z1/q)*z1**2-
     &          z1**2/2+B*log(k/z1/q)*z1**2+z1**2*B/2+2*z1*log(k/z1/q)+
     &          2*z1+2.D0/3.D0*z1**3*log(k/z1/q)+2.D0/9.D0*z1**3-
     &          2.D0/3.D0*B*z1**3*log(k/z1/q)-2.D0/9.D0*z1**3*B
                t0 = t0 + t1
	        endif
     
 	      Elseif(Insq.eq.1) then
c now gavins version ( for full splitting fct )
	          if(zc.gt.z1) then      
                   t1 = 2*log(k/q0)*log(zc)-2*log(k/q0)*zc-
     &	       2.D0/3.D0*log(k/q0)*zc**3+
     &             log(k/q0)*zc**2+2.D0/3.D0*log(k/q0)*B*zc**3-
     &             log(k/q0)*B*zc**2-2*log(k/q0)*log(z1)+
     &             2*log(k/q0)*z1+2.D0/3.D0*log(k/q0)*z1**3-
     &             log(k/q0)*z1**2-2.D0/3.D0*log(k/q0)*B*z1**3+
     &             log(k/q0)*B*z1**2
                   t0 = t0 + t1
c		     write(6,*) ' reached here ',t0,t1
	          endif
	      endif
	   endif
         Nonsud1 =  -asbmy(scal*k)*t0
ccc         Nonsud1 =  -asbmy(2.*k)*t0
c         Nonsud1 =  -asbmy(4.*k)*t0
c         write(6,*) ' nonsud ',asbmy(k),asbmy(2*k),k
	elseif(ns_sel.eq.2) then
c	zt = max(z1,zc)
c	   write(6,*) ' nonsudakov as(q)  with cut : ',ns_sel
c
c	write(6,*) ' nonsud ',asbmy(k),asbmy(2*k),k
c the simple non - sudakov

c here non sudakov with alphas(q) 
	   if(tk.lt.tc) then
	      Nonsud1 = 0.
	   else
            Nonsud1 = 0.5*log(tk/(tq+2.*log(z0)))*tq+
     &          log(tk/(tq+2*log(z0)))*log(z0)+0.5*tq+log(z0) -
     &         (0.5*log(tk/(tq+2.*log(zt)))*tq+
     &          log(tk/(tq+2*log(zt)))*log(zt)+0.5*tq+log(zt))
c            Nonsud1 = log(tk/(tq+2*log(z0)))*tq/2+
c     &	          log(tk/(tq+2*log(z0)))*log(z0)+
c     &                log(z0)-log(tk/(tq+2*log(zt)))*tq/2-
c     &                log(tk/(tq+2*log(zt)))*log(zt)-log(zt)
c            if(Nonsud1.ne.Nonsud2) write(6,*) Nonsud1,Nonsud2
            If(Insq.eq.0) then
               Nonsud1 =  Nonsud1 * Capbp 
	         if(zc.gt.z1) then      
   	            z01=zc

                  if(ktqt.ge.1.) then
	              z0c= 1.
	              z0c=z01
	            Elseif(z01.lt.ktqt.and.ktqt.le.1.) then
	              z0c=min(ktqt,z01)
	            elseif(ktqt.le.z01) then
	              z0c=z01
	            else
	              write(6,*) ' Nonsud : this should never happen '
	            Endif 
	            Nonsud1 = Nonsud1 + asbmy(q0)*
     &	       log(z0c/z1)*log(k**2/z0c/z1/q**2)
	         endif
 	      Elseif(Insq.eq.1) then
	         if(zc.gt.z1) then      
                 Nonsud1 = Nonsud1+(log(zc)-log(z1))*log(tk/tc)
	         endif
               Nonsud1 = Nonsud1 * Capbp  
	      endif
            if(Nonsud1.ne.Nonsud1) then
	        write(6,*) ' problem with new Non Sudakov '
	        write(6,*) ' z0,z1 ',z0,z1,k,q
	        write(6,*) ' log(tk/(tq+2*log(z0)))',
     &	     log(tk/(tq+2.*log(z0)))
	        write(6,*) ' log(tk/(tq+2*log(z1)))',
     &           log(tk/(tq+2.*log(z1)))
	        write(6,*) ' tq ',tq,' tk ',tk,log(z0),log(z1)
	        write(6,*) '  Non - Suadkov: ',Nonsud1
            endif
	   endif
         Nonsud1 = -Nonsud1  
	elseif(ns_sel.eq.3) then
c here non sudakov with fullspiltting and alphas(q) 
	   if(tk.lt.tc) then
	     t0 = 0.
	     Nonsud1 = 0.
	   else
           one = 1.d0
           t1 = log(tk)
           t3 = tq/2
           t4 = log(zt)
           t5 = t3+t4
           t6 = log(t5)
           t7 = zt**2
           t8 = t7*zt
           t9 = t6*t8
           t12 = t6*t7
           t14 = log(z0)
           t16 = 3./2.*tq
           t17 = exp(-t16)
           a0=-3*t4-t16
           t20 = Ei(one,a0)
           t24 = z0**2
           t25 = t24*z0
           t28 = 2*tq
           t31 = log(t28+4*t14)
	     b0=-3*t14-t16
           t36 = Ei(one,b0)
           t39 = exp(-t3)
           t40 = -t14-t3
           t41 = Ei(one,t40)
           t45 = log(-t40)
           t47 = exp(-tq)
	     c0=-2*t14-tq
           t50 = Ei(one,c0)
           t56 = t45*t25
           t58 = -z0*t1-t9/3-t6*zt-t4+t12/2+t1*t14+t14-t17*t20/3+
     &          t6*t4-t25*t1/3
     &          -tq*t31/2+t17*t36/3+t39*t41+t24*t1/2+t45*z0-t47*t50/2+
     &          zt*t1-t7*t1/2+t56/3
           t62 = log(t28+4*t4)
           t65 = t45*t24
           t70 = B*t17
           t80 = B*t47
	     d0=-2*t4-tq
           t92 = Ei(one,d0)
           t95 = Ei(one,-t5)
           t107 = -t1*t4+tq*t62/2-t65/2+t8*t1/3-t45*t14+t70*t20/3-
     &           t24*B*t1/2-
     &           t56*B/3-t12*B/2+t80*t50/2+t25*B*t1/3-t70*t36/3+t9*B/3+
     &           t47*t92/2-t39
     &           *t95+t7*B*t1/2+t65*B/2-t8*B*t1/3-t80*t92/2
           t108 = t58+t107

           t0 = t108
           Nonsud1 = -t0 * Capbp 
	     t0 = 0
           If(Insq.eq.0) then
c              Nonsud1 =  -t0 * Capbp 
	        if(zc.gt.z1) then      
   	           z01=zc

                 if(ktqt.ge.1.) then
	              z0c= 1.
	              z0c=z01
	           Elseif(z01.lt.ktqt.and.ktqt.le.1.) then
	              z0c=min(ktqt,z01)
	           elseif(ktqt.le.z01) then
	              z0c=z01
	           else
	             write(6,*) ' Nonsud : this should never happen '
	        Endif 
           
              t0 = -log(k/z0c/q)**2+log(k/z0c/q)*z0c**2+z0c**2/2-B*
     &	   log(k/z0c/q)*
     &         z0c**2-z0c**2*B/2-2*log(k/z0c/q)*z0c-2*z0c-2.D0/3.D0*
     &         z0c**3*log(k/z0c/q)
     &         -2.D0/9.D0*z0c**3+2.D0/3.D0*B*z0c**3*log(k/z0c/q)+
     &         2.D0/9.D0*z0c**3*B+log(k/z1/q)**2-log(k/z1/q)*z1**2-
     &         z1**2/2+B*log(k/z1/q)*z1**2+z1**2*B/2
     &         +2*log(k/z1/q)*z1+2*z1+2.D0/3.D0*z1**3*log(k/z1/q)+
     &         2.D0/9.D0*z1**3
     &         -2.D0/3.D0*B*z1**3*log(k/z1/q)-2.D0/9.D0*z1**3*B
              t0 = - asbmy(q0)* t0 
              if(t0.gt.0) write(6,*) 'exp(t0) > 0 ',t0
	        Nonsud1 = Nonsud1+ t0
              if(Nonsud1.gt.0) write(6,*) 'exp(Nonsud1) > 0 ',Nonsud1
		  endif
c	        Nonsud1 = Nonsud1/2
c		  write(6,*) ' here ',Nonsud1,t0,insq
 	     Elseif(Insq.eq.1) then
	        if(zc.gt.z1) then      
                 t0 = 2*log(k/q0)*log(zc)-2*log(k/q0)*zc-
     &	     2.D0/3.D0*log(k/q0)*zc**3+
     &           log(k/q0)*zc**2+2.D0/3.D0*log(k/q0)*B*zc**3-
     &           log(k/q0)*B*zc**2-2*log(k/q0)*log(z1)+
     &           2*log(k/q0)*z1+2.D0/3.D0*log(k/q0)*z1**3-
     &           log(k/q0)*z1**2-2.D0/3.D0*log(k/q0)*B*z1**3+
     &           log(k/q0)*B*z1**2
                 t0=-asbmy(scal*k)*t0
                 Nonsud1 = Nonsud1 + t0 
c		  write(6,*) ' here ',Nonsud1,t0,insq
		     
	        endif
	     Endif

           if(Nonsud1.ne.Nonsud1) then
	        write(6,*) ' problem with new Non Sudakov '
	        write(6,*) ' z0,z1 ',z0,z1,k,q
	        write(6,*) '  Non - Suadkov: ',Nonsud1,t0,t108
           endif
	   endif
	else 
	  write(6,*) '  Non - Suadkov: ',ns_sel,' not yet implemented '
        stop
      endif


      Nonsud1 = exp(Nonsud1)
     
     
	If(Nonsud1.lt.0.0.or.Nonsud1.gt.1.0001) Then
	  ncall=ncall+1
	  If(ncall.le.1000) then
	  write(6,*) ' Problem with new Non - Suadkov: ',Nonsud1,ns_sel
         Nonsud1 =  log(z0/z1)*log(k**2/z0/z1/q**2)
         write(6,*) Nonsud1,z0,z1,k,q,ktqt
         write(6,*) 'logs',log(z0/z1),log(k**2/z0/z1/q**2)
         Nonsud1 =  -asbmy(scal*k)*Nonsud1
         write(6,*) Nonsud1,asbmy(scal*k)
         write(6,*) 'exp(Nonsud1)',exp(Nonsud1) 
c	  write(6,*) ' QQt,Qtem ',QQt,Qtem
c        pause
	  endif 
	  Nonsud1 = 1.
	Endif
c	Nonsud1 = 1.
      Nonsud = Nonsud1
	endif
	Return
	End 

      Function Ei(x1,x2)
	Implicit None
	Double Precision x1,x2,Dexpin,Ei
	External Dexpin
	if(x1.eq.1.) then
	  Ei=dexpin(x2)
	  else
	  write(6,*) ' Ei: i= ',x1,' not implemented '
	  Endif
c	  write(6,*) ' Ei = ',ei
	Return
	End 
	Function Splitt(z1)
	Implicit None
	Double Precision Splitt,z1

	Double Precision  BB,ktqt,Qcdlam,scal,qtd,ktd
      Common /splitting/ BB,ktqt,Qcdlam,scal,qtd,ktd
	Integer  ikincut,Ipgg,ns_sel
	Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel

	If(Ipgg.eq.0.or.Ipgg.eq.2) then
c first the simple 1/(1-z)=1/z1	 
c	   write(6,*) ' splitt sing. part Ipgg: ',Ipgg
	  Splitt=1./z1
	  else
c now gavins version
c 	   write(6,*) ' splitt full  Ipgg: ',Ipgg
       Splitt = ((1.-z1)/z1 +BB*(1.-z1)*z1)
	endif
	Return
	End 
	 
	
      Function Splitt0(z)
	Implicit None
	Double Precision Splitt0,z,z1

	Double Precision  BB,ktqt,Qcdlam,scal,qtd,ktd
      Common /splitting/ BB,ktqt,Qcdlam,scal,qtd,ktd
	Integer  ikincut,Ipgg,ns_sel
	Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel

	z1 = z
c first the simple 1/z	 
	If(Ipgg.eq.0.or.Ipgg.eq.2) then
c	   write(6,*) ' splitt0 sing. part Ipgg: ',Ipgg
         Splitt0=1./z1
	   else
c now gavins version
c	   write(6,*) ' splitt0 full  Ipgg: ',Ipgg
         Splitt0 = (1.-z1)/z1 + (1.-BB)*z1*(1.-z1)
	 endif
	 Return
	 End 
	 
	function asbmy(xkt)
      Double Precision TQtp,xkt,Qlam0,Nf,Ca,Capb
	Double Precision asbmy
	Double Precision  BB,ktqt,Qcdlam,scal,qtd,ktd
      Common /splitting/ BB,ktqt,Qcdlam,scal,qtd,ktd
	Integer  ikincut,Ipgg,ns_sel
	Double Precision QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
	integer nff 
	Double Precision mass(6)
	data mass/0.,0.,0.,1.5,5.0,175./
      Integer Ialps
c      data Ialps/1/
      data Ialps/2/
      Integer nfmin,nfmax
      Double Precision ALAM2,Q2EFF,Q2THR,pi,B0,B1,alphas,algq
#include "casludat1.inc"      
      

c      Data Qlam0/0.71/

      Ialps = max(1,MSTU(111))
	
	Qlam0 = Qg0	
c change that to 4 flavors
      Nf = 4.D0
      Ca = 3.D0
      Capb=12.d0*Ca/(11.d0*Ca-2.d0*nf)	
      TQtp=Log((dmax1(xkt,Qlam0))**2/QCDLam**2)
      asbmy = Capb/TQtp

ccc now using the same as in smzvalgg for evolution
      
	pi=4.*atan(1.)
      q2eff = xkt**2     
      nfmin = 3
      nfmax = 5
      nff = 4
      alam2 = qcdlam**2
c alphas from PYTHIA      
C...Evaluate first or second order alpha_strong.
  100 IF(NFf.GT.MAX(2,nfmin)) THEN
        Q2THR=mass(NFf)**2
        IF(Q2EFF.LT.Q2THR) THEN
          NFf=NFf-1
c          ALAM2=ALAM2*(Q2THR/ALAM2)**(2D0/(33D0-2D0*NFf))
c          write(6,*) nff
          if(nff.eq.3) then
             ALAM2=ALAM2*(Q2THR/ALAM2)**(2./27.)
           if(Ialps.eq.2) ALAM2=ALAM2*(log(Q2THR/ALAM2))**(107./2025)
          endif
          GOTO 100
        ENDIF
      ENDIF
  110 IF(NFf.LT.MIN(8,nfmax)) THEN
        Q2THR=mass(NFf+1)**2
        IF(Q2EFF.GT.Q2THR) THEN
          NFf=NFf+1
c          ALAM2=ALAM2*(ALAM2/Q2THR)**(2D0/(33D0-2D0*NF))
          if(nff.eq.5) then 
             ALAM2=ALAM2*(ALAM2/Q2THR)**(2./23.)
             if(Ialps.eq.2) ALAM2=ALAM2*(log(Q2THR/ALAM2)**(-2.*963./13225.))
          endif
          if(nff.eq.6) then 
             ALAM2=ALAM2*(ALAM2/Q2THR)**(2./21.)
             if(Ialps.eq.2) ALAM2=ALAM2*(log(Q2THR/ALAM2)**(-2.*321./3381.)) 
          endif
          GOTO 110
        ENDIF
      ENDIF
	if(nff.gt.6.or.nff.le.0) write(6,*) nff,q2
	
      B0=(33D0-2D0*NFf)/6D0
      ALGQ=LOG(MAX(1.0001D0,Q2EFF/ALAM2))
      IF(Ialps.EQ.1) THEN
        alphas=MIN(10.0,2.*pi/(B0*ALGQ))
      ELSE
        B1=(153D0-19D0*NFf)/6D0
        alphas=MIN(10.0,2.*pi/(B0*ALGQ)*(1.-B1*LOG(ALGQ)/
     &  (B0**2*ALGQ)))
      ENDIF

c      asbmy = alphas
c      write(6,*) 'asbmy ', alphas,xkt
      asbmy = 2.*Ca*alphas/2./pi ! hannes new      
c      write(6,*) 'asbmy ', asbmy,xkt,sqrt(alam2),qcdlam
	return
	end
