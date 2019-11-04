      FUNCTION MULTSCATT(WEIGHT)

      IMPLICIT none
#include "caslujets.inc"
#include "casparas.inc"
#include "caspart.inc"
#include "casluco.inc"
#include "casmult.inc"
#include "caspar.inc"

	Integer Len
      Double Precision X(20),XX(20)
      DOUBLE PRECISION WEIGL,WEIGC,WEIGB,WEIGT
      COMMON /weiglh/weigl,weigc,weigb,weigt
      Integer IPROO
      Double Precision WEIG,MULTSCATT,MSCATT,XSEC,WEIGHT,Rn,ratio
	Integer NCALL,NMAX,I,J
	Double Precision GEV2NB
	Double precision test,xsec_max(10,10)
	Double Precision xlim(10)
c	Data xlim/0.d0,1d-7,1d-6,1d-5,1d-4,1d-3,1d-2,1d-1,0.5d0,1.d0/
c	Data xlim/1.d-12,1d-7,1d-6,1d-5,1d-4,1d-3,1d-2,1d-1,0.5d0,1.d0/
	Data xlim/0.00,0.01,0.05,0.1,0.2,0.4,0.6,0.8,0.99d0,1.d0/
	Double Precision tiny
	data tiny/1d-16/
	Double Precision  pyvar_31,pyvar_32
	logical first
	Integer Ierr
	Integer j1,j2,jj,ii
      Double Precision dcasrn
      EXTERNAL dcasrn
	data first/.true./
	data test/0./
      DATA NCALL/0/
	DATA NMAX/5000/
      DATA GEV2NB/.3893857E+6/
      SAVE xsec_max  
      
      WEIG=0.D0
      MULTSCATT=0.D0
	ncall = 0
	ierr = 0
      Len=8
	if(first) then 
c get maximun weight first....	
	   first=.false.
	   pyvar_31 = pyvar(31)
	   pyvar_32 = pyvar(32)
	   do j1=1,10
	   do j2=1,10
	   xsec_max(j1,j2) = 0
	   enddo
	   enddo
         pyvar(31) = 0
         pyvar(32) = 0
	   Do I =1,20
	   XX(I)=0
	   Enddo 
         do I = 1,500
            CALL dcasrnv(X,LEN)
		do ii=1,len
		  xx(ii) = x(ii)
		enddo
	      do j1=1,10
		  xx(1) = xlim(j1)+tiny
	   	  do j2= 1,10
		    xx(4) = xlim(j2)+tiny
c		    write(6,*) x(1),xx(1),x(4),xx(4)
	          test=MSCATT(XX,WEIGHT,ierr)
c		    if(test.eq.0)  
c     &	    write(6,*) ' multscatt test = ',test,(xx(i),i= 1,len)
	          if(test.ge.xsec_max(j1,j2)) xsec_max(j1,j2)=test
		  enddo
		enddo
         enddo
         pyvar(31) = pyvar_31
         pyvar(32) = pyvar_32
	   do j1=1,10
	   do j2=1,10
	     if(xsec_max(j1,j2).le.0.) then
	       write(6,*) '  multscatt: fatal xsec_max = ',
     &       xsec_max(j1,j2),j1,j2
	       xsec_max(j1,j2)=1d-10
	      endif
	   enddo
	   enddo
	endif
400   Xsec = 0.D0
      ncall = ncall + 1
	Rn = dcasrn()
	Do I =1,20
	  XX(i)=0
	Enddo 
      CALL dcasrnv(X,LEN)
      do ii=1,len
         xx(ii) = x(ii)
	enddo
	xsec=MSCATT(XX,WEIGHT,ierr)
	If(Ierr.ne.0) then
	  Weig = 0.
	  goto 500
	Endif
	Weig = Xsec
c     check values x(1),x(4)
      j1 = 0
	j2 = 0
      jj=1
101   if(jj.gt.9) then
         write(6,*) ' multscatt fatal jj>9 x(1)',x(1),Xlim(jj)
	   stop
	endif	
	if(x(1).gt.xlim(jj).and.x(1).le.xlim(jj+1)) then
	  j1 = jj
	  else 
	  jj = jj + 1
	  goto 101
	endif	
      jj=1
102   if(jj.gt.9) then
         write(6,*) ' multscatt fatal jj>9 x(4)',x(4),Xlim(jj)
	   stop
	endif	
	if(x(4).gt.xlim(jj).and.x(4).le.xlim(jj+1)) then
	  j2 = jj
	  else 
	  jj = jj + 1
	  goto 102
	endif	
	if(j1.le.0.or.j2.le.0) then
	  write(6,*) ' fatal j1,j2 ',j1,j2
	  stop
	endif
      ratio = xsec/xsec_max(j1,j2)
	if(ncall.gt.nmax) then
	  Weig = 0
c	  write(6,*) ' multscatt: too many trials.... ',Ncall,Nmax
	  goto 500
	endif
	if(ratio.gt.1d0) then
	   write(6,*) ' multscatt: ratio > 1 .... ',ratio,' j1,j2 ',j1,j2
	   write(6,*) ' multscatt: xsec_max = xsec ',xsec,xsec_max(j1,j2)
	   xsec_max(j1,j2) = xsec
	   goto 400
	elseif(ratio.lt.0d0) then  
	   write(6,*) ' multscatt: ratio < 0... FATAL !!!! '
	   stop
	elseif(ratio.ne.ratio) then  
	   write(6,*) ' multscatt: ratio NAN ',xsec,xsec_max(j1,j2),j1,j2
	   weig = 0
	   ratio = 0
	   goto 500
	else
	endif
c	write(6,*) ' multscatt  Rn=',Rn,ratio,ncall	
	if(Rn.ge.ratio) goto 400
c	write(6,*) ' multscatt event accepted. xsec=',xsec,xsec_max(j1,j2)	
c	write(6,*) ' multscatt event accepted. Rn=',Rn,ratio,ncall	

500   continue
      
      MULTSCATT=WEIG
	
      RETURN
      END
