	subroutine mcgen(fxnb,mxtry)
      Implicit None
      Double Precision fxnb
      double precision dcasrn
      External fxnb,dcasrn
      Integer mxtry
      Integer Mxdim
      Parameter (Mxdim = 50 )
      Double precision xlb,xub
      Integer ndim,ncall
      Common /mcint1/xlb(mxdim),xub(mxdim),ndim,ncall 
      Double precision fmax
      Common /mcint2/fmax 
      Double precision xrn(mxdim),x(mxdim)
      Double precision ff,ftest,ftrn
      Integer i,icall,itry
      Double Precision XG
      COMMON/XFXNB/XG(20)
      data icall /0/
      icall = icall + 1
      if(icall.le.1) then
      write(6,*) ' ****************************************************'
      write(6,*) ' *                                                  *'
      write(6,*) ' *       Simple MC generation     (MCGEN)           *'
      write(6,*) ' *                                                  *'
c      write(6,*) ' *       Author: H. Jung                            *'
      write(6,*) ' *                                                  *'
c      write(6,100) NCALL,NDIM
c 100  FORMAT('  *        NCALL = ',I6,' NDIM = ',I5, '               *')
      write(6,*) ' ****************************************************'      
      endif
      itry = 0
200   itry = itry + 1
      call dcasrnv(xrn,ndim)
      do i=1,ndim
           x(i) = xlb(i) + xrn(i)*(xub(i)-xlb(I))
           xg(i) = x(i)
      end do
      ff = fxnb(x)
      ftest=ff/fmax
      ftrn = dcasrn(icall)
      if(ftest.lt.ftrn) then
         if(itry.lt.mxtry) then 
           goto 200
           else
           write(6,*) ' MCGEN: mxtry = ',mxtry,' reached '
           endif
         
      endif
      
      
c      write(6,*) ' in MCGEN ',itry,icall
      
      return
      end
