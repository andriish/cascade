	subroutine mcint(fxnb,avgi,sd)
      Implicit None
      double precision fxnb,avgi,sd
      External fxnb
      Integer Mxdim
      Parameter (Mxdim = 50 )
      Double precision xlb,xub
      Integer ndim,ncall
      Common /mcint1/xlb(mxdim),xub(mxdim),ndim,ncall 
      Double precision fmax
      Common /mcint2/fmax 
      double precision ff,fsum,fsum2,sigma2,err
      Double precision xrn(mxdim),x(mxdim)
      Integer i,icall
      
      fmax=-1D20
      
      do icall = 1,ncall 
         call dcasrnv(xrn,ndim)
         do i=1,ndim
           x(i) = xlb(i) + xrn(i)*(xub(i)-xlb(I))
         end do
         ff = fxnb(x)
         if(ff.ge.fmax) fmax = ff
         fsum = fsum+ ff
         fsum2 = fsum2 + ff*ff
      end do
      
      fsum = fsum/ncall
      fsum2 = fsum2/ncall
      sigma2 = fsum2 - fsum*fsum
      err = sqrt(sigma2/ncall)
      write(6,*) ' ****************************************************'
      write(6,*) ' *                                                  *'
      write(6,*) ' *       Simple MC integration    (MCINT)           *'
      write(6,*) ' *                                                  *'
      write(6,*) ' *                                                  *'
      write(6,100) NCALL,NDIM
 100  FORMAT('  *        NCALL = ',I6,' NDIM = ',I5, '               *')
      write(6,*) ' ****************************************************'
      
      write(6,*) ' fsum = ',fsum,' err = ',err
      avgi = fsum
      sd = err
      return
      end
