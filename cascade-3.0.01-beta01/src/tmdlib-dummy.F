      Subroutine TMDinit(iset)
      Integer Iset,icount
      Data icount/0/
      icount = icount + 1
      if(icount.le.10) write(6,*) ' CASCADE: dummy TMDinit called '
      return
      end
      
      Subroutine TMDpdf(kf,x,xbar,kt,t,upv,usea,dnv,dsea,sea,charm,bottom,glu)
      Integer icount
      Data icount/0/
      icount = icount + 1
      if(icount.le.10) write(6,*) ' CASCADE: dummy TMDpdf called '
      return
      end
      
