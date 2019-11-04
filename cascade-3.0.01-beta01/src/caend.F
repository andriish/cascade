      SUBROUTINE caend(IK)
      IMPLICIT None
      CHARACTER* 16 CNAM
      CHARACTER* 6 CBEAM1,CBEAM2
      Integer IMIX,IK,LERR
      COMMON /OALPHAS/IMIX
      COMMON/ERR/LERR(1000)

      Integer Nqqb,Nccb,Nbbb,Npr1,Npr2,Npr3
      Common/NEVOUT/Nqqb,Nccb,Nbbb,Npr1,Npr2,Npr3
      Double Precision ACC1,ACC2
      Integer IINT,NCB
      COMMON /INTEGR/ ACC1,ACC2,IINT,NCB
	Integer n1min,n1max,n2min,n2max,n3min,n3max,ncall
	common/caerrstf/n1min,n1max,n2min,n2max,n3min,n3max,ncall

#include "caspart.inc"
#include "caseffic.inc"
#include "casparas.inc"
#include "caslujets.inc"
#include "casluco.inc"
      Double Precision xs_l,xs_c,xs_b
      IF(IK.EQ.10) WRITE(6,*) 'TIME LIMIT REACHED........'
      CALL PYNAME(KBEAM(1,2),CNAM)
      CBEAM1 = CNAM(1:6)
      CALL PYNAME(KBEAM(2,2),CNAM)
      CBEAM2 = CNAM(1:6)
      CALL PYNAME(KINT(2,2),CNAM)

      IF(IPRO.lt.15.and.ipro.gt.0) then
         WRITE(6,*) ' x - section for ',CBEAM1,CBEAM2,' --> Q Q_bar '
     +   //'X"'
      elseif(ipro.eq.18) then
         WRITE(6,*) ' x - section for ',CBEAM1,CBEAM2,' --> QCD  '
     +   //'X"'
      elseif(ipro.eq.500) then
         WRITE(6,*) ' x - section for ',CBEAM1,CBEAM2,' --> Z + X '
      elseif(ipro.eq.503) then
         WRITE(6,*) ' x - section for ',CBEAM1,CBEAM2,' --> Z+qqbar +'
     +   //'X"'
      elseif(ipro.eq.504) then
         WRITE(6,*) ' x - section for ',CBEAM1,CBEAM2,' --> Z+ccbar +'
     +   //'X"'
      elseif(ipro.eq.505) then
         WRITE(6,*) ' x - section for ',CBEAM1,CBEAM2,' --> Z+bbbar +'
     +   //'X"'
      elseif(ipro.eq.510) then
         WRITE(6,*) ' x - section for ',CBEAM1,CBEAM2,' --> W + X'
      elseif(ipro.eq.513) then
         WRITE(6,*) ' x - section for ',CBEAM1,CBEAM2,' --> W+q_iq_j +'
     +   //'X"'
      elseif(ipro.eq.514) then
         WRITE(6,*) ' x - section for ',CBEAM1,CBEAM2,' --> W+ cs  +'
     +   //'X"'
      elseif(ipro.eq.515) then
         WRITE(6,*) ' x - section for ',CBEAM1,CBEAM2,' --> W+bbbar +'
     +   //'X"'
      endif
      if(Ipro.ne.-1) WRITE(6,10000) AVGI,SD
10000 FORMAT('  sigma  = ',G10.3,' nb   +/- ',G10.3)
      IF(IK.EQ.20) THEN
         if(Ipro.eq.-1) then
            write(6,*) ' x - section for LHE process '
         else
            IF(iint.eq.0) CALL SPINFO( 6 )
         endif
         WRITE(6,10100) AVGI,SD
10100 FORMAT('  sigma visible = ',G10.3,' nb   +/- ',G10.3)
         if(Ipro.eq.10) then 
            Write(6,11001) dble(Npr1)/dble(nout)*avgi
11001       Format(' x-section for g*g* -> q qbar ',G9.3,' nb' )           
            Write(6,11002) dble(Npr2)/dble(nout)*avgi
11002       Format(' x-section for gg*  -> gg     ',G9.3,' nb')           
            Write(6,11003) dble(Npr3)/dble(nout)*avgi
11003       Format(' x-section for qg*  -> qg     ',G9.3,' nb')           
         endif
         xs_l = dble(Nqqb)/dble(nout)*avgi
         xs_c = dble(Nccb)/dble(nout)*avgi
         xs_b = dble(Nbbb)/dble(nout)*avgi
         write(6,11004) xs_l,Nqqb
11004    Format(' light  quarks : xsec = ',G9.3,' nb. Nqqb = ',I6 )           
         write(6,11005) xs_c,Nccb
11005    Format(' charm  quarks : xsec = ',G9.3,' nb. Nqqb = ',I6 )           
         write(6,11006) xs_b,Nbbb
11006    Format(' bottom quarks : xsec = ',G9.3,' nb. Nqqb = ',I6 )           

         WRITE(6,*) ' Nr of events generated : ',NIN
         WRITE(6,*) ' Nr of events written: ',NOUT
         WRITE(6,*) ' Error summary on event generation '
         WRITE(6,*) ' Errors and their meaning meaning:'
         if(LERR(200).gt.0) then
            WRITE(6,*) ' CAPS_LHEREAD scalup < kt  : ',LERR(200)
         endif
         if(LERR(201).gt.0) then
            WRITE(6,*) ' CAPS_LHEREAD x>1..........: ',LERR(201)
         endif
            WRITE(6,*) ' CAPS x > 0.999............: ',LERR(45)
         WRITE(6,*) ' CAPS boost PS error.......: ',LERR(46)
         WRITE(6,*) ' CAPS energy not conserved.: ',LERR(100)
         WRITE(6,*) ' CAREMN frag. cuts.........: ',LERR(48),LERR(49)
         WRITE(6,*) ' PYPREP error..............: ',LERR(50)
         WRITE(6,*) ' cascps check 1st..........: ',LERR(55)
         WRITE(6,*) ' cascps no more memory.....: ',LERR(51),LERR(52)
         WRITE(6,*) ' cascps boost error........: ',LERR(53)
         WRITE(6,*) ' cascps xfb(iflb)=0........: ',LERR(54)
         WRITE(6,*) ' cascps xsi < 0 ...........: ',LERR(58)
         WRITE(6,*) ' cascps NTRY > 15000.......: ',LERR(56)
         WRITE(6,*) ' cabran stopping < 0 ......: ',LERR(59)
         WRITE(6,*) ' cabran Wt > 1.............: ',LERR(60)
         WRITE(6,*) ' gg->gg rejection..........: ',LERR(90)
	   
	   write(6,*) ' calls to ccfmstfu:  ',ncall
	   write(6,*) ' ccfmstfu: x out of range (min/max) ',n2min,n2max
	   write(6,*) ' ccfmstfu: k2 out of range (min/max) ',n1min,n1max
	   write(6,*) ' ccfmstfu: p out of range (min/max) ',n3min,n3max

      ENDIF
      RETURN
      END
