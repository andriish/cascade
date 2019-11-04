      FUNCTION FXN1(X,WEIGHT)

      IMPLICIT none
      Double Precision X(20)
#include "caslujets.inc"
#include "casparas.inc"
#include "caspart.inc"
#include "casluco.inc"
      DOUBLE PRECISION WEIGL,WEIGC,WEIGB,WEIGT
      COMMON /weiglh/weigl,weigc,weigb,weigt
      DOUBLE PRECISION WeigW1,WeigW2,WeigW3,WeigW4,WeigW5,WeigW6
      COMMON /weigW/WeigW1,WeigW2,WeigW3,WeigW4,WeigW5,WeigW6
      Integer IGENFL,IMEPHI
      COMMON/GENWEI/IGENFL,IMEPHI
      Integer NCALL
      DATA NCALL/0/
      Integer IPROO
      Double Precision WEIG,FXN1,XSEC,WEIGHT
      Integer KWZ,KWZP1,KWZP2
      Common/CASWZ/KWZ,KWZP1,KWZP2
      IMEPHI = 0
      WEIG=0.D0
      FXN1=0.D0
      IPROO = IPRO

      Weig = 0.D0
      Weigl = 0.D0
      Weigc = 0.D0
      Weigb = 0.D0
      Weigt = 0.D0
      Xsec = 0.D0

      if(ipro.eq.2) then
         Call Xsecpsi(X,Xsec)
         Weig = Xsec
      elseif(ipro.eq.3) then
         Call Xsecchi(X,Xsec)
         Weig = Xsec
      elseif(ipro.eq.4) then
         Call Xsecqg(X,Xsec)
         Weig = Xsec
      elseif(ipro.eq.21) then
         Call Xsecjj(X,Xsec)
         Weig = Xsec
      elseif(ipro.eq.102) then
         Call Xhiggs(X,Xsec)
         Weig = Xsec
      elseif(ipro.eq.500) then 
         Call XV1(X,Xsec)
         Weig = Xsec
      elseif(ipro.eq.501) then
         Call XVq(X,Xsec)
         Weig = Xsec
      elseif(ipro.ge.503.and.ipro.le.509) then
         Call XVQQ(X,Xsec)
         Weig = Xsec
      elseif(ipro.eq.510) then
          write(6,*) ' fxn1: process not yet implemented: STOP '
          stop
c         Call XVq(X,Xsec)
c         Weig = Xsec
      elseif(ipro.eq.511) then
         KWZ = 24
         Call XVq(X,Xsec)
         WeigW1 = Xsec
         KWZ = -24
         Call XVq(X,Xsec)
         WeigW2 = Xsec
         Weig = WeigW1 + WeigW2
      elseif(ipro.ge.513.and.ipro.le.519) then
         If(Ipro.EQ.513) Then 
            KWZ = 24
            KWZP1 = -2
            KWZP2 = 1
            Call XVQQ(X,Xsec)
            WeigW1 = Xsec
            KWZP1 = -2
            KWZP2 = 3
            Call XVQQ(X,Xsec)
            WeigW2 = Xsec
            KWZP1 = -2
            KWZP2 = 5
            Call XVQQ(X,Xsec)
            WeigW3 = Xsec
            KWZ = -24
            KWZP1 = 2
            KWZP2 = -1
            Call XVQQ(X,Xsec)
            WeigW4 = Xsec
            KWZP1 = 2
            KWZP2 = -3
            Call XVQQ(X,Xsec)
            WeigW5 = Xsec
            KWZP1 = 2
            KWZP2 = -5
            Call XVQQ(X,Xsec)
            WeigW6 = Xsec
            Weig = WeigW1 + WeigW2 + WeigW3 + WeigW4 + WeigW5 + WeigW6
         Else If(Ipro.EQ.514) Then 
            KWZ = 24
            KWZP1 = -IHFLA
            KWZP2 = IHFLA-3
            Call XVQQ(X,Xsec)
            WeigW1 = Xsec
            KWZP1 = -IHFLA
            KWZP2 = IHFLA-1
            Call XVQQ(X,Xsec)
            WeigW2 = Xsec
            KWZP1 = -IHFLA
            KWZP2 = IHFLA+1
            Call XVQQ(X,Xsec)
            WeigW3 = Xsec
            KWZ = -24
            KWZP1 = IHFLA
            KWZP2 = -(IHFLA-3)
            Call XVQQ(X,Xsec)
            WeigW4 = Xsec
            KWZP1 = IHFLA
            KWZP2 = -(IHFLA-1)
            Call XVQQ(X,Xsec)
            WeigW5 = Xsec
            KWZP1 = IHFLA
            KWZP2 = -(IHFLA+1)
            Call XVQQ(X,Xsec)
            WeigW6 = Xsec
            Weig = WeigW1 + WeigW2 + WeigW3 + WeigW4 + WeigW5 + WeigW6
         Endif
      else

         Call Xsec1(X,Xsec)

         Weigl = Xsec
         Weig = Xsec
         ncall = ncall + 1
         IF(NFLAV.GT.3.and.IPRO.EQ.10) THEN
            IMEPHI = 1
            IPRO=11
            IHFLA = 4
            Call Xsec1(X,Xsec)

            Weigc = Xsec

            WEIG = WEIGL + WEIGC
            ncall = ncall + 1
            IF(NFLAV.EQ.5) THEN
               IMEPHI = 1
               IPRO=11
               IHFLA = 5
               Call Xsec1(X,Xsec)

               Weigb = Xsec

               WEIG = WEIGL + WEIGC + WEIGB
               ncall = ncall + 1
            ENDIF
            IF(NFLAV.EQ.6) THEN
               IMEPHI = 1
               IPRO=11
               IHFLA = 6
               Call Xsec1(X,Xsec)

               Weigt = Xsec

               WEIG = WEIGL + WEIGC + WEIGB + WEIGT
               ncall = ncall + 1
            ENDIF
         ENDIF
      endif
      IMEPHI = 0
      IPRO = IPROO
      FXN1=WEIG
      WEIGHT = WEIG
      RETURN
      END
