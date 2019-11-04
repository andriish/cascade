      SUBROUTINE COLORFLOW

      IMPLICIT NONE
#include "caslujets.inc"
#include "casgki.inc"
#include "caspara1.inc"
#include "casparas.inc"
#include "caspart.inc"
#include "caspara5.inc"
#include "casluco.inc"
#include "casludat1.inc"
#include "caspsipol.inc"

c..hju
      Integer ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON /CASCOL/ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
	Integer NVMAX
	Parameter (NVMAX=5000)
      Integer ks,ns
      COMMON /COLR/ KS(nvmax,5),NS
      Double Precision dcasrn

c      write(6,*) ' colres IRESPRO ',irespro,icolora
      IF(Ipro.eq.2) then
c color info for Quarkoinim production
        if(ires(1).eq.0) then
c color info for leptoproduction
           if(ipsiel1.eq.0) K(NS+5,1)=1
           K(NS+7,4)=(NS+3)*MSTU(5)+NS+6
           K(NS+7,5)=(NS+3)*MSTU(5)+NS+6
c           k(ns+5,1) = 1
        else
c color info for hadrproduction
           K(NS+1,1)=13
           K(NS+1,3)=1
           P(NS+2,3)=27
           P(NS+2,4)=27
c           K(NS+5,1)=1
           if(ipsiel1.eq.0) K(NS+5,1)=1
           K(NS+6,4)=(NS+7)*MSTU(5)
           K(NS+6,5)=(NS+7)*MSTU(5)
           K(NS+7,1)=14
           K(NS+7,3)=3
           K(NS+7,2)=21
           IF(dcasrn().GT.0.5) THEN
              K(NS+1,4)=27
              K(NS+1,5)=23
              K(NS+3,4)=21
              K(NS+3,5)=27
              K(NS+7,4)=(NS+1)*MSTU(5)+NS+6
              K(NS+7,5)=(NS+3)*MSTU(5)+NS+6
           ELSE
              K(NS+1,4)=23
              K(NS+1,5)=27
              K(NS+3,4)=27
              K(NS+3,5)=21
              K(NS+7,4)=(NS+3)*MSTU(5)+NS+6
              K(NS+7,5)=(NS+1)*MSTU(5)+NS+6
           ENDIF

        endif
      elseif(Ipro.eq.4) then
c  q from proton 1
         IF(K(21,1).GT.0.AND.K(23,2).EQ.21) THEN
c		write(6,*) ' here in colorflow '
                  K(NS+1,4) = 23
                  K(NS+1,5) = 0
                  K(NS+3,4)=21
                  K(NS+3,5)= 0 
                  K(NS+5,4)=0
                  K(NS+5,5)=0
                  K(NS+6,4)=0 
                  K(NS+6,5)= 0
                  K(NS+7,4)=(NS+3)*MSTU(5)+21
                  K(NS+7,5)=(NS+1)*MSTU(5)+23
                  K(NS+7,2) = K(NS+1,2)
c                  call pylist(2)
          endif

      elseif(Ipro.eq.21) then
         IF(K(21,2).eq.21.AND.K(23,2).EQ.21) THEN
c  		write(6,*) ' here in colorflow ',ipro
                  K(NS+1,4) = 23
                  K(NS+1,5) = 23
                  K(NS+3,4)=21
                  K(NS+3,5)=21
                  K(NS+5,4)=0
                  K(NS+5,5)=0
                  K(NS+6,4)=0 
                  K(NS+6,5)= 0
                  K(NS+7,4)=0
                  K(NS+7,5)= 0
                  K(NS+7,2) = 0
c                  call pylist(2)
          else
            write(6,*) ' colorflow: wrong configuration for ipro = ',ipro
          endif

      elseif(Ipro.eq.500.or.ipro.eq.501.or.ipro.eq.510.or.ipro.eq.511) then
c      write(6,*) ' colorflow: Z(W)q  start'
c  q from proton 1
         IF(K(23,2).GT.0.AND.K(21,2).EQ.21) THEN
            K(NS+1,4) = 27
            K(NS+1,5) = 23
            K(NS+3,4)=21
            K(NS+3,5)=21
            K(NS+5,1)=2
            K(NS+5,4)=0
            K(NS+5,5)=0
ctest            K(NS+6,4)=(NS+3)*MSTU(5)
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=(NS+7)*MSTU(5)
            K(NS+7,4)=(NS+1)*MSTU(5)+26
            K(NS+7,5)=(NS+1)*MSTU(5)+26
            K(NS+7,2) = K(NS+6,2)
c            write(6,*) ' colorflow: Z(W) 1 '
c            call pylist(2)
c  anti-q from proton 1
         ELSEIF(K(23,2).LT.0.AND.K(21,2).EQ.21) THEN
            K(NS+1,4) = 23
            K(NS+1,5) = 27
            K(NS+3,4)=21
            K(NS+3,5)=21
            K(NS+5,1)=2
            K(NS+5,4)=0
            K(NS+5,5)=0
ctest            K(NS+6,4)=(NS+3)*MSTU(5)
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=(NS+7)*MSTU(5)
            K(NS+7,4)=(NS+1)*MSTU(5)+26
            K(NS+7,5)=(NS+1)*MSTU(5)+26
            K(NS+7,2) = K(NS+6,2)
c            write(6,*) ' colorflow: Z(W) 2 '
c            call pylist(2)
c  q from proton 2
         ELSEIF(K(21,2).GT.0.AND.K(23,2).EQ.21) THEN
            K(NS+3,4) = 27
            K(NS+3,5) = 21
            K(NS+1,4)=23
            K(NS+1,5)=23
            K(NS+5,1)=2
            K(NS+5,4)=0
            K(NS+5,5)=0
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=(NS+7)*MSTU(5)
            K(NS+7,4)=(NS+3)*MSTU(5)+26
            K(NS+7,5)=(NS+3)*MSTU(5)+26
            K(NS+7,2) = K(NS+6,2)
c            write(6,*) ' colorflow: Z(W) 3 '
c            call pylist(2)
c  anti-q from proton 2
         ELSEIF(K(21,2).LT.0.AND.K(23,2).EQ.21) THEN
            K(NS+3,4) = 21
            K(NS+3,5) = 27
            K(NS+1,4)=23
            K(NS+1,5)=23
            K(NS+5,1)=2
            K(NS+5,4)=0
            K(NS+5,5)=0
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=(NS+7)*MSTU(5)
            K(NS+7,4)=(NS+3)*MSTU(5)+26
            K(NS+7,5)=(NS+3)*MSTU(5)+26
            K(NS+7,2) = K(NS+6,2)
c            write(6,*) ' colorflow: Z(W) 4 '
c            call pylist(2)
         endif
      else
      
      IF(IRESPRO.EQ.1) THEN
         K(NS+3,4)=21
         K(NS+3,5)=27
c qqbar event
         IF(K(NS+5,2).GT.0) THEN
            K(NS+5,4)=(NS+7)*MSTU(5)
            K(NS+5,5)=0
            K(NS+6,4)=0
            K(NS+6,5)=(NS+7)*MSTU(5)
            K(NS+1,4) = 27
            K(NS+1,5) = 23
            K(NS+7,4)=(NS+1)*MSTU(5)+25
            K(NS+7,5)=(NS+3)*MSTU(5)+26
         ELSE
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=0
            K(NS+5,4)=0
            K(NS+5,5)=(NS+7)*MSTU(5)
            K(NS+1,4) = 27
            K(NS+1,5) = 23
            K(NS+7,4)=(NS+1)*MSTU(5)+26
            K(NS+7,5)=(NS+3)*MSTU(5)+25
         ENDIF
      ELSEIF(IRESPRO.EQ.2) THEN
c glu glu event

         IF(ICOLORA.EQ.1) THEN
c color flow A
            K(NS+1,4) = 27
            K(NS+1,5) = 23
            K(NS+3,4)=21
            K(NS+3,5)=27
            K(NS+5,4)=(NS+7)*MSTU(5)
            K(NS+5,5)=(NS+6)*MSTU(5)
            K(NS+6,4)=(NS+5)*MSTU(5)
            K(NS+6,5)=(NS+7)*MSTU(5)
            K(NS+7,4)=(NS+1)*MSTU(5)+25
            K(NS+7,5)=(NS+3)*MSTU(5)+26
         ELSEIF(ICOLORA.EQ.2) THEN
c color flow B
            K(NS+1,4) = 27
            K(NS+1,5) = 23
            K(NS+3,4)=21
            K(NS+3,5)=27
            K(NS+5,4)=(NS+6)*MSTU(5)
            K(NS+5,5)=(NS+7)*MSTU(5)
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=(NS+5)*MSTU(5)
            K(NS+7,4)=(NS+1)*MSTU(5)+26
            K(NS+7,5)=(NS+3)*MSTU(5)+25
         ELSEIF(ICOLORA.EQ.3) THEN
c color flow C
            K(NS+1,4) = 26
            K(NS+1,5) = 25
            K(NS+3,4)=25
            K(NS+3,5)=26
            K(NS+5,4)=(NS+3)*MSTU(5)
            K(NS+5,5)=(NS+1)*MSTU(5)
            K(NS+6,4)=(NS+1)*MSTU(5)
            K(NS+6,5)=(NS+3)*MSTU(5)
            K(NS+7,4)=0
            K(NS+7,5)=0
            K(NS+7,1)=0
         ELSE
            write(6,*) ' wrong color configuartion ',ICOLORA
         ENDIF
c         call pylist(1)

      ELSEIF(IRESPRO.EQ.3) THEN
c qg event
C color configuartion A
         IF(ICOLORA.EQ.1) THEN
C antq from  proton 1
            IF(K(21,2).LT.0.AND.K(23,2).EQ.21) THEN
               IF(K(25,2).LT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+1,4) = 0
                  K(NS+1,5) = 23
                  K(NS+3,4)=21
                  K(NS+3,5)=27
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+6)*MSTU(5)
                  K(NS+6,4)=(NS+5)*MSTU(5)
                  K(NS+6,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+3)*MSTU(5)+25
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 1st '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).LT.0) THEN
                  K(NS+1,4) = 0
                  K(NS+1,5) = 23
                  K(NS+3,4)=21
                  K(NS+3,5)=27
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+5)*MSTU(5)
                  K(NS+5,4)=(NS+6)*MSTU(5)
                  K(NS+5,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 3rd '
               ENDIF
c  q from proton 1
            ELSEIF(K(21,2).GT.0.AND.K(23,2).EQ.21) THEN
               IF(K(25,2).GT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+1,4) = 23
                  K(NS+1,5) = 0
                  K(NS+3,4)=27
                  K(NS+3,5)=21
                  K(NS+5,4)=(NS+6)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=(NS+7)*MSTU(5)
                  K(NS+6,5)=(NS+5)*MSTU(5)
                  K(NS+7,4)=(NS+3)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 2nd '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).GT.0) THEN
                  K(NS+1,4) = 23
                  K(NS+1,5) = 0
                  K(NS+3,4)=27
                  K(NS+3,5)=21
                  K(NS+6,4)=(NS+5)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+5,4)=(NS+7)*MSTU(5)
                  K(NS+5,5)=(NS+6)*MSTU(5)
                  K(NS+7,4)=(NS+3)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 4th '
               ENDIF
C antq from proton 2
            ELSEIF(K(23,2).LT.0.AND.K(21,2).EQ.21) THEN
               IF(K(25,2).LT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+3,4) = 0
                  K(NS+3,5) = 21
                  K(NS+1,4)=23
                  K(NS+1,5)=27
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+6)*MSTU(5)
                  K(NS+6,4)=(NS+5)*MSTU(5)
                  K(NS+6,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+25
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 1st '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).LT.0) THEN
                  K(NS+3,4) = 0
                  K(NS+3,5) = 21
                  K(NS+1,4)=23
                  K(NS+1,5)=27
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+5)*MSTU(5)
                  K(NS+5,4)=(NS+6)*MSTU(5)
                  K(NS+5,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 3rd '
               ENDIF
c  q from proton 2
            ELSEIF(K(23,2).GT.0.AND.K(21,2).EQ.21) THEN
               IF(K(25,2).GT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+3,4) = 21
                  K(NS+3,5) = 0
                  K(NS+1,4)=27
                  K(NS+1,5)=23
                  K(NS+5,4)=(NS+6)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=(NS+7)*MSTU(5)
                  K(NS+6,5)=(NS+5)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 2nd '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).GT.0) THEN
                  K(NS+3,4) = 21
                  K(NS+3,5) = 0
                  K(NS+1,4)=27
                  K(NS+1,5)=23
                  K(NS+6,4)=(NS+5)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+5,4)=(NS+7)*MSTU(5)
                  K(NS+5,5)=(NS+6)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+25
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 4th '
               ENDIF
            ENDIF
C color configuartion B
         ELSE
C antq from proton 1
            IF(K(21,2).LT.0.AND.K(23,2).EQ.21) THEN
               IF(K(25,2).LT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+1,4) = 0
                  K(NS+1,5) = 26
                  K(NS+3,4)=26
                  K(NS+3,5)=25
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+3)*MSTU(5)
                  K(NS+6,4)=(NS+3)*MSTU(5)
                  K(NS+6,5)=(NS+1)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 1st '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).LT.0) THEN
                  K(NS+1,4) = 0
                  K(NS+1,5) = 25
                  K(NS+3,4)=25
                  K(NS+3,5)=26
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+3)*MSTU(5)
                  K(NS+5,4)=(NS+3)*MSTU(5)
                  K(NS+5,5)=(NS+1)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
                  K(NS+7,1)=0
c               write(6,*) ' 3rd '
               ENDIF
c  q from proton 1
            ELSEIF(K(21,2).GT.0.AND.K(23,2).EQ.21) THEN
               IF(K(25,2).GT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+1,4) = 26
                  K(NS+1,5) = 0
                  K(NS+3,4)=25
                  K(NS+3,5)=26
                  K(NS+5,4)=(NS+3)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=(NS+1)*MSTU(5)
                  K(NS+6,5)=(NS+3)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 2nd '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).GT.0) THEN
                  K(NS+1,4) = 25
                  K(NS+1,5) = 0
                  K(NS+3,4)=26
                  K(NS+3,5)=25
                  K(NS+6,4)=(NS+3)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+5,4)=(NS+1)*MSTU(5)
                  K(NS+5,5)=(NS+3)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 4th '
               ENDIF
C antq from proton 2
            ELSEIF(K(23,2).LT.0.AND.K(21,2).EQ.21) THEN
               IF(K(25,2).LT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+3,4) = 0
                  K(NS+3,5) = 26
                  K(NS+1,4)=26
                  K(NS+1,5)=25
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+1)*MSTU(5)
                  K(NS+6,4)=(NS+1)*MSTU(5)
                  K(NS+6,5)=(NS+3)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 1st '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).LT.0) THEN
                  K(NS+3,4) = 0
                  K(NS+3,5) = 25
                  K(NS+1,4)=25
                  K(NS+1,5)=26
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+1)*MSTU(5)
                  K(NS+5,4)=(NS+1)*MSTU(5)
                  K(NS+5,5)=(NS+3)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 3rd '
               ENDIF
c  q from proton 2
            ELSEIF(K(23,2).GT.0.AND.K(21,2).EQ.21) THEN
               IF(K(25,2).GT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+3,4) = 26
                  K(NS+3,5) = 0
                  K(NS+1,4)=25
                  K(NS+1,5)=26
                  K(NS+5,4)=(NS+1)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=(NS+3)*MSTU(5)
                  K(NS+6,5)=(NS+1)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 2nd '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).GT.0) THEN
                  K(NS+3,4) = 25
                  K(NS+3,5) = 0
                  K(NS+1,4)=26
                  K(NS+1,5)=25
                  K(NS+6,4)=(NS+1)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+5,4)=(NS+3)*MSTU(5)
                  K(NS+5,5)=(NS+1)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 4th '
               ENDIF
            ENDIF

         ENDIF
      ELSEIF(IRESPRO.EQ.4) THEN
c qqbar --> g g
         IF(ICOLORA.EQ.1) THEN
c color flow A
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 27
               K(NS+3,4)=27
               K(NS+3,5)=0
               K(NS+5,4)=(NS+6)*MSTU(5)
               K(NS+5,5)=(NS+7)*MSTU(5)
               K(NS+6,4)=(NS+7)*MSTU(5)
               K(NS+6,5)=(NS+5)*MSTU(5)
               K(NS+7,4)=(NS+3)*MSTU(5)+26
               K(NS+7,5)=(NS+1)*MSTU(5)+25
c            write(6,*) ' flow A, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 27
               K(NS+1,5) = 0
               K(NS+3,4)=0
               K(NS+3,5)=27
               K(NS+5,4)=(NS+7)*MSTU(5)
               K(NS+5,5)=(NS+6)*MSTU(5)
               K(NS+6,4)=(NS+5)*MSTU(5)
               K(NS+6,5)=(NS+7)*MSTU(5)
               K(NS+7,4)=(NS+1)*MSTU(5)+25
               K(NS+7,5)=(NS+3)*MSTU(5)+26
c            write(6,*) ' flow A, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSEIF(ICOLORA.EQ.2) THEN
c color flow B
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 27
               K(NS+3,4)=27
               K(NS+3,5)=0
               K(NS+5,4)=(NS+7)*MSTU(5)
               K(NS+5,5)=(NS+6)*MSTU(5)
               K(NS+6,4)=(NS+5)*MSTU(5)
               K(NS+6,5)=(NS+7)*MSTU(5)
               K(NS+7,4)=(NS+3)*MSTU(5)+25
               K(NS+7,5)=(NS+1)*MSTU(5)+26
c            write(6,*) ' flow B, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 27
               K(NS+1,5) = 0
               K(NS+3,4)=0
               K(NS+3,5)=27
               K(NS+5,4)=(NS+6)*MSTU(5)
               K(NS+5,5)=(NS+7)*MSTU(5)
               K(NS+6,4)=(NS+7)*MSTU(5)
               K(NS+6,5)=(NS+5)*MSTU(5)
               K(NS+7,4)=(NS+1)*MSTU(5)+26
               K(NS+7,5)=(NS+3)*MSTU(5)+25
c            write(6,*) ' flow B, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSE
            write(6,*) ' wrong color configuartion ',ICOLORA
         ENDIF
      ELSEIF(IRESPRO.EQ.5) THEN
c qqbar --> q qbar
         IF(ICOLORA.EQ.1) THEN
c color flow A
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 23
               K(NS+3,4)=21
               K(NS+3,5)=0
               K(NS+5,4)=0
               K(NS+5,5)=(NS+6)*MSTU(5)
               K(NS+6,4)=(NS+5)*MSTU(5)
               K(NS+6,5)=0
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow A, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 23
               K(NS+1,5) = 0
               K(NS+3,4)=0
               K(NS+3,5)=21
               K(NS+5,4)=(NS+6)*MSTU(5)
               K(NS+5,5)=0
               K(NS+6,4)=0
               K(NS+6,5)=(NS+5)*MSTU(5)
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow A, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSEIF(ICOLORA.EQ.2) THEN
c color flow B
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 27
               K(NS+3,4)=27
               K(NS+3,5)=0
               IF(K(25,2).GT.0) THEN
                  K(NS+5,4)=(NS+7)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+3)*MSTU(5)+25
                  K(NS+7,5)=(NS+1)*MSTU(5)+26
               ELSE
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+7)*MSTU(5)
                  K(NS+6,4)=(NS+7)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+7,4)=(NS+3)*MSTU(5)+26
                  K(NS+7,5)=(NS+1)*MSTU(5)+25
               ENDIF
c            write(6,*) ' flow B, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 27
               K(NS+1,5) = 0
               K(NS+3,4)=0
               K(NS+3,5)=27
               IF(K(25,2).GT.0) THEN
                  K(NS+5,4)=(NS+7)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+25
                  K(NS+7,5)=(NS+3)*MSTU(5)+26
               ELSE
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+7)*MSTU(5)
                  K(NS+6,4)=(NS+7)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+7,4)=(NS+1)*MSTU(5)+26
                  K(NS+7,5)=(NS+3)*MSTU(5)+25
               ENDIF
c            write(6,*) ' flow B, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSE
            write(6,*) ' wrong color configuartion ',ICOLORA
         ENDIF
      ELSEIF(IRESPRO.EQ.6) THEN
c q q --> q q
         IF(ICOLORA.EQ.1) THEN
c color flow A
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 26
               K(NS+1,5) = 0
               K(NS+3,4)=25
               K(NS+3,5)=0
               K(NS+5,4)=(NS+3)*MSTU(5)
               K(NS+5,5)=0
               K(NS+6,4)=(NS+1)*MSTU(5)
               K(NS+6,5)=0
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow A, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 26

               K(NS+3,4)=0
               K(NS+3,5)=25
               K(NS+5,4)=0
               K(NS+5,5)=(NS+3)*MSTU(5)
               K(NS+6,4)=0
               K(NS+6,5)=(NS+1)*MSTU(5)
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow A, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSEIF(ICOLORA.EQ.2) THEN
c color flow B
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 25
               K(NS+1,5) = 0

               K(NS+3,4)=26
               K(NS+3,5)=0
               K(NS+5,4)=(NS+1)*MSTU(5)
               K(NS+5,5)=0
               K(NS+6,4)=(NS+3)*MSTU(5)
               K(NS+6,5)=0
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow B, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 25

               K(NS+3,4)=0
               K(NS+3,5)=26
               K(NS+5,4)=0
               K(NS+5,5)=(NS+1)*MSTU(5)
               K(NS+6,4)=0
               K(NS+6,5)=(NS+3)*MSTU(5)
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow B, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSE
            write(6,*) ' wrong color configuartion ',ICOLORA
         ENDIF
      ELSE
         write(6,*) ' colres : irespro not implemented ',irespro
      ENDIF
      ENDIF
      RETURN
      END
