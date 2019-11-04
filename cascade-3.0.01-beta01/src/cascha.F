      SUBROUTINE cascha
*#**********************************************************************
*#
*#    SUBROUTINE cascha
*#
*#**********************************************************************
      character *4   param 
*
	Integer nloop,nmax
	Parameter (nmax=1000)
      character *4   para 
      character*512 Cdef,Cea
	Real Rval 
      character*512 Cval 
	Integer Ival,Ld,Le
	Common/steering/Nloop,Ld(nmax),Le(nmax),Ival(nmax),
     &	para(nmax),Rval(nmax),Cval(nmax)

#include "casparas.inc"
#include "casluco.inc"
#include "caspara5.inc"
#include "caspart.inc"
#include "casvalues.inc"
#include "casshwr.inc"
#include "casprkt.inc"
#include "casscal.inc"
#include "casglu.inc"
#include "casmult.inc"
#include "casherup.inc"
#include "caspsipol.inc"
#include "caslhe.inc"
      Integer Iremn
      common/casflags/Iremn	

      character*512 PDFPATH
      common/CASPDF/PDFPATH
      Integer ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON /CASCOL/ ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG

      Double Precision ACC1,ACC2
      COMMON /INTEGR/ ACC1,ACC2,IINT,NCB
      Integer IGLU
	Common/CAGLUON/Iglu
      DATA LND /0/ 

C---Determine the machine accuacy EPSMAC (D= 2.E-06)
      EPSMAC = 0.5
      DO 10 I= 1, 100
         EPSMAC = EPSMAC * 0.5
         IF ((1.0+EPSMAC) .EQ. 1.0) GO TO 20
   10 CONTINUE
      EPSMAC = 1.0E-6
   20 EPSMAC = 2.0 * EPSMAC
      WRITE(6,10000) EPSMAC
10000 FORMAT(1X,/,' computed machine accuracy = ',G10.2,/)

C---Interpret the text containing the parameters to be changed
      If(Nloop.ge.1) Then
	      Do 30 I=1,Nloop

            PARAM = Para(I)
            INT = Ival(I)
            REA = Rval(I)
            CEA = Cval(I)
            

            IF (PARAM.EQ.'KBE1') THEN
               IDEF = KE
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  KE = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRE1') THEN
               IDEF = Ires(1)
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Ires(1) = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'KBE2') THEN
               IDEF = KP
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  KP = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRE2') THEN
               IDEF = Ires(2)
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Ires(2) = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'NFRA') THEN
               IDEF = NFRAG
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  NFRAG = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IPST') THEN
               IDEF = IPST
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IPST = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IPOL') THEN
               IDEF = Ipsipol
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Ipsipol = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'I23S') THEN
               IDEF = I23S
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  I23S = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IFPS') THEN
               IDEF = IFPS
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IFPS = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ITIM') THEN
               IDEF = ITIMSHR
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ITIMSHR = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRAS') THEN
               IDEF = IRUNA
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRUNA = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRAM') THEN
               IDEF = IRUNAEM
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRUNAEM = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IQ2S') THEN
               IDEF = IQ2
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IQ2 = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IPRO') THEN
               IDEF = IPRO
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IPRO = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'NFLA') THEN
               IDEF = NFLAV
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  NFLAV = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'INTE') THEN
               IDEF = INTER
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  INTER = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'IHFL') THEN
               IDEF = IHFLA
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IHFLA = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'IRPA') THEN
               IDEF = IRPA
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRPA = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRPB') THEN
               IDEF = IRPB
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRPB = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRPC') THEN
               IDEF = IRPC
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRPC = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRPD') THEN
               IDEF = IRPD
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRPD = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRPE') THEN
               IDEF = IRPE
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRPE = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRPF') THEN
               IDEF = IRPF
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRPF = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRPG') THEN
               IDEF = IRPG
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRPG = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'ISEH') THEN
               IDEF = ISEMIH
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ISEMIH = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ICCF') THEN
               IDEF = ICCFM
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ICCFM = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IGLU') THEN
               IDEF = IGLU
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IGLU = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ICGL') THEN
               IDEF = ICOLGLU
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ICOLGLU = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IINT') THEN
               IDEF= IINT
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IINT = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'NCAL') THEN
               IDEF= NCB
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  NCB = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IREM') THEN
               IDEF= IRspl
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRspl = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'MULT') THEN
               IDEF= Mult
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Mult = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'ILHA') THEN
               IDEF = Ilha
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Ilha = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IFIN') THEN
               IDEF = IFINAL
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IFINAL = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRMN') THEN
               IDEF = Iremn
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Iremn = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
c               write(6,*) ' cascha ',Iremn
            ENDIF
            IF (PARAM.EQ.'ILHE') THEN
               IDEF = Ilhe
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Ilhe = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ITMW') THEN
               IDEF = iTMDw
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  iTMDw = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ISCA') THEN
               IDEF = Iscale
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Iscale = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'PBE2') THEN
               RDEF = PPIN
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PPIN = REA
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'PBE1') THEN
               RDEF = PLEPIN
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PLEPIN = REA
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'QMIN') THEN
               RDEF = SNGL(QMI)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  QMI = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'QMAX') THEN
               RDEF = SNGL(QMA)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  QMA = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'YMIN') THEN
               RDEF = SNGL(YMI)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  YMI = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'YMAX') THEN
               RDEF = SNGL(YMA)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  YMA = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'THEM') THEN
               RDEF = SNGL(THEMA)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  THEMA = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'THEI') THEN
               RDEF = SNGL(THEMI)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  THEMI = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'PTCU') THEN
               RDEF = SNGL(PT2CUT(Iabs(IPRO)))
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PT2CUT(IPRO) = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF


            IF (PARAM.EQ.'ACC1') THEN
               RDEF = ACC1
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  ACC1 = REA
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ACC2') THEN
               RDEF = ACC2
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  ACC2 = REA
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'XKT1') THEN
               RDEF = SNGL(PRKT1)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PRKT1 = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'XKTR') THEN
               RDEF = SNGL(PRKTREM)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PRKTREM = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'XKT2') THEN
               RDEF = SNGL(PRKT2)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PRKT2 = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'XKTR') THEN
               RDEF = SNGL(PRKTREM)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PRKTREM = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'SCAL') THEN
               RDEF = SNGL(SCALFA)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  SCALFA = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'SCAF') THEN
               RDEF = SNGL(SCALFAF)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  SCALFAF = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'UPDF') THEN
               CDEF = PDFPATH
               IF (CDEF.NE.CEA) THEN
                  PDFPATH = CEA
                  PRINT 10500, PARAM,LND,CDEF,CEA
               ELSE
                  PRINT 10600, PARAM,LND,CDEF,CEA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'OUT ') THEN
               CDEF = hepmcout
               IF (CDEF.NE.CEA) THEN
                  hepmcout = CEA
                  PRINT 10500, PARAM,LND,CDEF,CEA
               ELSE
                  PRINT 10600, PARAM,LND,CDEF,CEA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'CLHE ') THEN
               CDEF = lhefile
               IF (CDEF.NE.CEA) THEN
                  lhefile = CEA
                  PRINT 10500, PARAM,LND,CDEF,CEA
               ELSE
                  PRINT 10600, PARAM,LND,CDEF,CEA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'CWGT ') THEN
               CDEF = lhefile
               IF (CDEF.NE.CEA) THEN
                  Collfile = CEA
                  PRINT 10500, PARAM,LND,CDEF,CEA
               ELSE
                  PRINT 10600, PARAM,LND,CDEF,CEA
               ENDIF
            ENDIF

   30    CONTINUE
      ENDIF

10100 FORMAT(' cascha: ',A4,'(',I4,')',
     +       ' changed from the default ',I8,' TO',I8)
10200 FORMAT(' cascha (WARNING):',A4,'(',I4,')',
     +       ' not changed from the default ',I8,' TO',I8)
10300 FORMAT(' cascha: ',A4,'(',I4,')',
     +       ' changed from the default ',E10.3,' TO',E10.3)
10400 FORMAT(' cascha (WARNING): ',A4,'(',I4,')',
     +       ' not changed from the default ',E10.3,' TO',E10.3)
10500 FORMAT(' cascha: ',A4,'(',I4,')',
     +       ' changed from the default ',A50,' TO ',A50)
10600 FORMAT(' cascha (WARNING): ',A4,'(',I4,')',
     +       ' not changed from the default ',A50,' TO ',A50)

      RETURN
      END
