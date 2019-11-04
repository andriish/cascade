
      SUBROUTINE PYTCHA

*#**********************************************************************
*#
*#    SUBROUTINE PYTCHA
*#
*# PURPOSE: Interpret the seering to change the default
*#          parameters of PYTHIA 5.7.
*#          note: only those parameters can be changed which the user
*#                may change according to the PYTHIA manual.
*#
*#**********************************************************************
      character *4   param 
*
	Integer nloop,nmax
	Parameter (nmax=1000)
      character *4   para 
	Real Rval 
      character*50 Cval 
	Integer Ival,Ld,Le
	Common/steering/Nloop,Ld(nmax),Le(nmax),Ival(nmax),
     &	para(nmax),Rval(nmax),Cval(nmax)


      Double Precision paru,parj,pmas,parf,vckm,brat,VINT
      Integer MINT
      CHARACTER CHAF*16
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
 	COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT4/CHAF(500,2)
      COMMON/PYINT1/MINT(400),VINT(400)

C---additional PYTHIA commons (not part of standard PYTHIA)
      CHARACTER*8     FRAME,BEAM,TARGET
      Double Precision        win,pbeam,ptarge
      COMMON /PBEAM1/ FRAME,BEAM,TARGET
      COMMON /PBEAM2/ WIN,PBEAM(3),PTARGE(3)
      SAVE /PBEAM1/,/PBEAM2/

      Double Precision ckin,parp,pari
C---PYTHIA commons
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      SAVE /PYSUBS/,/PYPARS/

      CHARACTER*4     CHAINT


C---Determine the machine accuacy EPSMAC (D= 2.E-06)
      EPSMAC = 0.5
      DO 33 I= 1, 100
      EPSMAC = EPSMAC * 0.5
      IF ((1.0+EPSMAC) .EQ. 1.0)  GO TO 35
   33 CONTINUE
      EPSMAC = 1.0E-6
   35 EPSMAC = 2.0 * EPSMAC
      WRITE(6,330) EPSMAC
 330  FORMAT(1X,/,' COMPUTED MACHINE ACCURACY = ',G10.2,/)

C---Interpret the text containing the parameters to be changed
      If(Nloop.ge.1) Then
	      Do 30 I=1,Nloop

            LND = Ld(I)
            LNE = Le(I)

            PARAM = Para(I)
            INT = Ival(I)
            REA = Rval(I)

            IF (PARAM.EQ.'MSEL') THEN
                IDEF = MSEL
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    MSEL = INT
                    PRINT 91001, PARAM,LND,IDEF,INT
                ELSE
                    PRINT 92001, PARAM,LND,IDEF,INT
                ENDIF
            ENDIF

            IF (PARAM.EQ.'MSUB') THEN
                IDEF = MSUB(LND)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    MSUB(LND) = INT
                    PRINT 91001, PARAM,LND,IDEF,INT
                ELSE
                    PRINT 92001, PARAM,LND,IDEF,INT
                ENDIF
            ENDIF

            IF (PARAM.EQ.'MSTP') THEN
                IDEF = MSTP(LND)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    MSTP(LND) = INT
                    PRINT 91001, PARAM,LND,IDEF,INT
                ELSE
                    PRINT 92001, PARAM,LND,IDEF,INT
                ENDIF
            ENDIF

            IF (PARAM.EQ.'CKIN') THEN
                RDEF = CKIN(LND)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    CKIN(LND) = REA
                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'WIN ') THEN
                RDEF = WIN
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    WIN = dble(REA)
                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'PBEA') THEN
                RDEF = PBEAM(LND)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    PBEAM(LND) = REA
                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'PTAR') THEN
                RDEF = PTARGE(LND)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    PTARGE(LND) = REA
                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'PARP') THEN
                RDEF = PARP(LND)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    PARP(LND) = dble (REA)
                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'KFIN') THEN
                IDEF = KFIN(LND,LNE)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    KFIN(LND,LNE) = REA
                    PRINT 91003, PARAM,LND,LNE,IDEF,INT
                ELSE
                    PRINT 92003, PARAM,LND,LNE,IDEF,INT
                ENDIF
            ENDIF


            IF (PARAM.EQ.'MSTU') THEN
                IDEF = MSTU(LND)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    MSTU(LND) = INT
                    PRINT 91001, PARAM,LND,IDEF,INT
                ELSE
                    PRINT 92001, PARAM,LND,IDEF,INT
                ENDIF
            ENDIF

            IF (PARAM.EQ.'MINT') THEN
                IDEF = MINT(LND)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    MINT(LND) = INT
                    PRINT 91001, PARAM,LND,IDEF,INT
                ELSE
                    PRINT 92001, PARAM,LND,IDEF,INT
                ENDIF
            ENDIF
            IF (PARAM.EQ.'MSTJ') THEN
                IDEF = MSTJ(LND)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    MSTJ(LND) = INT
                    PRINT 91001, PARAM,LND,IDEF,INT
                ELSE
                    PRINT 92001, PARAM,LND,IDEF,INT
                ENDIF
            ENDIF

            IF (PARAM.EQ.'PARU') THEN
                RDEF = PARU(LND)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    PARU(LND) = dble(REA)
                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'PARJ') THEN
                RDEF = PARJ(LND)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    PARJ(LND)= dble(REA)
                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'KCHG') THEN
                IDEF = KCHG(LND,LNE)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    KCHG(LND,LNE) = INT
                    PRINT 91003, PARAM,LND,LNE,IDEF,INT
                ELSE
                    PRINT 92003, PARAM,LND,LNE,IDEF,INT
                ENDIF
            ENDIF

            IF (PARAM.EQ.'PMAS') THEN
                RDEF = PMAS(LND,LNE)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    PMAS(LND,LNE)=dble(REA)
                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'PARF') THEN
                RDEF = PARF(LND)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    PARF(LND) = dble(REA)
                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'VCKM') THEN
                RDEF = VCKM(LND,LNE)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    VCKM(LND,LNE) = dble(REA)
                    PRINT 91002, PARAM,LND,LNE,RDEF,REA
                ELSE
                    PRINT 92002, PARAM,LND,LNE,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'MDCY') THEN
                IDEF = MDCY(LND,LNE)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    MDCY(LND,LNE) = INT
                    PRINT 91004, PARAM,LND,LNE,IDEF,INT
                ELSE
                    PRINT 92004, PARAM,LND,LNE,IDEF,INT
                ENDIF
            ENDIF

            IF (PARAM.EQ.'MDME') THEN
                IDEF = MDME(LND,LNE)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    MDME(LND,LNE) = INT
                    PRINT 91004, PARAM,LND,LNE,IDEF,INT
                ELSE
                    PRINT 92004, PARAM,LND,LNE,IDEF,INT
                ENDIF
            ENDIF
            IF (PARAM.EQ.'BRAT') THEN
                RDEF = BRAT(LND)
                IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                    BRAT(LND) = dble(REA)
                    PRINT 91002, PARAM,LND,RDEF,REA
                ELSE
                    PRINT 92002, PARAM,LND,RDEF,REA
                ENDIF
            ENDIF

            IF (PARAM.EQ.'KFDP') THEN
                IDEF = KFDP(LND,LNE)
                IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                    KFDP(LND,LNE) = INT
                    PRINT 91004, PARAM,LND,LNE,IDEF,INT
                ELSE
                    PRINT 92004, PARAM,LND,LNE,IDEF,INT
                ENDIF
            ENDIF


30      CONTINUE
      ENDIF

91001 FORMAT(' pytcha: ',A4,'(',I4,')',
     &       ' changed from the default ',I8,' TO',I8)
92001 FORMAT(' pytcha (warning):',A4,'(',I4,')',
     &       ' not changed from the default ',I8,' TO',I8)
91004 FORMAT(' pytcha: ',A4,'(',I4,')',
     &       ' changed from the default ',I8,I8,' TO',I8,I8)
92004 FORMAT(' pytcha (warning):',A4,'(',I4,')',
     &       ' not changed from the default ',I8,I8,' TO',I8,I8)
91002 FORMAT(' pytcha: ',A4,'(',I4,')',
     &       ' changed from the default ',E10.3,' TO',E10.3)
92002 FORMAT(' pytcha (warning):',A4,'(',I4,')',
     &       ' not changed from the default ',E10.3,' TO',E10.3)

91003 FORMAT(' pytcha:',
     &'CHAF(',I4,') changed from the default:',A4,A4)
92003 FORMAT(' pytcha:',
     &'CHAF(',I4,') not changed from the default :',A4,A4)
*91004 FORMAT('                                     TO:',A4,A4)


      RETURN
      END
