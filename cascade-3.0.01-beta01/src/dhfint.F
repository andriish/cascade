*
* $Id: dhfint.F,v 1.1.1.1 2008/10/12 14:18:09 jung Exp $
*
* $Log: dhfint.F,v $
* Revision 1.1.1.1  2008/10/12 14:18:09  jung
* importing sources
*
* Revision 1.1.1.1  2008/10/12 13:56:39  jung
* importing sources
*
* Revision 1.1.1.1  1996/02/15 17:48:36  mclareni
* Kernlib
*
*
          FUNCTION DHFINT(NARG,ARG,NENT,ENT,TABLE)
          Implicit None
C
C   INTERPOLATION ROUTINE. AUTHOR C. LETERTRE.
C   MODIFIED BY B. SCHORR, 1.07.1982.
C

          INTEGER   NENT(500)
          Double Precision     ARG(500),   ENT(500),TABLE(1000000)
          INTEGER   INDEX(32)
          Double Precision       WEIGHT(32)
          Double Precision DHFINT
          Integer NARG,LMAX,ISTEP,KNOTS,N,NDIM,ISHIFT,I,K
          Integer LMIN,LOCA,LOCB,LOCC,LGFILE
          Double Precision X,H,ETA
          LOGICAL   MFLAG,    RFLAG
          DHFINT  =  0.
          IF(NARG .LT. 1  .OR.  NARG .GT. 5)  GOTO 300
          LMAX      =  0
          ISTEP     =  1
          KNOTS     =  1
          INDEX(1)  =  1
          WEIGHT(1) =  1.
          DO 100    N  =  1, NARG
             X     =  ARG(N)
             NDIM  =  NENT(N)
             LOCA  =  LMAX
             LMIN  =  LMAX + 1
             LMAX  =  LMAX + NDIM
             IF(NDIM .GT. 2)  GOTO 10
             IF(NDIM .EQ. 1)  GOTO 100
             H  =  X - ENT(LMIN)
             IF(H .EQ. 0.)  GOTO 90
             ISHIFT  =  ISTEP
             IF(X-ENT(LMIN+1) .EQ. 0.)  GOTO 21
             ISHIFT  =  0
             ETA     =  H / (ENT(LMIN+1) - ENT(LMIN))
             GOTO 30
  10         LOCB  =  LMAX + 1
  11         LOCC  =  (LOCA+LOCB) / 2
             IF(X-ENT(LOCC))  12, 20, 13
  12         LOCB  =  LOCC
             GOTO 14
  13         LOCA  =  LOCC
  14         IF(LOCB-LOCA .GT. 1)  GOTO 11
             LOCA    =  MIN0( MAX0(LOCA,LMIN), LMAX-1 )
             ISHIFT  =  (LOCA - LMIN) * ISTEP
             ETA     =  (X - ENT(LOCA)) / (ENT(LOCA+1) - ENT(LOCA))
             GOTO 30
  20         ISHIFT  =  (LOCC - LMIN) * ISTEP
  21         DO 22  K  =  1, KNOTS
                INDEX(K)  =  INDEX(K) + ISHIFT
  22            CONTINUE
             GOTO 90
  30         DO 31  K  =  1, KNOTS
                INDEX(K)         =  INDEX(K) + ISHIFT
                INDEX(K+KNOTS)   =  INDEX(K) + ISTEP
                WEIGHT(K+KNOTS)  =  WEIGHT(K) * ETA
                WEIGHT(K)        =  WEIGHT(K) - WEIGHT(K+KNOTS)
  31            CONTINUE
             KNOTS  =  2*KNOTS
  90         ISTEP  =  ISTEP * NDIM
 100         CONTINUE
          DO 200    K  =  1, KNOTS
             I  =  INDEX(K)
             DHFINT  =  DHFINT + WEIGHT(K) * TABLE(I)
 200         CONTINUE
          RETURN
 300      Continue
c          CALL KERMTR('E104.1',LGFILE,MFLAG,RFLAG)
          IF(MFLAG) THEN
             IF(LGFILE .EQ. 0) THEN
                WRITE(*,1000) NARG
             ELSE
                WRITE(LGFILE,1000) NARG
             ENDIF
          ENDIF
          IF(.NOT. RFLAG) CALL ABEND
          RETURN
1000      FORMAT('  HFUNCTION DHFINT ... NARG =',I6,
     +              '  NOT WITHIN RANGE')
          END
      
