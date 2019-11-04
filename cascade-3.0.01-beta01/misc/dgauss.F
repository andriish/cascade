*
* $Id: dgauss.F,v 1.1.1.1 2008/10/12 14:18:02 jung Exp $
*
* $Log: dgauss.F,v $
* Revision 1.1.1.1  2008/10/12 14:18:02  jung
* importing sources
*
* Revision 1.1.1.1  2008/10/12 13:56:30  jung
* importing sources
*
* Revision 1.1.1.1  1996/02/15 17:48:17  mclareni
* Kernlib
*
*
      DOUBLE PRECISION FUNCTION DGAUSS(F,A,B,EPS)
      DOUBLE PRECISION F,A,B,EPS
      DOUBLE PRECISION W(12),X(12),AA,BB,C1,C2,U,S8,S16,CONST
      LOGICAL MFLAG,RFLAG
      EXTERNAL F
C
C     ******************************************************************
C
C     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.
C
C     DGAUSS IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER
C     EPS.
C
C     ******************************************************************
C
      DATA W / 0.10122 85362 90376 259D0,
     1         0.22238 10344 53374 471D0,
     2         0.31370 66458 77887 287D0,
     3         0.36268 37833 78361 983D0,
     4         0.27152 45941 17540 949D-1,
     5         0.62253 52393 86478 929D-1,
     6         0.95158 51168 24927 848D-1,
     7         0.12462 89712 55533 872D0,
     8         0.14959 59888 16576 732D0,
     9         0.16915 65193 95002 538D0,
     A         0.18260 34150 44923 589D0,
     B         0.18945 06104 55068 496D0/
 
      DATA X / 0.96028 98564 97536 232D0,
     1         0.79666 64774 13626 740D0,
     2         0.52553 24099 16328 986D0,
     3         0.18343 46424 95649 805D0,
     4         0.98940 09349 91649 933D0,
     5         0.94457 50230 73232 576D0,
     6         0.86563 12023 87831 744D0,
     7         0.75540 44083 55003 034D0,
     8         0.61787 62444 02643 748D0,
     9         0.45801 67776 57227 386D0,
     A         0.28160 35507 79258 913D0,
     B         0.95012 50983 76374 402D-1/
C
C     ******************************************************************
C
C  START.
      DGAUSS=0.0D0
      IF(B.EQ.A) RETURN
      CONST=0.005D0/(B-A)
      BB=A
C
C  COMPUTATIONAL LOOP.
    1 AA=BB
      BB=B
    2    C1=0.5D0*(BB+AA)
         C2=0.5D0*(BB-AA)
         S8=0.0D0
         DO 3 I=1,4
            U=C2*X(I)
            S8=S8+W(I)*(F(C1+U)+F(C1-U))
    3    CONTINUE
         S8=C2*S8
         S16=0.0D0
         DO 4 I=5,12
            U=C2*X(I)
            S16=S16+W(I)*(F(C1+U)+F(C1-U))
    4    CONTINUE
         S16=C2*S16
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5
         BB=C1
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2
      DGAUSS=0.0D0
ccc      CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG)
      IF(MFLAG) THEN
         IF(LGFILE.EQ.0) THEN
            WRITE(*,6)
         ELSE
            WRITE(LGFILE,6)
         ENDIF
      ENDIF
      IF(.NOT. RFLAG) CALL ABEND
      RETURN
    5 DGAUSS=DGAUSS+S16
      IF(BB.NE.B) GO TO 1
      RETURN
C
    6 FORMAT( 4X, 'FUNCTION DGAUSS ... TOO HIGH ACCURACY REQUIRED')
      END
