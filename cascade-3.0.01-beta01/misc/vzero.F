*
* $Id: vzero.F,v 1.1.1.1 2008/10/12 14:18:02 jung Exp $
*
* $Log: vzero.F,v $
* Revision 1.1.1.1  2008/10/12 14:18:02  jung
* importing sources
*
* Revision 1.1.1.1  2008/10/12 13:56:30  jung
* importing sources
*
* Revision 1.1.1.1  1996/02/15 17:50:18  mclareni
* Kernlib
*
*
      SUBROUTINE VZERO (A,N)
C
C CERN PROGLIB# F121    VZERO           .VERSION KERNFOR  4.40  940929
C ORIG. 01/07/71, modif. 24/05/87 to set integer zero
C                 modif. 25/05/94 to depend on QINTZERO
C
      DIMENSION A(*)
C
      IF (N.LE.0)  RETURN
      DO 9 I= 1,N
#if defined(CERNLIB_QINTZERO)
    9 A(I)= 0
#endif
#if !defined(CERNLIB_QINTZERO)
    9 A(I)= 0.
#endif
C
      RETURN
      END
