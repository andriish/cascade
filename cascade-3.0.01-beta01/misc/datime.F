*
* $Id: datime.F,v 1.1.1.1 2008/10/12 14:18:02 jung Exp $
*
* $Log: datime.F,v $
* Revision 1.1.1.1  2008/10/12 14:18:02  jung
* importing sources
*
* Revision 1.1.1.1  2008/10/12 13:56:30  jung
* importing sources
*
* Revision 1.1.1.1  1996/02/15 17:50:39  mclareni
* Kernlib
*
*
      SUBROUTINE DATIME (ID,IT)
C
C CERN PROGLIB# Z007    DATIME          .VERSION KERNMPW  1.00  891208
C
      integer*4 today(3)
      call idate(today)   ! today(1)=day, (2)=month, (3)=year 
      iday=today(1)
      imonth=today(2)
      iyear=today(3)
      ID = IYEAR * 10000 + IMONTH * 100 +IDAY
      IT = IHH * 100 + IMM
      IYEAR = IYEAR + 1900
c      write(6,*) ' datime ',iday,imonth,iyear
      RETURN
      END
