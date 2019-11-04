
      FUNCTION MINUTS()

************************************************************************
* Returns (not exactly)  number of minutes since 6.6.91                *
* to get different random number generator seed when h1 database is    *
* not written to.                                                      *
************************************************************************

      CALL DATIME(ID,IT)

      IDD=MOD(ID,100)
      ID=ID/100
      IMM=MOD(ID,100)
      IYY=ID/100
      JMM=MOD(IT,100)
      IHH=IT/100

      MINUTS=(((IYY*12+IMM)*31+IDD)*24+IHH)*60+JMM - 49024019

      RETURN
      END
