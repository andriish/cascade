C...CAUPINIT
 
      SUBROUTINE CAUPINIT
      Implicit None
#include "caslujets.inc"
#include "casluco.inc"
#include "caseffic.inc"
#include "caspara.inc"
#include "casherup.inc"
C...CASCADE process initialization commonblock.
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &LPRUP(MAXPUP)
      SAVE /HEPRUP/
	
	Integer MSTP,MSTI
	Double Precision PARP,PARI
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      LOGICAL ex
      write(6,*) ' CASCADE caupinit'	
      If(ILHA.eq.-1.or.ILHA.ge.1) then
C Decide on weighting strategy, events have weight = 1 
	IDWTUP = 3
C Nr of external processes
      NPRUP = 1
c set incoming beams
      IDBMUP(1) = K(2,2)
      IDBMUP(2) = K(1,2)
	
	LPRUP(1)=IPRO
	EBMUP(1) = P(2,4)+1D-5*P(2,4)
	EBMUP(2) = P(1,4)+1D-5*P(1,4)
      PDFGUP(1) = -1
      PDFGUP(2) = -1
	PDFSUP(1) =-1.
	PDFSUP(2) =-1.
	XSECUP(1) = Avgi*1000.D0
	XERRUP(1) = SD*1000.D0
	XMAXUP(1) = 1.
      Endif
C write xsection information
      IF(ILHA.ge.10) then
c switch off color reconnections 
         MSTP(95) = 0
c set maximum number of partons
         MSTP(126) = 200
c call directly PYTHIA to do final state radiation and fragementation	
	   MSTP(61)=0 ! initial state PS
c         mstp(71)=1
c	   if(ifps.le.1) mstp(71)=0
c      MSTP(71)=0 ! final state PS
         MSTP(81)=0 ! multiple interactions
c      MSTP(111)=0 ! fragmentation 
	 IF(NFRAG.eq.0) MSTP(111)=0
C...Initialize.
         CALL PYINIT('USER',' ',' ',0D0)
        IDWTUP   = 1
        write(6,*) ' end of pythia ini ',mstp(81)
	else
         inquire(FILE='cascade.gen',EXIST=ex)
         if(ex) then
	   Write(6,*) ' File cascade.gen exists, will be overwritten '
         OPEN(21,FILE='cascade.gen',FORM='unformatted',STATUS=
     +         'OLD',ERR=100)
	   else
	   Write(6,*) ' Open new file cascade.gen '
         OPEN(21,FILE='cascade.gen',FORM='unformatted',STATUS=
     +         'NEW',ERR=100)
         endif
         if(ILHA.ne.-2) then
         WRITE(21,ERR=100) IDBMUP(1),IDBMUP(2),EBMUP(1),EBMUP(2),
     &   PDFGUP(1),PDFGUP(2),PDFSUP(1),PDFSUP(2),
     &   IDWTUP,NPRUP,XSECUP(1),XERRUP(1),XMAXUP(1),
     &   LPRUP(1)
         else
         READ(21,ERR=100) IDBMUP(1),IDBMUP(2),EBMUP(1),EBMUP(2),
     &   PDFGUP(1),PDFGUP(2),PDFSUP(1),PDFSUP(2),
     &   IDWTUP,NPRUP,XSECUP(1),XERRUP(1),XMAXUP(1),
     &   LPRUP(1)
         endif
      endif 
      	
      RETURN
100	Write(6,*) ' Error! File open failed for cascade.gen  STOP '
      STOP
      END 
