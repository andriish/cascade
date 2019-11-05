      PROGRAM CASMAIN
      Implicit None
c	Integer N1,N2
	Integer K1,K2,LUX
#include "casluco.inc"
#include "casparas.inc"
#include "caspara5.inc"
#include "casvalues.inc"
#include "caseffic.inc"
#include "caslujets.inc"
#include "casludat1.inc"
#include "casludat2.inc"
#include "casherup.inc"
#include "caslhe.inc"
      

      integer iev,isub
      double precision xsec,xsece

      Integer j
      real pta0, eta0, dna0, ptam, ptrms, etam, etrms, dnam, dnrms
      real ptrm0, etrm0, dnrm0
      real av, etas, pts, cost, eta, pt 
      real sec0, secrm0, sec, secrms
      real val, errval 
      
      
      Integer IGLU
	Common/CAGLUON/Iglu

      Integer ISEED,I
		
      Integer IVEC(25)
      Integer Irr	
	real timeleft
	Integer Minuts
	External Minuts
      External Pydata

      character*255 caseed
      character*255 test
      
      Double Precision weighthepmc,sumOfWeights
      
      character *512 hepmcout_env
      Integer II
      
      Real start,finish
C-----------------------
C...HepMC2
      character *132 hm2filnam
        Common/steerhm2/hm2filnam
      LOGICAL         hm2OUT
      COMMON /Qhm2LO/ hm2OUT
      INTEGER inithepmc,convhepmc
      INTEGER HEPMCSTATUS
C-----------------------
C...Rivet
      character *132 rvtfilnam
      character *132 rvtrun
      character *132 rvtana
        Integer Nrvtana
        Common/steerrvt/Nrvtana,rvtfilnam,rvtrun,rvtana(50)
      LOGICAL         rvtOUT
      COMMON /QrvtLO/ rvtOUT
      INTEGER rivetinit,rivetdone,rivetadd,rivetinitfirstevent,rivetrun
      INTEGER rvtHM,rivetinterfaceversion
      INTEGER RIVETSTATUS
C-----------------------
C...HepMC3
      character *132 hm3filnam
      Integer hm3filfor
        Common/steerhm3/hm3filnam,hm3filfor
      LOGICAL         hm3OUT
      COMMON /Qhm3LO/ hm3OUT

      INTEGER NMXHEP
      PARAMETER (NMXHEP=4000)
      COMMON /HEPEVT/  NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &                 JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),
     &                 VHEP(4,NMXHEP)
      INTEGER          NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
      DOUBLE PRECISION PHEP,VHEP

C...These variables will be used for conversion. Block size is different
C...in Pythia6 and in HepMC3, so the most simple portable way is to have
C... a second block of same size as in HepMC3 and  copy the content of 
C...block directly.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...Glue functions
      INTEGER HepMC3_new_writer,HepMC3_delete_writer
      INTEGER HepMC3_set_hepevt_address
      INTEGER HepMC3_convert_event,HepMC3_write_event,HepMC3_clear_event
      INTEGER HepMC3_set_attribute_int,HepMC3_set_attribute_double
      INTEGER HepMC3_set_pdf_info,HepMC3_set_cross_section
      INTEGER HepMC3_new_weight
      INTEGER HepMC3_set_weight_by_index
      EXTERNAL HepMC3_new_writer,HepMC3_delete_writer
      EXTERNAL HepMC3_set_hepevt_address
      EXTERNAL HepMC3_convert_event,HepMC3_write_event
      EXTERNAL HepMC3_clear_event
      EXTERNAL HepMC3_set_attribute_int,HepMC3_set_attribute_double
      EXTERNAL HepMC3_set_pdf_info,HepMC3_set_cross_section
      EXTERNAL HepMC3_new_weight
      EXTERNAL HepMC3_set_weight_by_index
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER NMXHEPL
      PARAMETER (NMXHEPL=10000)
      COMMON /HEPEVTL/  NEVHEPL,NHEPL,ISTHEPL(NMXHEPL),IDHEPL(NMXHEPL),
     &           JMOHEPL(2,NMXHEPL),JDAHEPL(2,NMXHEPL),PHEPL(5,NMXHEPL),
     &                 VHEPL(4,NMXHEPL)
      INTEGER          NEVHEPL,NHEPL,ISTHEPL,IDHEPL,JMOHEPL,JDAHEPL
      DOUBLE PRECISION PHEPL,VHEPL
      
      INTEGER OUTID, HEPMC3STATUS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      
      sumOfWeights = 0d0
      call cpu_time(start)

c  read default paramters and input steering file	
      call input_steer
            
c
c Start integration of x-section
c

      CALL CASCADE
      
c      call PYLIST(12)
CAV      CALL GETENV("HEPMCOUT", hepmcout_env)
CAV      if(hepmcout_env.ne.' ') then
CAV        hepmcout = trim(hepmcout_env)//char(0)
CAV        write(6,*) ' hepmcout from environment variable = ',trim(hepmcout)
CAV        else
CAV        hepmcout=trim(hepmcout)//char(0)
CAV        write(6,*) ' hepmcout from steering file = ',trim(hepmcout)
CAV      endif
C-- Init HepMC2      
      if (hm2OUT) then      
      HEPMCSTATUS=inithepmc(trim(hm2filnam)//char(0))
      endif
C-- Init HepMC3      
      if (hm3OUT) then
      OUTID=HepMC3_new_writer(0,hm3filfor,trim(hm3filnam)//char(0))
      HEPMC3STATUS=HepMC3_new_weight(OUTID,'Default'//char(0))
      NEVHEP=-123456
      HEPMC3STATUS=HepMC3_set_hepevt_address(NEVHEPL)
      endif
C-- Init Rivet
      if (rvtOUT) then
      rvtHM=rivetinterfaceversion()
      do i=1,Nrvtana
      RIVETSTATUS=rivetadd(trim(rvtana(i))//char(0))
      end do
      if (len(rvtrun).GT. 0) then
      RIVETSTATUS=rivetinit(trim(rvtrun)//char(0))
      end if
      endif      
      
      CALL GETENV("CASEED", caseed)
c note CASEED can have only 10 digits
c i.e   export CASEED='1234567890
      if(CASEED.eq.' ') then
         write(6,*) ' CASEED set = 12345 '
         iseed = 12345
      else
         read(caseed,'(I30)') iseed  
         write(6,*) ' CASEED from environment variable ', iseed
c         write(6,*) 'after read caseed ',   
      endif
c      WRITE (6,*) ' caseed: ',Iseed
c      LUX = 3
      LUX = 4
      K1=0
      K2=0
      CALL RLUXGO(LUX,ISEED,K1,K2)


c
c Print out result of integration of x-section
c
 
      CALL CAEND(1)
c
c Start event loop

      write(6,*) ' casmain ',nevent

      Do I=1,max(4,Nevent)

c generate an event
         CALL EVENT
         IF(MSTU(24).NE.0.or.LST(21).NE.0) goto 111

        iev = I
        isub = ipro
c xsec in CASCADE in nb, hepmc need pb
        xsec = avgi *1000.
        xsece = sd *1000.
        weighthepmc = 1d0
        if(ilha.ne.0) then
          weighthepmc = normhepmc*weight
c          weighthepmc = weight
c          write(6,*) 'weight ' ,weight,normhepmc
          sumOfWeights = sumOfWeights + weighthepmc
        endif
        write(6,*) ' casmain: event weight = ',weighthepmc,ilha 
        if ( hm2OUT .OR. rvtHM .EQ. 2 .AND.rvtOUT) then          
         HEPMCSTATUS=convhepmc(iev, isub, xsec, xsece, 
     +idP1,idP2,xidP1,xidP2,scale_process,xidP1,xidP2,pdf1,pdf2 )
C
C              HEPMCSTATUS=convhepmc(iev, isub, xsec, xsece,
C     &  KE, KP,x1,x2,q2pdfeval,xf1mom,xf2mom ,pdf1,pdf2  )
     
         end if
         if (rvtOUT) then
         if (I.EQ.1) then 
         write(*,*) ' once'
         RIVETSTATUS=rivetinitfirstevent(OUTID)
         end if
         RIVETSTATUS=rivetrun(OUTID)
         end if
                  if (hm3OUT .OR. rvtHM.EQ.3 .AND.rvtOUT) then
          NEVHEPL=IEV
          NHEPL=NHEP
           DO 500 J=1,NHEP
          ISTHEPL(J)=ISTHEP(J)
          IDHEPL(J)=IDHEP(J)
          JMOHEPL(1,J)=JMOHEP(1,J)
          JMOHEPL(2,J)=JMOHEP(2,J)
          JDAHEPL(1,J)=JDAHEP(1,J)
          JDAHEPL(2,J)=JDAHEP(2,J)
          PHEPL(1,J)=PHEP(1,J)
          PHEPL(2,J)=PHEP(2,J)
          PHEPL(3,J)=PHEP(3,J)
          PHEPL(4,J)=PHEP(4,J)
          PHEPL(5,J)=PHEP(5,J)
          VHEPL(1,J)=VHEP(1,J)
          VHEPL(2,J)=VHEP(2,J)
          VHEPL(3,J)=VHEP(3,J)
          VHEPL(4,J)=VHEP(4,J)
  500     CONTINUE  
          HEPMC3STATUS=HepMC3_convert_event(OUTID)
          HEPMC3STATUS=HepMC3_set_cross_section(OUTID,
     &  xsec, xsece,iev,iev)
          HEPMCSTATUS=hepmc3_set_weight_by_index(OUTID, weighthepmc,0)
C FIXME!
          HEPMC3STATUS=HepMC3_set_pdf_info(OUTID,
     &  idP1,idP2,xidP1,xidP2,scale_process,xidP1,xidP2,pdf1,pdf2)
          HEPMC3STATUS=HepMC3_set_attribute_int(OUTID,
     & -1,'mpi'//char(0))
          HEPMC3STATUS=HepMC3_set_attribute_int(OUTID,isub,
     &    'set_signal_process_id'//char(0))
          HEPMC3STATUS=HepMC3_set_attribute_int(OUTID,ISEED,
     &  'random_states1'//char(0))
          HEPMC3STATUS=HepMC3_set_attribute_double(OUTID,-1.0,
     &   'alphaEM'//char(0))
          HEPMC3STATUS=HepMC3_set_attribute_double(OUTID,-1.0,
     &   'alphaQCD'//char(0))
          HEPMC3STATUS=HepMC3_set_attribute_double(OUTID,scale_process,
     &  'event_scale'//char(0))
C...Note: no XS uncertainty
          HEPMC3STATUS=HepMC3_write_event(OUTID)
C-- This is for Rivet with HepMC3 interface          
         if (rvtOUT) then
         if (I.EQ.1) then 
         RIVETSTATUS=rivetinitfirstevent(OUTID)
         end if
         RIVETSTATUS=rivetrun(OUTID)
         end if
         HEPMC3STATUS=HepMC3_clear_event(OUTID)  
         endif
         
        
        if(N.le.0) write(6,*) "N=0  !!!!!"
111     Continue        
      Enddo
123   Continue
c
c Print out of generated events summary
c 
c      write(6,*) ' ilha ',ilha
      if(ilha.ne.0) then
         avgi = sumOfWeights/float(Nevent)*dabs(avgi)
         sd = sumOfWeights/float(Nevent)*dabs(sd)
c         avgi = sumOfWeights/float(Nevent)*avgi
c         sd = sumOfWeights/float(Nevent)*sd
      endif
c      write(6,*) 'sumofweights', sumOfWeights,float(Nevent),dabs(avgi)
c xsec in CASCADE in nv, hepmc need pb
CAV      xsec = avgi * 1000.
CAV      xsece = sd *1000.
CAV   WRITE(6,10000) xsec, xsece
CAV      call convhepmc(iev, isub, weighthepmc, xsec, xsece, idP1,idP2,xidP1,xidP2,scale_process,pdf1,pdf2 )
   
      if (hm3OUT) HEPMC3STATUS=HepMC3_delete_writer(OUTID)
      if (rvtOUT) then
      RIVETSTATUS=rivetdone(trim(rvtfilnam)//char(0))
      endif
   
      CALL CAEND(20)
c
10000 FORMAT(' final x-sections: sigma  = ',G10.3,' pb   +/- ',G10.3)
      sec = AVGI
      secrms = SD
      	
      if(ilha.ne.0) then          
         write(6,*) ' sum of weights ',sumOfWeights   
      endif   
      call cpu_time(finish)
      write(6,*) ' CASCADE running time = ',finish-start,' seconds.  Time/event =',(finish-start)/Nevent
      STOP

      END



