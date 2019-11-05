      subroutine steer_new

      implicit none
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
#include "casparam.inc"
      Integer Iremn
      common/casflags/Iremn	

      character*512 PDFPATH
      common/CASPDF/PDFPATH
      Integer ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON /CASCOL/ ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG

      Integer IINT,NCB
      Double Precision ACC1,ACC2
      COMMON /INTEGR/ ACC1,ACC2,IINT,NCB
      Integer IGLU
	Common/CAGLUON/Iglu



      Integer NrEvents,Process_Id,Hadronisation,
     & SpaceShower,SpaceShowerOrderAlphas,
     & TimeShower,ScaleTimeShower,
     & PartonEvolution,
     & PartonDensity,
     & EnergyShareRemnant,Remnant,
     & lheHasOnShellPartons,lheReweightTMD,lheScale,lheWeightId,
     & pythia6
      Double Precision ScaleFactorFinalShower
      character *512 TMDDensityPath,lheInput
C
      character *132 hm2filnam
      Common/steerhm2/hm2filnam
      LOGICAL         hm2OUT
      COMMON /Qhm2LO/ hm2OUT

      character *132 hm3filnam
      integer hm3filfor
      Common/steerhm3/hm3filnam,hm3filfor
      LOGICAL         hm3OUT
      COMMON /Qhm3LO/ hm3OUT

      character *132 rvtfilnam
      character *132 rvtrun
      character *132 rvtana
      Integer Nrvtana
      Common/steerrvt/Nrvtana,rvtfilnam,rvtrun,rvtana(50)
      LOGICAL         rvtOUT,rvtOUT1,rvtOUT2
      COMMON /QrvtLO/ rvtOUT
C      
      
      
      namelist/CASCADE_input/
     & NrEvents,Process_Id,Hadronisation,
     & SpaceShower,SpaceShowerOrderAlphas,
     & TimeShower,ScaleTimeShower,ScaleFactorFinalShower,
     & PartonEvolution,
     & PartonDensity,
     & EnergyShareRemnant,Remnant,
     & lheHasOnShellPartons,lheReweightTMD,lheWeightId,lheScale,
     & pythia6,
     & TMDDensityPath,lheInput
      
      NrEvents=Nevent
      Process_Id  =Ipro 
      Hadronisation=Nfrag
      if(IFPS.eq.1.or.IFPS.eq.3) Spaceshower=1
      if(IFPS.eq.2.or.IFPS.eq.3) Timeshower=1
      IFPS = 0
      SpaceShowerOrderAlphas = IordasPS
c scalfa = scale factor for alphas
c scalfaf = scale factor for final state scale
c ifin = defintion of scale for final state shower
      ScaleTimeShower=ifinal
      PartonEvolution=Iccfm
      EnergyShareRemnant=IRspl
      Remnant=IRemn
      PartonDensity =IGLU
      lheHasOnShellPartons=ILHE
      lheReweightTMD=iTMDw
      lheWeightId=lheWid
      lheScale=Iscale
      ScaleFactorFinalShower=Scalfaf
      TMDDensityPath=PDFPATH
      lheInput=lhefile
      pythia6=p6frag

c      read (5,NML=CASCADE_input,END=141,ERR=42)
      read (5,NML=CASCADE_input,END=141)
      WRITE ( *, CASCADE_input ) 
      
      Nevent = NrEvents
      Ipro = Process_Id      
      Nfrag = Hadronisation     
      if(SpaceShower.eq.1) IFPS = 1
      if(TimeShower.eq.1) IFPS = 2
      if(SpaceShower.eq.1.and.TimeShower.eq.1) IFPS = 3
      Ifinal = ScaleTimeShower
      Iccfm = PartonEvolution
      IRspl = EnergyShareRemnant
      IRemn = Remnant
      IGLU = PartonDensity    
      if(iglu.gt. 102100 .and.Iglu.le.102290 ) IordasPS = 2
      IordasPS = SpaceShowerOrderAlphas
      ILHE = lheHasOnShellPartons
      iTMDw = lheReweightTMD
      lheWid = lheWeightId
      Iscale = lheScale
      p6frag = pythia6
      
      Scalfaf = ScaleFactorFinalShower
      PDFPATH = trim(TMDDensityPath)
      lhefile = lheInput
c      Collfile = 
      
      
      close (5)
      return
      
 141  continue
      close (5)
      
      return
      
 42   continue
      print '(''Error reading namelist &CASCADE_input, STOP'')'
      write( *, CASCADE_input ) 
      STOP
      end       
      
