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
include "HepMCSteering.inc"
C      
      
      LOGICAL HepMC2OUT,  HepMC3OUT,  RivetOUT
      character *132 HepMC2File,  HepMC3File  
      character *132 RivetFile,RivetAnalyses(50), RivetRun
      integer HepMC3Format
      namelist/CASCADE_input/
     & NrEvents,Process_Id,Hadronisation,
     & SpaceShower,SpaceShowerOrderAlphas,
     & TimeShower,ScaleTimeShower,ScaleFactorFinalShower,
     & PartonEvolution,
     & PartonDensity,
     & EnergyShareRemnant,Remnant,
     & lheHasOnShellPartons,lheReweightTMD,lheWeightId,lheScale,
     & pythia6,
     & TMDDensityPath,lheInput,
     & HepMC2File,  HepMC3File,  RivetFile,
     & HepMC2OUT,  HepMC3OUT,  RivetOUT,
     & HepMC3Format,
     & RivetAnalyses, RivetRun
       

      
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
      
       hm2OUT=HepMC2OUT
       hm3OUT=HepMC3OUT
       hm2filnam=HepMC2File
       hm3filnam=HepMC3File
       hm3filfor=HepMC3Format
       rvtana=RivetAnalyses
       rvtrun=RivetRun
       rvtOUT=RivetOUT
       rvtfilnam=RivetFile
       Nrvtana=size(rvtana)
       
      
      
      
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
      
