      SUBROUTINE casini

      IMPLICIT NONE
#include "casparas.inc"
#include "caspara5.inc"
#include "casparam.inc"
#include "casluco.inc"
#include "casludat1.inc"
#include "casvalues.inc"
#include "casshwr.inc"
#include "casprkt.inc"
#include "casscal.inc"
#include "casglu.inc"
#include "casmult.inc"
#include "casherup.inc"
#include "caspsipol.inc"
#include "caspar.inc"
#include "caslhe.inc"

      Integer Iremn
      common/casflags/Iremn	
      
      character*512 PDFPATH
      common/CASPDF/PDFPATH

      CHARACTER CHIN*100

	Integer IDEBUG
      COMMON/ INTERN/IDEBUG
      Double Precision ACC1,ACC2
	Integer IINT,NCB
      COMMON /INTEGR/ ACC1,ACC2,IINT,NCB
      Integer IGLU
	Common/CAGLUON/Iglu
      DOUBLE PRECISION SIN2W,XMW2
      COMMON/ELWEAK/SIN2W,XMW2
	DOUBLE PRECISION XMW
	Double Precision Qscal
	Integer Iqqbar
	Common/GLUDAT2/Qscal,Iqqbar
      Double Precision PARP,PARI
      Integer MSTP,MSTI
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      Integer I
	Double Precision  BB,ktqt,Qcdlam,scal,qtd,ktd
      Common /splitting/ BB,ktqt,Qcdlam,scal,qtd,ktd
      CALL BCASIN
c select unintegrated gluon : 	
c IGLU = 1: CCFM old set JS2001
c IGLU =  1001: CCFM J2003 set 1 
c IGLU =  1002: CCFM J2003 set 2 
c IGLU =  1003: CCFM J2003 set 3 
c IGLU =  1010: CCFM set A0
c IGLU =  1011: CCFM set A0+
c IGLU =  1012: CCFM set A0-
c IGLU =  1013: CCFM set A1
c IGLU =  1020: CCFM set B0
c IGLU =  1021: CCFM set B0+
c IGLU =  1022: CCFM set B0-
c IGLU =  1023: CCFM set B1
c
c IGLU=2 --> derivative GRV
c IGLU=3 --> Bluemlein
c IGLU=4 --> KMS
c IGLU=5 --> GBW
c IGLU=6 --> KMR
c IGLU=6 --> RS
      call getenv('PDFPATH',pdfpath)
      if(pdfpath.ne.'') then
         call getenv('PDFPATH',pdfpath)
      else
         PDFPATH='../share/'
      endif
      IGLU = 1201
      IDEBUG = 0
c collinear gluon(for tests...)
      call setlhaparm('LHAPDF')
c      call setlhaparm('SILENT')
      call setlhaparm('LOWKEY')
      ICOLGLU=10042 	
c select integration procedure
      NCB = 20000
      ACC1 = 1.
      ACC2 = 0.5
C FINAL STATE PARTON SHOWER
      IFPS = 3
      ILEPTO = 1
      IFINAL = 1
c     scale factor for final state PS
      SCALFAF = 1.      
C LEPIN KE = 11
      KE = 11
	Ires(1) = 0
      KGL = 21
      KPH=22
      KEB = KPH
      KPA = 1
      NFLAV = 5
      IHF = 0
	
C PPIN KP = 2212
      KP = 2212
	Ires(2) = 1
C select interaction type
C INTER = 0 --> photon exchange
      INTER = 0
C select BGF in semihard approach a la Catani et al. ISEMIH = 1
      ISEMIH = 1
C RUNNING ALPHA_S IRUNA = 1
      IRUNA=1
C RUNNING ALPHA_EM IRUNAEM = 1
      IRUNAEM=0
C SELECT SCALE FOR STRUCTURE FKT AND ALPHA_S
C IQ2 =1 SCALE Q2 = MASS**2
C IQ2 =2 SCALE Q2 = SHAT
C IQ2 =3 SCALE Q2 = MASS**2 + PT**2
C IQ2 =4 SCALE Q2 = Q2
C IQ2 =5 SCALE Q2 = Q2 + pt**2
      IQ2=3
C SELECT PROCESS
      IPRO= 10
	PT2CUT(2) = 0
	PT2CUT(3) = 0
	PT2CUT(10) = 0
	PT2CUT(11) = 0
	PT2CUT(500) = 0
	PT2CUT(503) = 0
	PT2CUT(504) = 0
	PT2CUT(505) = 0
	PT2CUT(510) = 0
	PT2CUT(513) = 0
	PT2CUT(514) = 0
	PT2CUT(515) = 0
C Minimum Q^2 of electron to be generated
      QMI = 5.d0
C Maximum Q^2 of electron to be generated
      QMA = 10D8
C Minimum y of electron to be generated
      YMI=0.0d0
C Minimum y of electron to be generated
      YMA=1.0d0
C Maximium theta angle of scattered electron
      THEMA = 180.0D0
C Minimum  theta angle of scattered electron
      THEMI =   0.0D0
C ELECTRON MOMENTUM
      PLEPIN =-30.
C PROTON MOMENTUM
      PPIN   = 820.
C PERFORM FRAGMENTATION NFRAG=0/1
      NFRAG = 1
c printout of intermediate/virtual partons during initial state cascade
      IPST = 0
c calculate with jpsi polaristion
      ipsipol = 0
c 2s/3s state for vectormeson
      i23s = 0      
c safety for fragmentation
      PYPAR(12) = 1.5
      MSTU(12) = 0
      CHIN='PMAS(4,1)=1.5'
      CALL PYGIVE(CHIN)
c      CALL PYGIVE('PMAS(1,1)=0.450')
c      CALL PYGIVE('PMAS(2,1)=0.450')
      IHFLA = 3
c

	ICCFM = 1
	
	IORDER = 1
      BB = 0.5
c primordial kt
      PRKT1 = 0d0
      PRKT2 = 0d0
c width of Gaussian distribution in transverse momentum
c when a non-trivial target remnant is split into two particles
c (PARL(14) in LEPTO )
      PRKTREM = 0.35d0
c energy sharing of remnant:
      IRspl = 4	
	SCALFA = 1d0	
c multiple scattering
      Mult = 0	
C SET ELECTROWEAK PARAMETERS
      SIN2W = 0.23D0
      XMW = 80.D0
      XMW2 = XMW**2
c parameters for scales
      Qscal = 1d0
	Iqqbar = 1

c use les houches accord for frag and final state ps 
      ILHA = 0	
      
c      MSTP(128) =1
c         write(6,*) ' CASINI : Printing PYTHIA params '
c         call pylist(13)
c      write(6,*) ' mstp(120-130) ',(MSTP(I),I=120,130)
c      write(6,*) ' mstp(128) ',MSTP(128)

c      hepmcout='cascade.hepmcout'//char(0)
      hepmcout='/dev/null'//char(0)
      lhefile = 'cascade.lhe'//char(0)
      Ilhe = 0
      iTMDw = 0
      Iscale = 1
      lheWid = 0

      Iremn = 1 ! keep remnant for hadronisation
      
      p6frag = 1 ! use Pythia6 for fragmentation
      
      IordasPS = 1 
      
      scale_process = -99999.
      
      WRITE(6,*) ' casini: initialisation of CASCADE '

      RETURN
      END
C...  END INITIALIZE
