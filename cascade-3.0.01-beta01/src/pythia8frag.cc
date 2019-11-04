#include <stdio.h>
#include <string.h>
#include <cstring>
#include <iostream>
#include "HepMC/PythiaWrapper.h"
#include "HepMC/IO_HEPEVT.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/HEPEVT_Wrapper.h"
#include "PythiaHelper.h"
#include "boost/algorithm/string/trim.hpp"
//#include "HepMC/WeightContainer.h"

using namespace std;


#include <string>
#include <iostream>

#include "Pythia8/Pythia.h"
using namespace Pythia8;



extern "C" {   

int ncountp8 =0 ;

void pythia8frag_(){
    
   cout << " pythia8frag=" << endl;
  // Generator; shorthand for event and particleData.
  static Pythia pythia;
  static Event& event      = pythia.event;
  static ParticleData& pdt = pythia.particleData;
  // Set true to also see space-time information in event listings.
  static int type=0;
  static bool showScaleAndVertex = (type == 0) ? true : false;

  // Initialize.
  pythia.init();
  // Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
  //pythia.readString("ProcessLevel:all = off");
  // Reset event record to allow for new event.
  //double ee= 100;

  //double mm = pdt.mSel(id);
  //double pp = sqrtpos(ee*ee - mm*mm);
  event.reset();
  
  cout<< " PYJETS "<< pyjets.n<< endl;
  //event.list();
   
  int icol_index = 0 ;
  int Istr = 0 ;
  int icol1 = 0 ;
  int icol2 = 0;
  
   
  for ( int i=0; i<= pyjets.n - 1 ; i++ ){
//     cout << " filling event record " << i << " " <<pyjets.n << endl;
  int iorig = pyjets.k[2][i] ;
  int iflavor = pyjets.k[1][i] ;
  int istat = pyjets.k[0][i] ;
  int istat_old = pyjets.k[0][i-1] ;
  //cout << " i = " << i << " istat = " << istat << " istat_old " << istat_old << " iflavor = " << iflavor << " iorig = "<< iorig << endl;
//  cout << " i " << i << " k4 " << pyjets.k[3][i] << " k5 " << pyjets.k[4][i] << endl;
  int istat_p8 =istat ;
  
  if (istat == 21 ) { istat_p8 = -12 ;}
  else if(istat == 2 && Istr == 0 ) { 
         istat_p8 = 71;
         Istr = 1 ;
         icol_index = icol_index +1  ;
         //cout << " icol " << icol_index << "  " << iflavor << endl;
         if( iflavor != 21 && iflavor > 0 ) { icol1 = 500+icol_index ;}
         else if( iflavor != 21 && iflavor < 0 ) { icol2 = 500+icol_index ;}
         } 
         //cout << "  gluon"  << icol1 << " " << icol2 << endl; 
  else if (istat == 2 && icol1 > 0 && icol2 == 0 ){ 
            istat_p8 = 71;
            icol2 = 500 + icol_index;  
            icol_index= icol_index +1 ;
            icol1 = 500 + icol_index; }
  else if (istat == 2 && icol1 == 0 && icol2 > 0) { 
            istat_p8 = 71;
            icol1 = 500 + icol_index ; 
            icol_index= icol_index +1 ;
            icol2 = 500 + icol_index; }
         //else  {icol1 = icol2 ; }     
         //cout << "  parton "  << iflavor << " icol " << icol1 << " " << icol2 << endl; 
  else if (istat == 2 && icol1 != 0 && icol2 != 0) { 
            istat_p8 = 71;
            icol1 = 500 + icol_index ; 
            icol_index= icol_index +1 ;
            icol2 = 500 + icol_index; 
            icol_index= icol_index +1 ;
  }
  else if( istat == 2 ) {
     istat_p8 = 71;
     Istr = 0;  }
  else if(istat == 1 && istat_old == 2 ) {
      istat_p8 = 71; 
      if (Istr == 0) {
         Istr = 1;
         }
      }  
  else if(istat == 1 ) {istat_p8 = 23; }  
  else if(istat == 11 ) {istat_p8 = -21; }
  else if(istat == 13 ) {istat_p8 = -21; }
  else if(istat == -1) {istat_p8 = -299; }  // treat old continuation lines
 // append( id, status, col, acol, px, py, pz, e, m)  
  event.append( iflavor, istat_p8, icol1, icol2 ,pyjets.p[1][i], pyjets.p[2][i],  pyjets.p[3][i], pyjets.p[4][i], pyjets.p[5][i]);
 // cout << "  parton "  << iflavor << " icol " << icol1 << " " << icol2 << endl; 
  
  }
  //event.list(showScaleAndVertex);
  event.list();

   ++ncountp8  ; 
      

   }
}
