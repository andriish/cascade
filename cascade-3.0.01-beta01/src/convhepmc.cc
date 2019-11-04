//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#include <cmath>
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

extern "C" {
    extern struct {
      char  hepmcout[];
    } cahepmcout_;
}
#define cahepmcout cahepmcout_

extern "C" {   
    // Instantiate an IO strategy for reading from HEPEVT.
       HepMC::IO_HEPEVT hepevtio;
       
       int ncount = 0;   
       HepMC::GenCrossSection cross;

//void convhepmc_(int & ievent, int & iproc, double & xsec){
void convhepmc_(int & ievent, int & iproc, double & weight, double & xsec, double & xsece, int & id1pdf, int & id2pdf, double & x1pdf, double & x2pdf, double & QFac, double & pdf1, double & pdf2 ){
      
       // with the static command, called only once
       static char * outfile = cahepmcout.hepmcout+'\0' ;
       static HepMC::IO_GenEvent ascii_io(outfile ,std::ios::out); 

      
//       cout << " convhepmc=" << outfile << "$end"<<endl;
      if ( ncount < 10) {
        if (outfile!=NULL) { cout << " convhepmc: filename = " <<  outfile << endl;}
          else { cout << "  convhepmc: NO filename set " <<endl; 
          return ;}
          ++ncount;
      }
      

   // pythia pyhepc routine convert common PYJETS in common HEPEVT
	call_pyhepc( 1 );
      //HepMC::HEPEVT_Wrapper::print_hepevt();
      HepMC::GenEvent* evt = hepevtio.read_next_event();
      //HepMC::IO_HEPEVT::print_inconsistency_errors();
      
	// HepMC::HEPEVT_Wrapper::check_hepevt_consistency();
      // HepMC::IO_HEPEVT::set_trust_mothers_before_daughters( true );
      //   from version 2.06.09 on:      
      //      evt->define_units( HepMC::Units::GEV, HepMC::Units::MM );
      evt->use_units(HepMC::Units::GEV, HepMC::Units::MM);
      evt->set_event_number(ievent);
      evt->set_signal_process_id(iproc);
      evt->weights().push_back(weight);
      //      std::cout << " ievent " << ievent << " iproc " << iproc << " xsec " <<xsec<< std::endl;
      // set cross section information set_cross_sectio( xsec, xsec_err)
      const double xsecval = xsec;
      const double xsecerr = xsece ;
      cross.set_cross_section( xsecval, xsecerr );
	evt->set_cross_section( cross );
      // Store PDF information.
      // evt->set_pdf_info(id1pdf, id2pdf, x1pdf, x2pdf, QFac, pdf1, pdf2 );
      HepMC::PdfInfo pdf( id1pdf, id2pdf, x1pdf, x2pdf, QFac, pdf1, pdf2, 230, 230);
      evt->set_pdf_info(pdf);      // write the event out to the ascii file
	ascii_io << evt;

    delete evt;
   }
}
