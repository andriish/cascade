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

using namespace std;

extern "C" {
    extern struct {
      char  hepmcout[];
    } cahepmcout_;
}
#define cahepmcout cahepmcout_

extern "C" {   
    HepMC::GenEvent* event_hepmc2=NULL;
    HepMC::IO_HEPEVT hepevent_hepmc2io;
    HepMC::IO_GenEvent* ascii_io=NULL;

    int convhepmc_(int & ievent, int & iproc, double & weight, double & xsec, double & xsece, int & id1pdf, int & id2pdf, double & x1pdf, double & x2pdf, double & QFac, double & pdf1, double & pdf2 ){
      
       // with the static command, called only once
      static char * outfile = cahepmcout.hepmcout+'\0' ;
      ascii_io= new  HepMC::IO_GenEvent(outfile ,std::ios::out); 
    

      hepevent_hepmc2io.set_trust_mothers_before_daughters( true );
   // pythia pyhepc routine convert common PYJETS in common HEPEVT
      call_pyhepc( 1 );
      HepMC::GenEvent* event_hepmc2 = hepevtio.read_next_event();
      event_hepmc2->use_units(HepMC::Units::GEV, HepMC::Units::MM);
      event_hepmc2->set_event_number(ievent);
      event_hepmc2->set_signal_process_id(iproc);
      event_hepmc2->set_mpi(-1);
//      event_hepmc2->set_event_scale(q2pdfeval);
      event_hepmc2->set_alphaQED(-1);
      event_hepmc2->set_alphaQCD(-1);
      //Set beams
      event_hepmc2->barcode_to_particle(1)->set_status(4);
      event_hepmc2->barcode_to_particle(2)->set_status(4);
      // Store PDF information.
      HepMC::PdfInfo pdf( id1pdf, id2pdf, x1pdf, x2pdf, QFac, pdf1, pdf2, 230, 230);
      event_hepmc2->set_pdf_info(pdf);      // write the event out to the ascii file

        event_hepmc2->weights().push_back(weight);

        // set cross section information set_cross_sectio( xsec, xsec_err)
        const double xsecval = xsec;
        const double xsecerr = xsece;
        HepMC::GenCrossSection cross;
        cross.set_cross_section( xsecval, xsecerr );
        event_hepmc2->set_cross_section( cross );
        // write the event out to the ascii file
        if (ascii_io)
        (*ascii_io) << event_hepmc2;
    return 0;
    }
}
