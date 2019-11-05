#include <iostream>
#include "HepMC/PythiaWrapper.h"
#include "HepMC/IO_HEPEVT.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/HEPEVT_Wrapper.h"
#include "PythiaHelper.h"

using namespace std;

extern "C" {
    HepMC::GenEvent* event_hepmc2=NULL;
    HepMC::IO_HEPEVT hepevent_hepmc2io;
    HepMC::IO_GenEvent* ascii_io=NULL;
    int  inithepmc_(char * outfile )
    {
        if (strlen(outfile)>0)
            ascii_io=new HepMC::IO_GenEvent(outfile ,std::ios::out);
        else
            ascii_io=new HepMC::IO_GenEvent("rapgap.hepmc" ,std::ios::out);
        return 0;
    }
    int getorig_(int &a);


    int convhepmc_(int & ievent, int & iproc, double & xsec, double & xsece,
                    int& flav1, int& flav2,double &  x1,double &  x2,double &  q2pdfeval,
                    double &  xf1mom, double & xf2mom ,int&  pdf1,int& pdf2
                   ) {


        hepevent_hepmc2io.set_trust_mothers_before_daughters( true );

        // pythia pyhepc routine convert common PYJETS in common HEPEVT
        call_pyhepc( 1 );
        if (event_hepmc2) delete event_hepmc2;
        event_hepmc2 = hepevent_hepmc2io.read_next_event();
        event_hepmc2->use_units(HepMC::Units::GEV, HepMC::Units::MM);
        event_hepmc2->set_event_number(ievent);
        event_hepmc2->set_signal_process_id(iproc);
        event_hepmc2->set_mpi(-1);
        event_hepmc2->set_event_scale(q2pdfeval);
        event_hepmc2->set_alphaQED(-1);
        event_hepmc2->set_alphaQCD(-1);
        std::vector<long> rs;
        rs.push_back(1);
        event_hepmc2->set_random_states(rs);

        //Set beams
        event_hepmc2->barcode_to_particle(1)->set_status(4);
        event_hepmc2->barcode_to_particle(2)->set_status(4);
//        for (HepMC::GenEvent::particle_iterator  p=event_hepmc2->particles_begin(); p!=event_hepmc2->particles_end(); p++)
//        (*p)->set_barcode_(10000+(*p)->barcode());
        //Set PDF info
        HepMC::PdfInfo pdf( flav1, flav2, x1, x2, q2pdfeval, xf1mom, xf2mom , pdf1,pdf2);
        event_hepmc2->set_pdf_info(pdf);


        //Fix problems with broken record
        //1) Detached FSR photons and stable electrons with end vertex
        bool created_proper_final_electron=false;
        for (HepMC::GenEvent::particle_iterator  p=event_hepmc2->particles_begin(); p!=event_hepmc2->particles_end(); p++)
        {
            //Stable particles with end vertex should not exist
            if ((*p)->end_vertex()&&(*p)->status()==1)
            {
                //In case that is a photon we know that FSR/ISR photons enter hadronisation. We detach their end vertices.
                if (std::abs((*p)->pdg_id())==22) (*p)->end_vertex()->remove_particle((*p));
                //In case that is an electron we have to do more manipulations.
                if (std::abs((*p)->pdg_id())!=11) continue;
                //Manipulations were done
                if (created_proper_final_electron) continue;
                //We assume to find detached photon
                if ((*p)->end_vertex()->particles_out_size()!=1) continue;
                HepMC::GenParticle* G=*((*p)->end_vertex()->particles_out_const_begin());
                HepMC::GenParticle* P=new HepMC::GenParticle();
                HepMC::FourVector XX1=(*p)->momentum();
                HepMC::FourVector XX2=G->momentum();
                HepMC::FourVector XX(XX1.px()+XX2.px(),XX1.py()+XX2.py(),XX1.pz()+XX2.pz(),XX1.t()+XX2.t());
                (*p)->set_status(11); //set status intermediate
                (*p)->set_momentum(XX); //set momentum which is a sum of final  state e+gamma
                //create proper final state electron.
                P->set_status(1);
                P->set_momentum(XX1);
                P->set_pdg_id((*p)->pdg_id());
                //add proper final state electron to end vertex of original electron.
                (*p)->end_vertex()->add_particle_out(P);
                created_proper_final_electron=true;// Doing that just once. Not sure if there are cases where nore than one is needed.
            }
        
        }
//        printf("conv %i\n",event_hepmc2->particles_size());

        //2) ISR photons in interaction vertex
        HepMC::GenParticle* EIN=event_hepmc2->barcode_to_particle(1);      //e Beam
        HepMC::GenVertex* VI=EIN->end_vertex();                   //Interaction vertex
        HepMC::GenParticle* G=0;                                  //ISR
        for (HepMC::GenEvent::particle_iterator  p=event_hepmc2->particles_begin(); p!=event_hepmc2->particles_end(); p++)
        {
            if ((*p)->production_vertex()!=VI) continue;
            int bc=-1;
            if (std::abs((*p)->pdg_id())==22) bc=(*p)->barcode();
            else continue; //Found photon in the interaction vertex
            int orig=getorig_(bc);                                               //Look in Pythia where it came from
            //printf("Found gamma in VI; bc=%i %i\n",bc, orig);
            if (orig!=1) continue;                                               //Check it is ISR
            G=*p;
        }

        //If ISR photon was found replace
        // e_beam -> Interaction_vertex -> ISRgamma
        // with
        //  e_beam -> ISRgamma e -> Interaction_vertex
        if (G) {
            HepMC::FourVector XX1=EIN->momentum();
            HepMC::FourVector XX2=G->momentum();
            HepMC::FourVector XX(XX1.px()-XX2.px(),XX1.py()-XX2.py(),XX1.pz()-XX2.pz(),XX1.t()-XX2.t());

            HepMC::GenParticle* ER= new HepMC::GenParticle(XX,EIN->pdg_id(),11);
            HepMC::GenVertex*  VR= new HepMC::GenVertex();


            event_hepmc2->add_vertex(VR);
            VI->remove_particle(EIN);
            VI->remove_particle(G);
            VR->add_particle_in(EIN);
            VR->add_particle_out(ER);
            VR->add_particle_out(G);
            VI->add_particle_in(ER);
        }

        event_hepmc2->weights().push_back(1.0);

        //      std::cout << " ievent " << ievent << " iproc " << iproc << " xsec " <<xsec<< std::endl;
        // set cross section information set_cross_sectio( xsec, xsec_err)
        const double xsecval = xsec;
        const double xsecerr = xsece ;
        HepMC::GenCrossSection cross;
        cross.set_cross_section( xsecval, xsecerr );
        event_hepmc2->set_cross_section( cross );
        // write the event out to the ascii file
        if (ascii_io)
        (*ascii_io) << event_hepmc2;
    return 0;
    }
}
