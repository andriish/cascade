// -*- C++ -*-
//
// This file is part of HepMC3
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#ifndef Pythia6_Pythia6ToHepMC2_H
#define Pythia6_Pythia6ToHepMC2_H
#ifdef  DUMMYPYTHIA6TOHEPMC2
extern "C" {

    int hepmc2_delete_writer_(const int & position)
    {
        return 0;
    }
    int hepmc2_convert_event_(const int & position)
    {
        return 0;
    }
    int hepmc2_write_event_(const int & position)
    {
        return 0;
    }
    int hepmc2_clear_event_(const int & position)
    {
        return 0;
    }
    int hepmc2_set_cross_section_(const int & position, const double& x,const double& xe, const int& n1,const int& n2)
    {
        return 0;
    }

    int hepmc2_set_pdf_info_(const int & position,const int& parton_id1, const int& parton_id2, const double& x1, const double& x2,
                             const double& scale_in, const double& xf1,const double& xf2,
                             const int& pdf_id1, const int& pdf_id2)
    {
        return 0;
    }
    int hepmc2_set_hepevt_address_(int* a)
    {
        return 0;
    }
    int hepmc2_set_attribute_int_(const int & position,const int & attval,const char* attname)
    {
        return 0;
    }
    int hepmc2_set_attribute_double_(const int & position,const double & attval,const char* attname)
    {
        return 0;
    }
    int hepmc2_new_writer_(const int & position,const int & mode,const char* ffilename)
    {
        return  0;
    }
    int hepmc2_new_weight_(const int & position, const char* name)
    {
        return 0;
    }
    int hepmc2_set_weight_by_index_(const int & position,const double& val, const int & pos)
    {
        return 0;
    }
    int hepmc2_set_weight_by_name_(const int & position,const double& val, const char* name)
    {
        return 0;
    }
}


#else
#include <iostream>
#include "HepMC/PythiaWrapper.h"
#include "HepMC/IO_HEPEVT.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/HEPEVT_Wrapper.h"
#include "PythiaHelper.h"

std::map<int,std::pair<HepMC::IO_GenEvent*,HepMC::GenEvent*> > hepmc2_gWriters;
HepMC::IO_HEPEVT  hepmc2_gHEPEVT;
HepMC::GenEvent* hepmc2_gWriters_get_event(const int & position)
    {
     return    hepmc2_gWriters[position].second;
     }
using namespace HepMC;
extern "C" {
    int getorig_(int &a);
    int hepmc2_delete_writer_(const int & position)
    {
        return 0;
    }
    int hepmc2_convert_event_(const int & position)
    {
        hepmc2_gHEPEVT.set_trust_mothers_before_daughters( true );

        // pythia pyhepc routine convert common PYJETS in common HEPEVT
        call_pyhepc( 1 );
        GenEvent* event_hepmc2=hepmc2_gWriters[position].second;
        if (event_hepmc2) delete event_hepmc2;
        
        event_hepmc2= hepmc2_gHEPEVT.read_next_event();
        hepmc2_gWriters[position].second  = event_hepmc2;
        event_hepmc2->use_units(HepMC::Units::GEV, HepMC::Units::MM);


        //Set beams
        event_hepmc2->barcode_to_particle(1)->set_status(4);
        event_hepmc2->barcode_to_particle(2)->set_status(4);

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

        return 0;
    }
    int hepmc2_write_event_(const int & position)
    {
    (*hepmc2_gWriters[position].first) << hepmc2_gWriters[position].second;
        return 0;
    }
    int hepmc2_clear_event_(const int & position)
    {
        hepmc2_gWriters[position].second->clear();
        return 0;
    }
    int hepmc2_set_cross_section_(const int & position, const double& x,const double& xe, const int& n1,const int& n2)
    {
        const double xsecval = x;
        const double xsecerr = xe ;
        HepMC::GenCrossSection cross;
        cross.set_cross_section( xsecval, xsecerr );
        hepmc2_gWriters[position].second->set_cross_section( cross );
        return 0;
    }

    int hepmc2_set_pdf_info_(const int & position,const int& parton_id1, const int& parton_id2, const double& x1, const double& x2,
                      const double& scale_in, const double& xf1,const double& xf2,
                      const int& pdf_id1, const int& pdf_id2)
    {
        HepMC::PdfInfo pdf( parton_id1, parton_id2,x1,x2,scale_in,xf1,xf2,pdf_id1, pdf_id2);
        hepmc2_gWriters[position].second->set_pdf_info(pdf);
        return 0;
    }

    int hepmc2_set_hepevt_address_(int* a)
    {
            return 0;
    }
    int hepmc2_set_attribute_int_(const int & position,const int & attval,const char* attname)
    {
		std::string sta(attname);
		GenEvent* event_hepmc2=hepmc2_gWriters[position].second;
		if (sta==std::string("mpi")) { event_hepmc2->set_mpi(attval); return 0;}
		if (sta==std::string("signal_process_id")) { event_hepmc2->set_signal_process_id(attval); return 0;}
		//        std::vector<long> rs;
//        rs.push_back(1);
//        event_hepmc2->set_random_states(rs);
        return 0;
    }
    int hepmc2_set_attribute_double_(const int & position,const double & attval,const char* attname)
    {
		std::string sta(attname);
		GenEvent* event_hepmc2=hepmc2_gWriters[position].second;
		if (sta==std::string("alphaQED")) { event_hepmc2->set_alphaQED(attval); return 0;}
		if (sta==std::string("alphaQCD")) { event_hepmc2->set_alphaQCD(attval); return 0;}
		if (sta==std::string("event_scale")) { event_hepmc2->set_event_scale(attval); return 0;}
        return 0;
    }

    int hepmc2_new_writer_(const int & position,const int & mode,const char* ffilename)
    {
		
		        int r_position=position;
        if (r_position==0)
        {
            if (hepmc2_gWriters.size()==0) r_position=1;
            if (hepmc2_gWriters.size()!=0) r_position=hepmc2_gWriters.rend()->first+1;
        }
        if (hepmc2_gWriters.find(r_position)!=hepmc2_gWriters.end()) {
            printf("Error in %s: Writer at position %i already exists\n",__FUNCTION__,r_position);
            exit(1);
        }
		
		hepmc2_gWriters[r_position]=std::pair<IO_GenEvent*,GenEvent*>( new HepMC::IO_GenEvent(ffilename ,std::ios::out),new GenEvent(Units::GEV,Units::MM));
        return  r_position;
    }
    int hepmc2_new_weight_(const int & position, const char* name)
    {
        return 0;
    }
    int hepmc2_set_weight_by_index_(const int & position,const double& val, const int & pos)
    {
        hepmc2_gWriters[position].second->weights()[pos]=val;
        return 0;
    }
    int hepmc2_set_weight_by_name_(const int & position,const double& val, const char* name)
    {
        return 0;
    }    
}
#endif
#endif
