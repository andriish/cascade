#include <iostream>
#include <set>
#include <vector>
#include <string>
//Preparation for the future Rivet with HepMC3 
#ifdef USE_RIVET_HEPMC3
#include "HepMC3/GenEvent.h"
#else
#include "HepMC/GenEvent.h"
#endif
#include "Rivet/Rivet.hh"

extern "C" {

#ifdef USE_RIVET_HEPMC3
    int rivetinterfaceversion_()
    {
        return 3;
    }   
    /**  HepMC event to reads from*/
    HepMC3::GenEvent* event=NULL;    
    extern HepMC3::GenEvent* gWriters_get_event(const int & position);

#else
    int rivetinterfaceversion_()
    {
        return 2;
    }    
    /**  HepMC event to reads from*/
    HepMC::GenEvent* event=NULL;    
    extern HepMC::GenEvent* event_hepmc2;
#endif 
    /**The name of the file where the histograms are dumped.*/
    std::string filename;
    /** Analyses with optional analysis parameters.*/
    std::set<std::string> analyses;
    /**The Rivet object.*/
    Rivet::AnalysisHandler * rivet=NULL;
    /** Run name*/
    std::string rname;
    /**Ignore beams flag.*/
    bool igBeam=true;

    int rivetinit_(char* rname1) {
        if ( rivet ) return 0;
        rname=std::string(rname1);
        rivet = new Rivet::AnalysisHandler(rname);
        rivet->setIgnoreBeams(igBeam);
        Rivet::addAnalysisLibPath(".");
        for (std::set<std::string>::iterator it = analyses.begin();
                it != analyses.end(); ++it) {
            rivet->addAnalysis(*it);
        }
        return 0;
    }
    int rivetinitfirstevent_(const int &  id)
    {
#ifdef USE_RIVET_HEPMC3
        event=gWriters_get_event(id);
#else
        event=event_hepmc2;
#endif
        rivet->init(*event_hepmc2);
        return 0;
    }
    int rivetrun_(const int &  id) {
#ifdef USE_RIVET_HEPMC3
        event=gWriters_get_event(id);
#else
        event=event_hepmc2;
#endif
//         printf("rivet %i\n",event->particles_size());
        if (!event) { puts("Something is wrong with event!"); return 1;}
        
        if (!event->particles_size()) { puts("Something is wrong with particles!"); return 2;}
        if (!event->cross_section()) { puts("Something is wrong with cross-section!"); return 3;}
        rivet->analyze(*event_hepmc2);
        return 0;
    }
    int rivetadd_(char* ana)
    {
        analyses.insert(std::string(ana));
        return   analyses.size();
    }
    int rivetdone_(char* filename1) {
        if ( !rivet ) return 0;
        filename=std::string(filename1);
        rivet->finalize();
        rivet->writeData(filename);
        delete rivet;
        rivet = NULL;
        return 0;
    }
}
