#include <iostream>
#include <set>
#include <vector>
#include <string>
#include "HepMC/GenEvent.h"
#include "Rivet/Rivet.hh"

extern "C" {
    /**  HepMC event to reads from*/
    extern HepMC::GenEvent* event_hepmc2;
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
    int rivetinitfirstevent_()
    {
        rivet->init(*event_hepmc2);
        return 0;
    }
    int rivetrun_() {

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
