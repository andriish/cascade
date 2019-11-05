#include <stdio.h>
extern "C" {
   int ncount_hepmc3=0 ;
   int hepmc3_delete_writer_(const int &)
    {
        return 0;
        
    }
    int hepmc3_convert_event_(const int &)
    {
        return 0;
    }
    int hepmc3_write_event_(const int &)
    {
        if ( ncount_hepmc3 < 10) {
            printf(" CASCADE: dummy version of write_event (HEPMC3) is used\n" );
            ++ncount_hepmc3;
        }
        return 0;
    }
    int hepmc3_clear_event_(const int &)
    {
        return 0;
    }
    int hepmc3_set_cross_section_(const int & position, const double& x,const double& xe, const int& n1,const int& n2)
    {
        return 0;
    }
    int hepmc3_set_hepevt_address_(int*)
    {
        return 0;
    }
    int hepmc3_set_attribute_int_(const int &,const int &,const char*, size_t)
    {
        return 0;
    }
    int hepmc3_set_attribute_double_(const int &,const double &,const char*, size_t)
    {
        return 0;
    }
    int hepmc3_new_writer_(const int & position,const int & mode,const char* ffilename, size_t len)
    {
        return 0;
    }
    int hepmc3_set_pdf_info_(const int & position,const int& parton_id1, const int& parton_id2, const double& x1, const double& x2,
                      const double& scale_in, const double& xf1,const double& xf2,
                      const int& pdf_id1, const int& pdf_id2)
    {
        return 0;
    }
    
}
