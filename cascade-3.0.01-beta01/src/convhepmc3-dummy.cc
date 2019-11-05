extern "C" {
   int ncount_hepmc3=0 ;
   int delete_writer_(const int & position)
    {
        return 0;
    }
    int convert_event_(const int & position)
    {
        return 0;
    }
    int write_event_(const int & position)
    {
        if ( ncount_hepmc3 < 10) {
            printf(" CASCADE: dummy version of write_event (HEPMC3) is used\n" );
            ++ncount_hepmc3;
        }
        return 0;
    }
    int clear_event_(const int & position)
    {
        return 0;
    }
    int set_cross_section_(const int & position, const double& x,const double& xe, const int& n1,const int& n2)
    {
        return 0;
    }
    int set_hepevt_address_(int* a)
    {
        return 0;
    }
    int set_attribute_int_(const int & position,const int & attval,const char* attname, size_t len)
    {
        return 0;
    }
    int set_attribute_double_(const int & position,const double & attval,const char* attname, size_t len)
    {
        return 0;
    }
    int new_writer_(const int & position,const int & mode,const char* ffilename, size_t len)
    {
        return 0;
    }
}
