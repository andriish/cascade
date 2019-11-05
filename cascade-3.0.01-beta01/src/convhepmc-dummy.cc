extern "C" {
    int ncount_hepmc2=0;
    void convhepmc_(int & ievent, int & iproc, double & xsec, double & xsece,
                    int& flav1, int& flav2,double &  x1,double &  x2,double &  q2pdfeval,
                    double &  xf1mom, double & xf2mom ,int&  pdf1,int& pdf2
                   )
    {
        if ( ncount_hepmc2 < 10) {
            printf(" RAPGAP: dummy version of convhepmc (HEPMC2) is used\n" );
            ++ncount_hepmc2;
        }
    }
}
