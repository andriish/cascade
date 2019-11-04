//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#include <cmath>
#include <iostream>

using namespace std;

extern "C" {   
int ncount ;
void convhepmc_(int & ievent, int & iproc, double & xsec, double & xsece ){
      if ( ncount < 10) {
          cout << " CASCADE: dummy version of convhepmc is used "<< endl;    
          ++ncount;
      }
   }   
}
