      function getorig(a)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer a,getorig
      getorig=K(a,3)
      end 

C      function getq2q()
C#include "rgpara1.inc"
C      double precision getq2q
C      getq2q=Q2Q
C      end 