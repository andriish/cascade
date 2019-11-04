
      SUBROUTINE XSECMM(SUMSUM)
C*******************************************************************C
C      Partonic differential cross section                          C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      COMMON/MOMEN1/G1(4),A1(4),G2(4),A2(4),PJ1(4),PJ2(4)
      COMMON/SPINJ1/SPJ1(4,4),SPS1(4,4),SPJ2(4,4),SPS2(4,4)
      COMMON/AMPJJ1/AMP(4,4)
      SAVE
      SUMSUM=0.D0
      DO 1 m=1,4
      DO 1 n=1,4
      AMP(m,n)=0.D0
    1 CONTINUE
C
C... Polarization-summed J/psi spin density matrices
      CALL GAUGEJ1
C... Feynman diagrams with normal ordering of Gluons
      CALL SCALAR1(G1,A1, G2,A2, PJ1,PJ2)
      CALL QEDSUM1(G1,A1, G2,A2, PJ1,PJ2)
      CALL QCDSUM1(G1,A1, G2,A2, PJ1,PJ2)
C    Feynman diagrams with the permutation of Gluons
      CALL SCALAR1(G2,A2, G1,A1, PJ1,PJ2)
      CALL QEDSUM1(G2,A2, G1,A1, PJ1,PJ2)
      CALL QCDSUM1(G2,A2, G1,A1, PJ1,PJ2)
C
C... Convolute the production amplitudes
      DO 10 m1=1,4
      DO 10 m2=1,4
      DO 10 n1=1,4
      DO 10 n2=1,4
      SUMSUM = SUMSUM 
     .         + AMP(m1,m2)*AMP(n1,n2)*SPJ1(m1,n1)*SPJ2(m2,n2)
   10 CONTINUE
      RETURN
      END

      SUBROUTINE GAUGEJ1
C*******************************************************************C
C      J/psi polarization vectors and spin density matrix           C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/SPINJ1/SPJ1(4,4),SPS1(4,4),SPJ2(4,4),SPS2(4,4)
      COMMON/MOMEN1/G1(4),A1(4),G2(4),A2(4),PJ1(4),PJ2(4)
      COMMON/QUARK1/XC,XC2,XB,XB2,XJ,XJ2,CHC,CHB
      COMMON/GMUNU1/DF(4,4),DC(4)
      SAVE
C... Invariant definition of Spin density matrix
      DO 20 m=1,4
      DO 20 n=1,4
      SPJ1(m,n)= -(DF(m,n) - PJ1(m)*PJ1(n)/XJ2) *DC(m)*DC(n)
      SPJ2(m,n)= -(DF(m,n) - PJ2(m)*PJ2(n)/XJ2) *DC(m)*DC(n)
   20 CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION AF(X,Y,Z)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      AF = DSQRT(X*X + Y*Y + Z*Z - 2.*X*Y - 2.*X*Z - 2.*Y*Z)
      RETURN
      END

      DOUBLE PRECISION FUNCTION AF2(X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      AF2 = X*X + Y*Y + Z*Z - 2.*X*Y - 2.*X*Z - 2.*Y*Z
      RETURN
      END

      DOUBLE PRECISION FUNCTION GF(X,Y,Z,U,V,W)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      GF =          X*Z*W + X*U*V + Y*Z*V + Y*U*W
     . -X*Y*(Z+U+V+W-X-Y) -Z*U*(X+Y+V+W-Z-U) -V*W*(X+Y+Z+U-V-W)
      RETURN
      END

      DOUBLE PRECISION FUNCTION DOT_S(XX,YY)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XX(4),YY(4)
      DOT_S= XX(4)*YY(4) -XX(1)*YY(1) -XX(2)*YY(2) -XX(3)*YY(3)
      RETURN
      END

      SUBROUTINE METRIC1
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/GMUNU1/DF(4,4),DC(4)
      SAVE
      DO 1 I=1,4
      DO 1 J=1,4
      DF(I,J)= 0.D0
    1 CONTINUE !P(i) =(Beam, P_t, P_t, Energy)
      DO 2 I=1,3
      DF(I,I)=-1.D0
    2 DC(I)  =-1.D0
      DF(4,4)= 1.D0
      DC(4)  = 1.D0
      RETURN
      END
