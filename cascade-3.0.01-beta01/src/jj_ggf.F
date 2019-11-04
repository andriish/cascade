      SUBROUTINE QEDSUM1(G1,a1, G2,a2, P1,P2)
C********************************************************************C
C    Contributions from the diagrams 11,12,13,14,21,22,23,31,32,33   C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G1(4),a1(4),G2(4),a2(4),P1(4),P2(4)
      COMMON/DOTPR1/W1,W2,G1G2,G1P1,G1P2,G2P1,G2P2,P1P2
     .       A1A2,A1G1,A1G2,A1P1,A1P2,A2G1,A2G2,A2P1,A2P2
      COMMON/QUARK1/XC,XC2,XB,XB2,XM,XM2,CHC,CHB
      COMMON/OTHER1/XM4,R,Q,R2,Q2,RQ,St,H1,H2
      COMMON/BEAMS1/IB1,IB2
      COMMON/AMPJJ1/AMP(4,4)
      COMMON/GMUNU1/DF(4,4),DC(4)
      SAVE
C
      Z11 = 8*(H1-G1P1*R)*(H2-G2P1*Q)*(H1+St*R2-G1P1*R-G1P2*R)
      Z12 = 8*(H1-G1P2*R)*(H2-G2P1*Q)*(H1+St*R2-G1P1*R-G1P2*R)
      Z13 = 8*(H1-G1P1*R)*(H2-G2P2*Q)*(H1+St*R2-G1P1*R-G1P2*R)
      Z14 = 8*(H1-G1P2*R)*(H2-G2P2*Q)*(H1+St*R2-G1P1*R-G1P2*R)
      Z21 = 8*(H1-G1P1*R) *(H1+H2+G1G2-G1P1*R-G2P1*R)* (St*Q2)
      Z22 = 8*(H1-G1P1*R)                *(H2-G2P2*R)* (St*Q2)
      Z23 = 8*(H1+H2+G1G2-G1P2*R-G2P2*R) *(H2-G2P2*R)* (St*Q2)
      Z31 = 8*(H2-G2P1*Q) *(H1+H2+G1G2-G1P1*Q-G2P1*Q)* (St*R2)
      Z32 = 8*(H2-G2P1*Q)                *(H1-G1P2*Q)* (St*R2)
      Z33 = 8*(H1+H2+G1G2-G1P2*Q-G2P2*Q) *(H1-G1P2*Q)* (St*R2)
C
      DO 88 n1=1,4
      DO 88 n2=1,4
      SUM11 =
     &  - 2*G1(n1)*G2(n2)*XM2*A1A2 - 2*G1(n1)*a1(n2)*XM2*A2G2 + 4*
     & G1(n1)*a1(n2)*XM2*Q*A2P1 + 2*G1(n1)*a2(n2)*XM2*A1G2 - 2*G1(n2)*
     & G2(n1)*XM2*A1A2 + 2*G1(n2)*a1(n1)*XM2*A2G2 - 4*G1(n2)*a1(n1)*XM2
     & *Q*A2P1 + 2*G1(n2)*a2(n1)*XM2*A1G2 + 2*G2(n1)*a1(n2)*XM2*A2G1 - 
     & 2*G2(n1)*a2(n2)*XM2*A1G1 + 4*G2(n1)*a2(n2)*XM2*R*A1P1 + 2*G2(n2)
     & *a1(n1)*XM2*A2G1 + 2*G2(n2)*a2(n1)*XM2*A1G1 - 4*G2(n2)*a2(n1)*
     & XM2*R*A1P1 + 2*a1(n1)*a2(n2)*RQ*XM4 - 2*a1(n1)*a2(n2)*R2*XM4 - 2
     & *a1(n1)*a2(n2)*XM2*G1G2 + 2*a1(n2)*a2(n1)*RQ*XM4 - 2*a1(n2)*
     & a2(n1)*R2*XM4 - 2*a1(n2)*a2(n1)*XM2*G1G2 - 2*DF(n1,n2)*RQ*XM4*
     & A1A2 + 2*DF(n1,n2)*R2*XM4*A1A2 - 2*DF(n1,n2)*XM2*A1G2*A2G1 + 2*
     & DF(n1,n2)*XM2*A1G1*A2G2 + 2*DF(n1,n2)*XM2*G1G2*A1A2 + 8*DF(n1,n2
     & )*XM2*RQ*A1P1*A2P1 - 4*DF(n1,n2)*XM2*Q*A1G1*A2P1 - 4*DF(n1,n2)*
     & XM2*R*A1P1*A2G2
C
      SUM12 =
     &  - 4*P1(n2)*P2(n1)*A1G2*A2G1 + 4*P1(n2)*P2(n1)*G1G2*A1A2 - 4*
     & P1(n2)*P2(n1)*RQ*P1P2*A1A2 + 4*P1(n2)*P2(n1)*XM2*R2*A1A2 + 4*
     & P1(n2)*G1(n1)*A1G2*A2P2 - 4*P1(n2)*G1(n1)*G2P2*A1A2 + 4*P1(n2)*
     & a1(n1)*G2P2*A2G1 - 4*P1(n2)*a1(n1)*G1G2*A2P2 + 4*P1(n2)*a1(n1)*
     & RQ*P1P2*A2P2 - 4*P1(n2)*a1(n1)*XM2*R2*A2P2 + 4*P2(n1)*G2(n2)*
     & A1P1*A2G1 - 4*P2(n1)*G2(n2)*G1P1*A1A2 + 4*P2(n1)*a2(n2)*G1P1*
     & A1G2 - 4*P2(n1)*a2(n2)*G1G2*A1P1 + 4*P2(n1)*a2(n2)*RQ*P1P2*A1P1
     &  - 4*P2(n1)*a2(n2)*XM2*R2*A1P1 - 4*G1(n1)*G2(n2)*A1P1*A2P2 + 4*
     & G1(n1)*G2(n2)*P1P2*A1A2 + 2*G1(n1)*G2(n2)*XM2*A1A2 - 2*G1(n1)*
     & a1(n2)*XM2*A2G2 + 4*G1(n1)*a1(n2)*XM2*Q*A2P1 - 4*G1(n1)*a2(n2)*
     & P1P2*A1G2 + 4*G1(n1)*a2(n2)*G2P2*A1P1 - 2*G1(n1)*a2(n2)*XM2*A1G2
     &  + 4*G1(n2)*G2(n1)*A1P2*A2P1 - 2*G1(n2)*G2(n1)*XM2*A1A2 + 2*
     & G1(n2)*a1(n1)*XM2*A2G2 - 4*G1(n2)*a1(n1)*XM2*Q*A2P1 - 4*G1(n2)*
     & a2(n1)*G2P1*A1P2 + 2*G1(n2)*a2(n1)*XM2*A1G2 - 4*G2(n1)*a1(n2)*
     & G1P2*A2P1
      SUM12 = SUM12 + 2*G2(n1)*a1(n2)*XM2*A2G1 + 2*G2(n1)*a2(n2)*XM2*
     & A1G1 - 4*G2(n1)*a2(n2)*XM2*R*A1P2 - 4*G2(n2)*a1(n1)*P1P2*A2G1 + 
     & 4*G2(n2)*a1(n1)*G1P1*A2P2 - 2*G2(n2)*a1(n1)*XM2*A2G1 - 2*G2(n2)*
     & a2(n1)*XM2*A1G1 + 4*G2(n2)*a2(n1)*XM2*R*A1P2 - 4*a1(n1)*a2(n2)*
     & G1P1*G2P2 + 4*a1(n1)*a2(n2)*G1G2*P1P2 - 4*a1(n1)*a2(n2)*RQ*
     & P1P2**2 + 2*a1(n1)*a2(n2)*R2*XM4 + 2*a1(n1)*a2(n2)*XM2*G1G2 - 2*
     & a1(n1)*a2(n2)*XM2*RQ*P1P2 + 4*a1(n1)*a2(n2)*XM2*R2*P1P2 + 4*
     & a1(n2)*a2(n1)*G1P2*G2P1 - 2*a1(n2)*a2(n1)*R2*XM4 - 2*a1(n2)*
     & a2(n1)*XM2*G1G2 + 2*a1(n2)*a2(n1)*XM2*RQ*P1P2 + 4*DF(n1,n2)*P1P2
     & *A1G2*A2G1 - 4*DF(n1,n2)*G2P2*A1P1*A2G1 - 4*DF(n1,n2)*G1P1*A1G2*
     & A2P2 + 4*DF(n1,n2)*G1P1*G2P2*A1A2 + 4*DF(n1,n2)*G1G2*A1P1*A2P2
     &  - 4*DF(n1,n2)*G1G2*P1P2*A1A2 - 4*DF(n1,n2)*RQ*P1P2*A1P1*A2P2 + 
     & 4*DF(n1,n2)*RQ*P1P2**2*A1A2 - 2*DF(n1,n2)*R2*XM4*A1A2 + 2*DF(n1,
     & n2)*XM2*A1G2*A2G1 - 2*DF(n1,n2)*XM2*A1G1*A2G2 - 2*DF(n1,n2)*XM2*
     & G1G2*A1A2 - 8*DF(n1,n2)*XM2*RQ*A1P2*A2P1 + 2*DF(n1,n2)*XM2*RQ*
     & P1P2*A1A2
      SUM12 = SUM12 + 4*DF(n1,n2)*XM2*R2*A1P1*A2P2 - 4*DF(n1,n2)*XM2*R2
     & *P1P2*A1A2 + 4*DF(n1,n2)*XM2*Q*A1G1*A2P1 + 4*DF(n1,n2)*XM2*R*
     & A1P2*A2G2
C
      SUM13 =
     &  - 4*P1(n2)*P2(n1)*A1G2*A2G1 + 4*P1(n2)*P2(n1)*G1G2*A1A2 - 4*
     & P1(n2)*P2(n1)*RQ*P1P2*A1A2 + 4*P1(n2)*P2(n1)*XM2*R2*A1A2 + 4*
     & P1(n2)*G2(n1)*A1P2*A2G1 - 4*P1(n2)*G2(n1)*G1P2*A1A2 + 4*P1(n2)*
     & a2(n1)*G1P2*A1G2 - 4*P1(n2)*a2(n1)*G1G2*A1P2 + 4*P1(n2)*a2(n1)*
     & RQ*P1P2*A1P2 - 4*P1(n2)*a2(n1)*XM2*R2*A1P2 + 4*P2(n1)*G1(n2)*
     & A1G2*A2P1 - 4*P2(n1)*G1(n2)*G2P1*A1A2 + 4*P2(n1)*a1(n2)*G2P1*
     & A2G1 - 4*P2(n1)*a1(n2)*G1G2*A2P1 + 4*P2(n1)*a1(n2)*RQ*P1P2*A2P1
     &  - 4*P2(n1)*a1(n2)*XM2*R2*A2P1 + 4*G1(n1)*G2(n2)*A1P1*A2P2 - 2*
     & G1(n1)*G2(n2)*XM2*A1A2 + 2*G1(n1)*a1(n2)*XM2*A2G2 - 4*G1(n1)*
     & a1(n2)*XM2*Q*A2P2 - 4*G1(n1)*a2(n2)*G2P2*A1P1 + 2*G1(n1)*a2(n2)*
     & XM2*A1G2 - 4*G1(n2)*G2(n1)*A1P2*A2P1 + 4*G1(n2)*G2(n1)*P1P2*A1A2
     &  + 2*G1(n2)*G2(n1)*XM2*A1A2 - 2*G1(n2)*a1(n1)*XM2*A2G2 + 4*
     & G1(n2)*a1(n1)*XM2*Q*A2P2 - 4*G1(n2)*a2(n1)*P1P2*A1G2 + 4*G1(n2)*
     & a2(n1)*G2P1*A1P2 - 2*G1(n2)*a2(n1)*XM2*A1G2 - 4*G2(n1)*a1(n2)*
     & P1P2*A2G1
      SUM13 = SUM13 + 4*G2(n1)*a1(n2)*G1P2*A2P1 - 2*G2(n1)*a1(n2)*XM2*
     & A2G1 - 2*G2(n1)*a2(n2)*XM2*A1G1 + 4*G2(n1)*a2(n2)*XM2*R*A1P1 - 4
     & *G2(n2)*a1(n1)*G1P1*A2P2 + 2*G2(n2)*a1(n1)*XM2*A2G1 + 2*G2(n2)*
     & a2(n1)*XM2*A1G1 - 4*G2(n2)*a2(n1)*XM2*R*A1P1 + 4*a1(n1)*a2(n2)*
     & G1P1*G2P2 - 2*a1(n1)*a2(n2)*R2*XM4 - 2*a1(n1)*a2(n2)*XM2*G1G2 + 
     & 2*a1(n1)*a2(n2)*XM2*RQ*P1P2 - 4*a1(n2)*a2(n1)*G1P2*G2P1 + 4*
     & a1(n2)*a2(n1)*G1G2*P1P2 - 4*a1(n2)*a2(n1)*RQ*P1P2**2 + 2*a1(n2)*
     & a2(n1)*R2*XM4 + 2*a1(n2)*a2(n1)*XM2*G1G2 - 2*a1(n2)*a2(n1)*XM2*
     & RQ*P1P2 + 4*a1(n2)*a2(n1)*XM2*R2*P1P2 + 4*DF(n1,n2)*P1P2*A1G2*
     & A2G1 - 4*DF(n1,n2)*G2P1*A1P2*A2G1 - 4*DF(n1,n2)*G1P2*A1G2*A2P1
     &  + 4*DF(n1,n2)*G1P2*G2P1*A1A2 + 4*DF(n1,n2)*G1G2*A1P2*A2P1 - 4*
     & DF(n1,n2)*G1G2*P1P2*A1A2 - 4*DF(n1,n2)*RQ*P1P2*A1P2*A2P1 + 4*DF(
     & n1,n2)*RQ*P1P2**2*A1A2 - 2*DF(n1,n2)*R2*XM4*A1A2 + 2*DF(n1,n2)*
     & XM2*A1G2*A2G1 - 2*DF(n1,n2)*XM2*A1G1*A2G2 - 2*DF(n1,n2)*XM2*G1G2
     & *A1A2 - 8*DF(n1,n2)*XM2*RQ*A1P1*A2P2 + 2*DF(n1,n2)*XM2*RQ*P1P2*
     & A1A2
      SUM13 = SUM13 + 4*DF(n1,n2)*XM2*R2*A1P2*A2P1 - 4*DF(n1,n2)*XM2*R2
     & *P1P2*A1A2 + 4*DF(n1,n2)*XM2*Q*A1G1*A2P2 + 4*DF(n1,n2)*XM2*R*
     & A1P1*A2G2
C
      SUM14 =
     &  - 2*G1(n1)*G2(n2)*XM2*A1A2 + 2*G1(n1)*a1(n2)*XM2*A2G2 - 4*
     & G1(n1)*a1(n2)*XM2*Q*A2P2 + 2*G1(n1)*a2(n2)*XM2*A1G2 - 2*G1(n2)*
     & G2(n1)*XM2*A1A2 - 2*G1(n2)*a1(n1)*XM2*A2G2 + 4*G1(n2)*a1(n1)*XM2
     & *Q*A2P2 + 2*G1(n2)*a2(n1)*XM2*A1G2 + 2*G2(n1)*a1(n2)*XM2*A2G1 + 
     & 2*G2(n1)*a2(n2)*XM2*A1G1 - 4*G2(n1)*a2(n2)*XM2*R*A1P2 + 2*G2(n2)
     & *a1(n1)*XM2*A2G1 - 2*G2(n2)*a2(n1)*XM2*A1G1 + 4*G2(n2)*a2(n1)*
     & XM2*R*A1P2 + 2*a1(n1)*a2(n2)*RQ*XM4 - 2*a1(n1)*a2(n2)*R2*XM4 - 2
     & *a1(n1)*a2(n2)*XM2*G1G2 + 2*a1(n2)*a2(n1)*RQ*XM4 - 2*a1(n2)*
     & a2(n1)*R2*XM4 - 2*a1(n2)*a2(n1)*XM2*G1G2 - 2*DF(n1,n2)*RQ*XM4*
     & A1A2 + 2*DF(n1,n2)*R2*XM4*A1A2 - 2*DF(n1,n2)*XM2*A1G2*A2G1 + 2*
     & DF(n1,n2)*XM2*A1G1*A2G2 + 2*DF(n1,n2)*XM2*G1G2*A1A2 + 8*DF(n1,n2
     & )*XM2*RQ*A1P2*A2P2 - 4*DF(n1,n2)*XM2*Q*A1G1*A2P2 - 4*DF(n1,n2)*
     & XM2*R*A1P2*A2G2
C
      SUM21 =
     & 2*G1(n1)*G2(n2)*XM2*A1A2 + 2*G1(n1)*a1(n2)*XM2*A2G2 + 4*G1(n1)*
     & a1(n2)*XM2*A2G1 - 4*G1(n1)*a1(n2)*XM2*R*A2P1 - 2*G1(n1)*a2(n2)*
     & XM2*A1G2 - 2*G1(n2)*G2(n1)*XM2*A1A2 - 2*G1(n2)*a1(n1)*XM2*A2G2
     &  - 4*G1(n2)*a1(n1)*XM2*A2G1 + 4*G1(n2)*a1(n1)*XM2*R*A2P1 + 2*
     & G1(n2)*a2(n1)*XM2*A1G2 + 2*G2(n1)*a1(n2)*XM2*A2G1 + 2*G2(n1)*
     & a2(n2)*XM2*A1G1 - 4*G2(n1)*a2(n2)*XM2*R*A1P1 - 2*G2(n2)*a1(n1)*
     & XM2*A2G1 - 2*G2(n2)*a2(n1)*XM2*A1G1 + 4*G2(n2)*a2(n1)*XM2*R*A1P1
     &  + 2*a1(n1)*a2(n2)*XM2*G1G2 + 2*a1(n1)*a2(n2)*XM2*W1 - 4*a1(n1)*
     & a2(n2)*XM2*R*G1P1 - 2*a1(n2)*a2(n1)*XM2*G1G2 - 2*a1(n2)*a2(n1)*
     & XM2*W1 + 4*a1(n2)*a2(n1)*XM2*R*G1P1 - 2*DF(n1,n2)*XM2*A1G2*A2G1
     &  - 2*DF(n1,n2)*XM2*A1G1*A2G2 - 4*DF(n1,n2)*XM2*A1G1*A2G1 + 2*DF(
     & n1,n2)*XM2*G1G2*A1A2 + 2*DF(n1,n2)*XM2*W1*A1A2 - 8*DF(n1,n2)*XM2
     & *R2*A1P1*A2P1 + 4*DF(n1,n2)*XM2*R*A1P1*A2G2 + 8*DF(n1,n2)*XM2*R*
     & A1P1*A2G1 + 4*DF(n1,n2)*XM2*R*A1G1*A2P1 - 4*DF(n1,n2)*XM2*R*G1P1
     & *A1A2
C
      SUM22 =
     &  - 4*P1(n2)*P2(n1)*A1G2*A2G1 + 4*P1(n2)*P2(n1)*G1G2*A1A2 - 4*
     & P1(n2)*P2(n1)*R2*P1P2*A1A2 + 4*P1(n2)*P2(n1)*XM2*R2*A1A2 + 4*
     & P1(n2)*G2(n1)*A1P2*A2G1 - 4*P1(n2)*G2(n1)*G1P2*A1A2 + 4*P1(n2)*
     & a2(n1)*G1P2*A1G2 - 4*P1(n2)*a2(n1)*G1G2*A1P2 + 4*P1(n2)*a2(n1)*
     & R2*P1P2*A1P2 - 4*P1(n2)*a2(n1)*XM2*R2*A1P2 + 4*P2(n1)*G1(n2)*
     & A1G2*A2P1 - 4*P2(n1)*G1(n2)*G2P1*A1A2 + 4*P2(n1)*a1(n2)*G2P1*
     & A2G1 - 4*P2(n1)*a1(n2)*G1G2*A2P1 + 4*P2(n1)*a1(n2)*R2*P1P2*A2P1
     &  - 4*P2(n1)*a1(n2)*XM2*R2*A2P1 + 4*G1(n1)*G2(n2)*A1P1*A2P2 + 2*
     & G1(n1)*G2(n2)*XM2*A1A2 - 2*G1(n1)*a1(n2)*XM2*A2G2 + 4*G1(n1)*
     & a1(n2)*XM2*R*A2P2 - 4*G1(n1)*a2(n2)*G2P2*A1P1 - 2*G1(n1)*a2(n2)*
     & XM2*A1G2 - 4*G1(n2)*G2(n1)*A1P2*A2P1 + 4*G1(n2)*G2(n1)*P1P2*A1A2
     &  + 2*G1(n2)*G2(n1)*XM2*A1A2 + 2*G1(n2)*a1(n1)*XM2*A2G2 - 4*
     & G1(n2)*a1(n1)*XM2*R*A2P2 - 4*G1(n2)*a2(n1)*P1P2*A1G2 + 4*G1(n2)*
     & a2(n1)*G2P1*A1P2 - 2*G1(n2)*a2(n1)*XM2*A1G2 - 4*G2(n1)*a1(n2)*
     & P1P2*A2G1 - 4*DF(n1,n2)*XM2*R*A1P1*A2G2
      SUM22 = SUM22 + 4*G2(n1)*a1(n2)*G1P2*A2P1 - 2*G2(n1)*a1(n2)*XM2*
     & A2G1 + 2*G2(n1)*a2(n2)*XM2*A1G1 - 4*G2(n1)*a2(n2)*XM2*R*A1P1 - 4
     & *G2(n2)*a1(n1)*G1P1*A2P2 - 2*G2(n2)*a1(n1)*XM2*A2G1 - 2*G2(n2)*
     & a2(n1)*XM2*A1G1 + 4*G2(n2)*a2(n1)*XM2*R*A1P1 + 4*a1(n1)*a2(n2)*
     & G1P1*G2P2 + 2*a1(n1)*a2(n2)*R2*XM4 + 2*a1(n1)*a2(n2)*XM2*G1G2 - 
     & 2*a1(n1)*a2(n2)*XM2*R2*P1P2 - 4*a1(n2)*a2(n1)*G1P2*G2P1 + 4*
     & a1(n2)*a2(n1)*G1G2*P1P2 - 4*a1(n2)*a2(n1)*R2*P1P2**2 + 2*a1(n2)*
     & a2(n1)*R2*XM4 + 2*a1(n2)*a2(n1)*XM2*G1G2 + 2*a1(n2)*a2(n1)*XM2*
     & R2*P1P2 + 4*DF(n1,n2)*P1P2*A1G2*A2G1 - 4*DF(n1,n2)*G2P1*A1P2*
     & A2G1 - 4*DF(n1,n2)*G1P2*A1G2*A2P1 + 4*DF(n1,n2)*G1P2*G2P1*A1A2
     &  + 4*DF(n1,n2)*G1G2*A1P2*A2P1 - 4*DF(n1,n2)*G1G2*P1P2*A1A2 - 4*
     & DF(n1,n2)*R2*P1P2*A1P2*A2P1 + 4*DF(n1,n2)*R2*P1P2**2*A1A2 - 2*
     & DF(n1,n2)*R2*XM4*A1A2 + 2*DF(n1,n2)*XM2*A1G2*A2G1 + 2*DF(n1,n2)*
     & XM2*A1G1*A2G2 - 2*DF(n1,n2)*XM2*G1G2*A1A2 + 4*DF(n1,n2)*XM2*R2*
     & A1P2*A2P1 + 8*DF(n1,n2)*XM2*R2*A1P1*A2P2 - 2*DF(n1,n2)*XM2*R2*
     & P1P2*A1A2 - 4*DF(n1,n2)*XM2*R*A1G1*A2P2
C
      SUM23 =
     & 2*G1(n1)*G2(n2)*XM2*A1A2 - 2*G1(n1)*a1(n2)*XM2*A2G2 + 4*G1(n1)*
     & a1(n2)*XM2*R*A2P2 - 2*G1(n1)*a2(n2)*XM2*A1G2 - 2*G1(n2)*G2(n1)*
     & XM2*A1A2 + 2*G1(n2)*a1(n1)*XM2*A2G2 - 4*G1(n2)*a1(n1)*XM2*R*A2P2
     &  + 2*G1(n2)*a2(n1)*XM2*A1G2 + 2*G2(n1)*a1(n2)*XM2*A2G1 - 4*
     & G2(n1)*a2(n2)*XM2*A1G2 - 2*G2(n1)*a2(n2)*XM2*A1G1 + 4*G2(n1)*
     & a2(n2)*XM2*R*A1P2 - 2*G2(n2)*a1(n1)*XM2*A2G1 + 4*G2(n2)*a2(n1)*
     & XM2*A1G2 + 2*G2(n2)*a2(n1)*XM2*A1G1 - 4*G2(n2)*a2(n1)*XM2*R*A1P2
     &  + 2*a1(n1)*a2(n2)*XM2*G1G2 + 2*a1(n1)*a2(n2)*XM2*W2 - 4*a1(n1)*
     & a2(n2)*XM2*R*G2P2 - 2*a1(n2)*a2(n1)*XM2*G1G2 - 2*a1(n2)*a2(n1)*
     & XM2*W2 + 4*a1(n2)*a2(n1)*XM2*R*G2P2 - 4*DF(n1,n2)*XM2*A1G2*A2G2
     &  - 2*DF(n1,n2)*XM2*A1G2*A2G1 - 2*DF(n1,n2)*XM2*A1G1*A2G2 + 2*DF(
     & n1,n2)*XM2*G1G2*A1A2 + 2*DF(n1,n2)*XM2*W2*A1A2 - 8*DF(n1,n2)*XM2
     & *R2*A1P2*A2P2 + 4*DF(n1,n2)*XM2*R*A1P2*A2G2 + 8*DF(n1,n2)*XM2*R*
     & A1G2*A2P2 + 4*DF(n1,n2)*XM2*R*A1G1*A2P2 - 4*DF(n1,n2)*XM2*R*G2P2
     & *A1A2
C
      SUM31 =
     &  - 2*G1(n1)*G2(n2)*XM2*A1A2 + 2*G1(n1)*a1(n2)*XM2*A2G2 - 4*
     & G1(n1)*a1(n2)*XM2*Q*A2P1 + 2*G1(n1)*a2(n2)*XM2*A1G2 + 2*G1(n2)*
     & G2(n1)*XM2*A1A2 - 2*G1(n2)*a1(n1)*XM2*A2G2 + 4*G1(n2)*a1(n1)*XM2
     & *Q*A2P1 - 2*G1(n2)*a2(n1)*XM2*A1G2 - 2*G2(n1)*a1(n2)*XM2*A2G1 + 
     & 4*G2(n1)*a2(n2)*XM2*A1G2 + 2*G2(n1)*a2(n2)*XM2*A1G1 - 4*G2(n1)*
     & a2(n2)*XM2*Q*A1P1 + 2*G2(n2)*a1(n1)*XM2*A2G1 - 4*G2(n2)*a2(n1)*
     & XM2*A1G2 - 2*G2(n2)*a2(n1)*XM2*A1G1 + 4*G2(n2)*a2(n1)*XM2*Q*A1P1
     &  - 2*a1(n1)*a2(n2)*XM2*G1G2 - 2*a1(n1)*a2(n2)*XM2*W2 + 4*a1(n1)*
     & a2(n2)*XM2*Q*G2P1 + 2*a1(n2)*a2(n1)*XM2*G1G2 + 2*a1(n2)*a2(n1)*
     & XM2*W2 - 4*a1(n2)*a2(n1)*XM2*Q*G2P1 - 4*DF(n1,n2)*XM2*A1G2*A2G2
     &  - 2*DF(n1,n2)*XM2*A1G2*A2G1 - 2*DF(n1,n2)*XM2*A1G1*A2G2 + 2*DF(
     & n1,n2)*XM2*G1G2*A1A2 + 2*DF(n1,n2)*XM2*W2*A1A2 - 8*DF(n1,n2)*XM2
     & *Q2*A1P1*A2P1 + 4*DF(n1,n2)*XM2*Q*A1P1*A2G2 + 8*DF(n1,n2)*XM2*Q*
     & A1G2*A2P1 + 4*DF(n1,n2)*XM2*Q*A1G1*A2P1 - 4*DF(n1,n2)*XM2*Q*G2P1
     & *A1A2
C
      SUM32 =
     &  - 4*P1(n2)*P2(n1)*A1G2*A2G1 + 4*P1(n2)*P2(n1)*G1G2*A1A2 - 4*
     & P1(n2)*P2(n1)*Q2*P1P2*A1A2 + 4*P1(n2)*P2(n1)*XM2*Q2*A1A2 + 4*
     & P1(n2)*G1(n1)*A1G2*A2P2 - 4*P1(n2)*G1(n1)*G2P2*A1A2 + 4*P1(n2)*
     & a1(n1)*G2P2*A2G1 - 4*P1(n2)*a1(n1)*G1G2*A2P2 + 4*P1(n2)*a1(n1)*
     & Q2*P1P2*A2P2 - 4*P1(n2)*a1(n1)*XM2*Q2*A2P2 + 4*P2(n1)*G2(n2)*
     & A1P1*A2G1 - 4*P2(n1)*G2(n2)*G1P1*A1A2 + 4*P2(n1)*a2(n2)*G1P1*
     & A1G2 - 4*P2(n1)*a2(n2)*G1G2*A1P1 + 4*P2(n1)*a2(n2)*Q2*P1P2*A1P1
     &  - 4*P2(n1)*a2(n2)*XM2*Q2*A1P1 - 4*G1(n1)*G2(n2)*A1P1*A2P2 + 4*
     & G1(n1)*G2(n2)*P1P2*A1A2 + 2*G1(n1)*G2(n2)*XM2*A1A2 + 2*G1(n1)*
     & a1(n2)*XM2*A2G2 - 4*G1(n1)*a1(n2)*XM2*Q*A2P1 - 4*G1(n1)*a2(n2)*
     & P1P2*A1G2 + 4*G1(n1)*a2(n2)*G2P2*A1P1 - 2*G1(n1)*a2(n2)*XM2*A1G2
     &  + 4*G1(n2)*G2(n1)*A1P2*A2P1 + 2*G1(n2)*G2(n1)*XM2*A1A2 - 2*
     & G1(n2)*a1(n1)*XM2*A2G2 + 4*G1(n2)*a1(n1)*XM2*Q*A2P1 - 4*G1(n2)*
     & a2(n1)*G2P1*A1P2 - 2*G1(n2)*a2(n1)*XM2*A1G2 - 4*G2(n1)*a1(n2)*
     & G1P2*A2P1  - 4*DF(n1,n2)*XM2*Q*A1P2*A2G2
      SUM32 = SUM32 - 2*G2(n1)*a1(n2)*XM2*A2G1 - 2*G2(n1)*a2(n2)*XM2*
     & A1G1 + 4*G2(n1)*a2(n2)*XM2*Q*A1P2 - 4*G2(n2)*a1(n1)*P1P2*A2G1 + 
     & 4*G2(n2)*a1(n1)*G1P1*A2P2 - 2*G2(n2)*a1(n1)*XM2*A2G1 + 2*G2(n2)*
     & a2(n1)*XM2*A1G1 - 4*G2(n2)*a2(n1)*XM2*Q*A1P2 - 4*a1(n1)*a2(n2)*
     & G1P1*G2P2 + 4*a1(n1)*a2(n2)*G1G2*P1P2 - 4*a1(n1)*a2(n2)*Q2*
     & P1P2**2 + 2*a1(n1)*a2(n2)*Q2*XM4 + 2*a1(n1)*a2(n2)*XM2*G1G2 + 2*
     & a1(n1)*a2(n2)*XM2*Q2*P1P2 + 4*a1(n2)*a2(n1)*G1P2*G2P1 + 2*a1(n2)
     & *a2(n1)*Q2*XM4 + 2*a1(n2)*a2(n1)*XM2*G1G2 - 2*a1(n2)*a2(n1)*XM2*
     & Q2*P1P2 + 4*DF(n1,n2)*P1P2*A1G2*A2G1 - 4*DF(n1,n2)*G2P2*A1P1*
     & A2G1 - 4*DF(n1,n2)*G1P1*A1G2*A2P2 + 4*DF(n1,n2)*G1P1*G2P2*A1A2
     &  + 4*DF(n1,n2)*G1G2*A1P1*A2P2 - 4*DF(n1,n2)*G1G2*P1P2*A1A2 - 4*
     & DF(n1,n2)*Q2*P1P2*A1P1*A2P2 + 4*DF(n1,n2)*Q2*P1P2**2*A1A2 - 2*
     & DF(n1,n2)*Q2*XM4*A1A2 + 2*DF(n1,n2)*XM2*A1G2*A2G1 + 2*DF(n1,n2)*
     & XM2*A1G1*A2G2 - 2*DF(n1,n2)*XM2*G1G2*A1A2 + 8*DF(n1,n2)*XM2*Q2*
     & A1P2*A2P1 + 4*DF(n1,n2)*XM2*Q2*A1P1*A2P2 - 2*DF(n1,n2)*XM2*Q2*
     & P1P2*A1A2 - 4*DF(n1,n2)*XM2*Q*A1G1*A2P1
C
      SUM33 =
     &  - 2*G1(n1)*G2(n2)*XM2*A1A2 - 2*G1(n1)*a1(n2)*XM2*A2G2 - 4*
     & G1(n1)*a1(n2)*XM2*A2G1 + 4*G1(n1)*a1(n2)*XM2*Q*A2P2 + 2*G1(n1)*
     & a2(n2)*XM2*A1G2 + 2*G1(n2)*G2(n1)*XM2*A1A2 + 2*G1(n2)*a1(n1)*XM2
     & *A2G2 + 4*G1(n2)*a1(n1)*XM2*A2G1 - 4*G1(n2)*a1(n1)*XM2*Q*A2P2 - 
     & 2*G1(n2)*a2(n1)*XM2*A1G2 - 2*G2(n1)*a1(n2)*XM2*A2G1 - 2*G2(n1)*
     & a2(n2)*XM2*A1G1 + 4*G2(n1)*a2(n2)*XM2*Q*A1P2 + 2*G2(n2)*a1(n1)*
     & XM2*A2G1 + 2*G2(n2)*a2(n1)*XM2*A1G1 - 4*G2(n2)*a2(n1)*XM2*Q*A1P2
     &  - 2*a1(n1)*a2(n2)*XM2*G1G2 - 2*a1(n1)*a2(n2)*XM2*W1 + 4*a1(n1)*
     & a2(n2)*XM2*Q*G1P2 + 2*a1(n2)*a2(n1)*XM2*G1G2 + 2*a1(n2)*a2(n1)*
     & XM2*W1 - 4*a1(n2)*a2(n1)*XM2*Q*G1P2 - 2*DF(n1,n2)*XM2*A1G2*A2G1
     &  - 2*DF(n1,n2)*XM2*A1G1*A2G2 - 4*DF(n1,n2)*XM2*A1G1*A2G1 + 2*DF(
     & n1,n2)*XM2*G1G2*A1A2 + 2*DF(n1,n2)*XM2*W1*A1A2 - 8*DF(n1,n2)*XM2
     & *Q2*A1P2*A2P2 + 4*DF(n1,n2)*XM2*Q*A1P2*A2G2 + 8*DF(n1,n2)*XM2*Q*
     & A1P2*A2G1 + 4*DF(n1,n2)*XM2*Q*A1G1*A2P2 - 4*DF(n1,n2)*XM2*Q*G1P2
     & *A1A2
C
      IF(IB1.EQ.8 .AND. IB2.EQ.8) THEN         !For Gluon*Gluon
      AMP(n1,n2) = AMP(n1,n2)
     . + SUM11/Z11 *(16./3.) + SUM12/Z12 *(-2./3.)
     . + SUM13/Z13 *(-2./3.) + SUM14/Z14 *(16./3.)
     . + SUM21/Z21 *(16./3.) + SUM31/Z31 *(16./3.)
     . + SUM22/Z22 *(-2./3.) + SUM32/Z32 *(-2./3.)
     . + SUM23/Z23 *(16./3.) + SUM33/Z33 *(16./3.)
      ENDIF
      IF(IB1.EQ.1 .AND. IB2.EQ.1) THEN         !For Gamma*Gamma
      AMP(n1,n2) = AMP(n1,n2)
     . + SUM11/Z11 *4*CHC*CHB + SUM12/Z12 *4*CHC*CHB
     . + SUM13/Z13 *4*CHC*CHB + SUM14/Z14 *4*CHC*CHB
     . + SUM21/Z21 *4*CHC*CHC + SUM31/Z31 *4*CHB*CHB
     . + SUM22/Z22 *4*CHC*CHC + SUM32/Z32 *4*CHB*CHB
     . + SUM23/Z23 *4*CHC*CHC + SUM33/Z33 *4*CHB*CHB
      ENDIF
   88 CONTINUE
C
      RETURN
      END

      SUBROUTINE QCDSUM1(G1,a1, G2,a2, P1,P2)
C********************************************************************C
C    Contributions from the diagrams 41,42,45,46,62,64               C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G1(4),a1(4),G2(4),a2(4),P1(4),P2(4)
      COMMON/DOTPR1/W1,W2,G1G2,G1P1,G1P2,G2P1,G2P2,P1P2
     .       A1A2,A1G1,A1G2,A1P1,A1P2,A2G1,A2G2,A2P1,A2P2
      COMMON/QUARK1/XC,XC2,XB,XB2,XM,XM2,CHC,CHB
      COMMON/OTHER1/XM4,R,Q,R2,Q2,RQ,St,H1,H2
      COMMON/BEAMS1/IB1,IB2
      COMMON/AMPJJ1/AMP(4,4)
      COMMON/GMUNU1/DF(4,4),DC(4)
      SAVE
      IF(IB1.NE.8 .OR. IB2.NE.8) RETURN    !For Gluon*Gluon only
C
      Z41 =-8*(H1-G1P1*R) *(St*Q2)*(H2+St*Q2-G2P2*Q-G2P1*Q)
      Z42 =-8*(H1-G1P2*R) *(St*Q2)*(H2+St*Q2-G2P2*Q-G2P1*Q)
      Z45 =-8*(H1-G1P1*Q) *(St*R2)*(H2+St*R2-G2P2*R-G2P1*R)
      Z46 =-8*(H1-G1P2*Q) *(St*R2)*(H2+St*R2-G2P2*R-G2P1*R)
      Z62 =-8*(St*R2)*(St*Q2) *(H1+St*R2-G1P1*R-G1P2*R)
      Z64 =-4*(St*R2)*(St*Q2)
C
      DO 88 n1=1,4
      DO 88 n2=1,4
      SUM41 =
     &  - 3*P1(n2)*P2(n1)*A1G2*A2G1 - P1(n2)*P2(n1)*A1G1*A2G2 + 3*
     & P1(n2)*P2(n1)*G1G2*A1A2 - 2*P1(n2)*P2(n1)*Q*A1P1*A2G1 + 2*P1(n2)
     & *P2(n1)*Q*G1P1*A1A2 + 2*P1(n2)*P2(n1)*R*A1P1*A2G2 + P1(n2)*
     & G1(n1)*A1P2*A2G2 - P1(n2)*G1(n1)*A1G2*A2P2 - P1(n2)*G1(n1)*G2P2*
     & A1A2 + 2*P1(n2)*G1(n1)*Q*A1P1*A2P2 + 3*P1(n2)*G2(n1)*A1P2*A2G1
     &  + P1(n2)*G2(n1)*A1G1*A2P2 - 3*P1(n2)*G2(n1)*G1P2*A1A2 - 2*
     & P1(n2)*G2(n1)*R*A1P1*A2P2 + P1(n2)*a1(n1)*G2P2*A2G1 - P1(n2)*
     & a1(n1)*G1P2*A2G2 + P1(n2)*a1(n1)*G1G2*A2P2 - 2*P1(n2)*a1(n1)*Q*
     & G1P1*A2P2 + P1(n2)*a2(n1)*G2P2*A1G1 + 3*P1(n2)*a2(n1)*G1P2*A1G2
     &  - 3*P1(n2)*a2(n1)*G1G2*A1P2 + 2*P1(n2)*a2(n1)*Q*G1P2*A1P1 - 2*
     & P1(n2)*a2(n1)*Q*G1P1*A1P2 - 2*P1(n2)*a2(n1)*R*G2P2*A1P1 - P2(n1)
     & *G1(n2)*A1P1*A2G2 + 3*P2(n1)*G1(n2)*A1G2*A2P1 - 3*P2(n1)*G1(n2)*
     & G2P1*A1A2 + 2*P2(n1)*G1(n2)*Q*A1P1*A2P1 + P2(n1)*G2(n2)*A1P1*
     & A2G1 + P2(n1)*G2(n2)*A1G1*A2P1 - P2(n1)*G2(n2)*G1P1*A1A2 - 2*
     & P2(n1)*G2(n2)*R*A1P1*A2P1
      SUM41 = SUM41 + 3*P2(n1)*a1(n2)*G2P1*A2G1 + P2(n1)*a1(n2)*G1P1*
     & A2G2 - 3*P2(n1)*a1(n2)*G1G2*A2P1 - 2*P2(n1)*a1(n2)*Q*G1P1*A2P1
     &  + P2(n1)*a2(n2)*G2P1*A1G1 - P2(n1)*a2(n2)*G1P1*A1G2 + P2(n1)*
     & a2(n2)*G1G2*A1P1 - 2*P2(n1)*a2(n2)*R*G2P1*A1P1 - G1(n1)*G2(n2)*
     & A1P2*A2P1 + 3*G1(n1)*G2(n2)*A1P1*A2P2 + G1(n1)*G2(n2)*St*A1A2 + 
     & G1(n1)*a1(n2)*G2P2*A2P1 + G1(n1)*a1(n2)*G2P1*A2P2 - G1(n1)*
     & a1(n2)*St*A2G2 + 2*G1(n1)*a1(n2)*XM2*A2G2 - 4*G1(n1)*a1(n2)*XM2*
     & Q*A2P2 - 4*G1(n1)*a1(n2)*XM2*Q*A2P1 - 3*G1(n1)*a2(n2)*G2P2*A1P1
     &  - G1(n1)*a2(n2)*G2P1*A1P2 + G1(n1)*a2(n2)*St*A1G2 - 2*G1(n1)*
     & a2(n2)*Q*St*A1P1 - 3*G1(n2)*G2(n1)*A1P2*A2P1 + G1(n2)*G2(n1)*
     & A1P1*A2P2 + 3*G1(n2)*G2(n1)*St*A1A2 - G1(n2)*a1(n1)*G2P2*A2P1 - 
     & G1(n2)*a1(n1)*G2P1*A2P2 + G1(n2)*a1(n1)*St*A2G2 - 2*G1(n2)*
     & a1(n1)*XM2*A2G2 + 4*G1(n2)*a1(n1)*XM2*Q*A2P2 + 4*G1(n2)*a1(n1)*
     & XM2*Q*A2P1 + G1(n2)*a2(n1)*G2P2*A1P1 + 3*G1(n2)*a2(n1)*G2P1*A1P2
     &  - 3*G1(n2)*a2(n1)*St*A1G2 - 2*G1(n2)*a2(n1)*Q*St*A1P1 + 3*
     & G2(n1)*a1(n2)*G1P2*A2P1
      SUM41 = SUM41 - G2(n1)*a1(n2)*G1P1*A2P2 - 3*G2(n1)*a1(n2)*St*A2G1
     &  - G2(n1)*a2(n2)*G1P2*A1P1 + G2(n1)*a2(n2)*G1P1*A1P2 - G2(n1)*
     & a2(n2)*St*A1G1 + 2*G2(n1)*a2(n2)*R*St*A1P1 + G2(n2)*a1(n1)*G1P2*
     & A2P1 - 3*G2(n2)*a1(n1)*G1P1*A2P2 - G2(n2)*a1(n1)*St*A2G1 - 
     & G2(n2)*a2(n1)*G1P2*A1P1 + G2(n2)*a2(n1)*G1P1*A1P2 - G2(n2)*
     & a2(n1)*St*A1G1 + 2*G2(n2)*a2(n1)*R*St*A1P1 + a1(n1)*a2(n2)*G1P2*
     & G2P1 + 3*a1(n1)*a2(n2)*G1P1*G2P2 - a1(n1)*a2(n2)*St*G1G2 + 2*
     & a1(n1)*a2(n2)*Q*St*G1P1 - 3*a1(n2)*a2(n1)*G1P2*G2P1 - a1(n2)*
     & a2(n1)*G1P1*G2P2 + 3*a1(n2)*a2(n1)*St*G1G2 + 2*a1(n2)*a2(n1)*Q*
     & St*G1P1 - DF(n1,n2)*G2P2*A1P1*A2G1 - DF(n1,n2)*G2P2*A1G1*A2P1 - 
     & 3*DF(n1,n2)*G2P1*A1P2*A2G1 - DF(n1,n2)*G2P1*A1G1*A2P2 + DF(n1,n2
     & )*G1P2*A1P1*A2G2 - 3*DF(n1,n2)*G1P2*A1G2*A2P1 + 3*DF(n1,n2)*G1P2
     & *G2P1*A1A2 - DF(n1,n2)*G1P1*A1P2*A2G2 + DF(n1,n2)*G1P1*A1G2*A2P2
     &  + DF(n1,n2)*G1P1*G2P2*A1A2 + 3*DF(n1,n2)*G1G2*A1P2*A2P1 - DF(n1
     & ,n2)*G1G2*A1P1*A2P2 + 3*DF(n1,n2)*St*A1G2*A2G1 + DF(n1,n2)*St*
     & A1G1*A2G2
      SUM41 = SUM41 - 3*DF(n1,n2)*St*G1G2*A1A2 - 2*DF(n1,n2)*Q*G1P2*
     & A1P1*A2P1 + 2*DF(n1,n2)*Q*G1P1*A1P2*A2P1 + 2*DF(n1,n2)*Q*St*A1P1
     & *A2G1 - 2*DF(n1,n2)*Q*St*G1P1*A1A2 + 2*DF(n1,n2)*R*G2P2*A1P1*
     & A2P1 + 2*DF(n1,n2)*R*G2P1*A1P1*A2P2 - 2*DF(n1,n2)*R*St*A1P1*A2G2
     &  - 2*DF(n1,n2)*XM2*A1G1*A2G2 - 8*DF(n1,n2)*XM2*RQ*A1P1*A2P2 - 8*
     & DF(n1,n2)*XM2*RQ*A1P1*A2P1 + 4*DF(n1,n2)*XM2*Q*A1G1*A2P2 + 4*DF(
     & n1,n2)*XM2*Q*A1G1*A2P1 + 4*DF(n1,n2)*XM2*R*A1P1*A2G2
C
      SUM42 =
     & 3*P1(n2)*P2(n1)*A1G2*A2G1 + P1(n2)*P2(n1)*A1G1*A2G2 - 3*P1(n2)*
     & P2(n1)*G1G2*A1A2 + 2*P1(n2)*P2(n1)*Q*A1P2*A2G1 - 2*P1(n2)*P2(n1)
     & *Q*G1P2*A1A2 - 2*P1(n2)*P2(n1)*R*A1P2*A2G2 + P1(n2)*G1(n1)*A1P2*
     & A2G2 - 3*P1(n2)*G1(n1)*A1G2*A2P2 + 3*P1(n2)*G1(n1)*G2P2*A1A2 - 2
     & *P1(n2)*G1(n1)*Q*A1P2*A2P2 - P1(n2)*G2(n1)*A1P2*A2G1 - P1(n2)*
     & G2(n1)*A1G1*A2P2 + P1(n2)*G2(n1)*G1P2*A1A2 + 2*P1(n2)*G2(n1)*R*
     & A1P2*A2P2 - 3*P1(n2)*a1(n1)*G2P2*A2G1 - P1(n2)*a1(n1)*G1P2*A2G2
     &  + 3*P1(n2)*a1(n1)*G1G2*A2P2 + 2*P1(n2)*a1(n1)*Q*G1P2*A2P2 - 
     & P1(n2)*a2(n1)*G2P2*A1G1 + P1(n2)*a2(n1)*G1P2*A1G2 - P1(n2)*
     & a2(n1)*G1G2*A1P2 + 2*P1(n2)*a2(n1)*R*G2P2*A1P2 - P2(n1)*G1(n2)*
     & A1P1*A2G2 + P2(n1)*G1(n2)*A1G2*A2P1 + P2(n1)*G1(n2)*G2P1*A1A2 - 
     & 2*P2(n1)*G1(n2)*Q*A1P2*A2P1 - 3*P2(n1)*G2(n2)*A1P1*A2G1 - P2(n1)
     & *G2(n2)*A1G1*A2P1 + 3*P2(n1)*G2(n2)*G1P1*A1A2 + 2*P2(n1)*G2(n2)*
     & R*A1P2*A2P1 - P2(n1)*a1(n2)*G2P1*A2G1 + P2(n1)*a1(n2)*G1P1*A2G2
     &  - P2(n1)*a1(n2)*G1G2*A2P1
      SUM42 = SUM42 + 2*P2(n1)*a1(n2)*Q*G1P2*A2P1 - P2(n1)*a2(n2)*G2P1*
     & A1G1 - 3*P2(n1)*a2(n2)*G1P1*A1G2 + 3*P2(n1)*a2(n2)*G1G2*A1P1 + 2
     & *P2(n1)*a2(n2)*Q*G1P2*A1P1 - 2*P2(n1)*a2(n2)*Q*G1P1*A1P2 + 2*
     & P2(n1)*a2(n2)*R*G2P1*A1P2 - G1(n1)*G2(n2)*A1P2*A2P1 + 3*G1(n1)*
     & G2(n2)*A1P1*A2P2 - 3*G1(n1)*G2(n2)*St*A1A2 + G1(n1)*a1(n2)*G2P2*
     & A2P1 + G1(n1)*a1(n2)*G2P1*A2P2 - G1(n1)*a1(n2)*St*A2G2 + 2*
     & G1(n1)*a1(n2)*XM2*A2G2 - 4*G1(n1)*a1(n2)*XM2*Q*A2P2 - 4*G1(n1)*
     & a1(n2)*XM2*Q*A2P1 - 3*G1(n1)*a2(n2)*G2P2*A1P1 - G1(n1)*a2(n2)*
     & G2P1*A1P2 + 3*G1(n1)*a2(n2)*St*A1G2 + 2*G1(n1)*a2(n2)*Q*St*A1P2
     &  - 3*G1(n2)*G2(n1)*A1P2*A2P1 + G1(n2)*G2(n1)*A1P1*A2P2 - G1(n2)*
     & G2(n1)*St*A1A2 - G1(n2)*a1(n1)*G2P2*A2P1 - G1(n2)*a1(n1)*G2P1*
     & A2P2 + G1(n2)*a1(n1)*St*A2G2 - 2*G1(n2)*a1(n1)*XM2*A2G2 + 4*
     & G1(n2)*a1(n1)*XM2*Q*A2P2 + 4*G1(n2)*a1(n1)*XM2*Q*A2P1 + G1(n2)*
     & a2(n1)*G2P2*A1P1 + 3*G1(n2)*a2(n1)*G2P1*A1P2 - G1(n2)*a2(n1)*St*
     & A1G2 + 2*G1(n2)*a2(n1)*Q*St*A1P2 + 3*G2(n1)*a1(n2)*G1P2*A2P1 - 
     & G2(n1)*a1(n2)*G1P1*A2P2
      SUM42 = SUM42 + G2(n1)*a1(n2)*St*A2G1 - G2(n1)*a2(n2)*G1P2*A1P1
     &  + G2(n1)*a2(n2)*G1P1*A1P2 + G2(n1)*a2(n2)*St*A1G1 - 2*G2(n1)*
     & a2(n2)*R*St*A1P2 + G2(n2)*a1(n1)*G1P2*A2P1 - 3*G2(n2)*a1(n1)*
     & G1P1*A2P2 + 3*G2(n2)*a1(n1)*St*A2G1 - G2(n2)*a2(n1)*G1P2*A1P1 + 
     & G2(n2)*a2(n1)*G1P1*A1P2 + G2(n2)*a2(n1)*St*A1G1 - 2*G2(n2)*
     & a2(n1)*R*St*A1P2 + a1(n1)*a2(n2)*G1P2*G2P1 + 3*a1(n1)*a2(n2)*
     & G1P1*G2P2 - 3*a1(n1)*a2(n2)*St*G1G2 - 2*a1(n1)*a2(n2)*Q*St*G1P2
     &  - 3*a1(n2)*a2(n1)*G1P2*G2P1 - a1(n2)*a2(n1)*G1P1*G2P2 + a1(n2)*
     & a2(n1)*St*G1G2 - 2*a1(n2)*a2(n1)*Q*St*G1P2 + 3*DF(n1,n2)*G2P2*
     & A1P1*A2G1 + DF(n1,n2)*G2P2*A1G1*A2P1 + DF(n1,n2)*G2P1*A1P2*A2G1
     &  + DF(n1,n2)*G2P1*A1G1*A2P2 + DF(n1,n2)*G1P2*A1P1*A2G2 - DF(n1,
     & n2)*G1P2*A1G2*A2P1 - DF(n1,n2)*G1P2*G2P1*A1A2 - DF(n1,n2)*G1P1*
     & A1P2*A2G2 + 3*DF(n1,n2)*G1P1*A1G2*A2P2 - 3*DF(n1,n2)*G1P1*G2P2*
     & A1A2 + DF(n1,n2)*G1G2*A1P2*A2P1 - 3*DF(n1,n2)*G1G2*A1P1*A2P2 - 3
     & *DF(n1,n2)*St*A1G2*A2G1 - DF(n1,n2)*St*A1G1*A2G2 + 3*DF(n1,n2)*
     & St*G1G2*A1A2
      SUM42 = SUM42 - 2*DF(n1,n2)*Q*G1P2*A1P1*A2P2 + 2*DF(n1,n2)*Q*G1P1
     & *A1P2*A2P2 - 2*DF(n1,n2)*Q*St*A1P2*A2G1 + 2*DF(n1,n2)*Q*St*G1P2*
     & A1A2 - 2*DF(n1,n2)*R*G2P2*A1P2*A2P1 - 2*DF(n1,n2)*R*G2P1*A1P2*
     & A2P2 + 2*DF(n1,n2)*R*St*A1P2*A2G2 + 2*DF(n1,n2)*XM2*A1G1*A2G2 + 
     & 8*DF(n1,n2)*XM2*RQ*A1P2*A2P2 + 8*DF(n1,n2)*XM2*RQ*A1P2*A2P1 - 4*
     & DF(n1,n2)*XM2*Q*A1G1*A2P2 - 4*DF(n1,n2)*XM2*Q*A1G1*A2P1 - 4*DF(
     & n1,n2)*XM2*R*A1P2*A2G2
C
      SUM45 =
     & 3*P1(n2)*P2(n1)*A1G2*A2G1 + P1(n2)*P2(n1)*A1G1*A2G2 - 3*P1(n2)*
     & P2(n1)*G1G2*A1A2 - 2*P1(n2)*P2(n1)*Q*A1P1*A2G2 + 2*P1(n2)*P2(n1)
     & *R*A1P1*A2G1 - 2*P1(n2)*P2(n1)*R*G1P1*A1A2 - P1(n2)*G1(n1)*A1P2*
     & A2G2 + P1(n2)*G1(n1)*A1G2*A2P2 + P1(n2)*G1(n1)*G2P2*A1A2 - 2*
     & P1(n2)*G1(n1)*R*A1P1*A2P2 - 3*P1(n2)*G2(n1)*A1P2*A2G1 - P1(n2)*
     & G2(n1)*A1G1*A2P2 + 3*P1(n2)*G2(n1)*G1P2*A1A2 + 2*P1(n2)*G2(n1)*Q
     & *A1P1*A2P2 - P1(n2)*a1(n1)*G2P2*A2G1 + P1(n2)*a1(n1)*G1P2*A2G2
     &  - P1(n2)*a1(n1)*G1G2*A2P2 + 2*P1(n2)*a1(n1)*R*G1P1*A2P2 - 
     & P1(n2)*a2(n1)*G2P2*A1G1 - 3*P1(n2)*a2(n1)*G1P2*A1G2 + 3*P1(n2)*
     & a2(n1)*G1G2*A1P2 + 2*P1(n2)*a2(n1)*Q*G2P2*A1P1 - 2*P1(n2)*a2(n1)
     & *R*G1P2*A1P1 + 2*P1(n2)*a2(n1)*R*G1P1*A1P2 + P2(n1)*G1(n2)*A1P1*
     & A2G2 - 3*P2(n1)*G1(n2)*A1G2*A2P1 + 3*P2(n1)*G1(n2)*G2P1*A1A2 - 2
     & *P2(n1)*G1(n2)*R*A1P1*A2P1 - P2(n1)*G2(n2)*A1P1*A2G1 - P2(n1)*
     & G2(n2)*A1G1*A2P1 + P2(n1)*G2(n2)*G1P1*A1A2 + 2*P2(n1)*G2(n2)*Q*
     & A1P1*A2P1
      SUM45 = SUM45 - 3*P2(n1)*a1(n2)*G2P1*A2G1 - P2(n1)*a1(n2)*G1P1*
     & A2G2 + 3*P2(n1)*a1(n2)*G1G2*A2P1 + 2*P2(n1)*a1(n2)*R*G1P1*A2P1
     &  - P2(n1)*a2(n2)*G2P1*A1G1 + P2(n1)*a2(n2)*G1P1*A1G2 - P2(n1)*
     & a2(n2)*G1G2*A1P1 + 2*P2(n1)*a2(n2)*Q*G2P1*A1P1 + G1(n1)*G2(n2)*
     & A1P2*A2P1 - 3*G1(n1)*G2(n2)*A1P1*A2P2 - G1(n1)*G2(n2)*St*A1A2 - 
     & G1(n1)*a1(n2)*G2P2*A2P1 - G1(n1)*a1(n2)*G2P1*A2P2 + G1(n1)*
     & a1(n2)*St*A2G2 - 2*G1(n1)*a1(n2)*XM2*A2G2 + 4*G1(n1)*a1(n2)*XM2*
     & R*A2P2 + 4*G1(n1)*a1(n2)*XM2*R*A2P1 + 3*G1(n1)*a2(n2)*G2P2*A1P1
     &  + G1(n1)*a2(n2)*G2P1*A1P2 - G1(n1)*a2(n2)*St*A1G2 + 2*G1(n1)*
     & a2(n2)*R*St*A1P1 + 3*G1(n2)*G2(n1)*A1P2*A2P1 - G1(n2)*G2(n1)*
     & A1P1*A2P2 - 3*G1(n2)*G2(n1)*St*A1A2 + G1(n2)*a1(n1)*G2P2*A2P1 + 
     & G1(n2)*a1(n1)*G2P1*A2P2 - G1(n2)*a1(n1)*St*A2G2 + 2*G1(n2)*
     & a1(n1)*XM2*A2G2 - 4*G1(n2)*a1(n1)*XM2*R*A2P2 - 4*G1(n2)*a1(n1)*
     & XM2*R*A2P1 - G1(n2)*a2(n1)*G2P2*A1P1 - 3*G1(n2)*a2(n1)*G2P1*A1P2
     &  + 3*G1(n2)*a2(n1)*St*A1G2 + 2*G1(n2)*a2(n1)*R*St*A1P1 - 3*
     & G2(n1)*a1(n2)*G1P2*A2P1
      SUM45 = SUM45 + G2(n1)*a1(n2)*G1P1*A2P2 + 3*G2(n1)*a1(n2)*St*A2G1
     &  + G2(n1)*a2(n2)*G1P2*A1P1 - G2(n1)*a2(n2)*G1P1*A1P2 + G2(n1)*
     & a2(n2)*St*A1G1 - 2*G2(n1)*a2(n2)*Q*St*A1P1 - G2(n2)*a1(n1)*G1P2*
     & A2P1 + 3*G2(n2)*a1(n1)*G1P1*A2P2 + G2(n2)*a1(n1)*St*A2G1 + 
     & G2(n2)*a2(n1)*G1P2*A1P1 - G2(n2)*a2(n1)*G1P1*A1P2 + G2(n2)*
     & a2(n1)*St*A1G1 - 2*G2(n2)*a2(n1)*Q*St*A1P1 - a1(n1)*a2(n2)*G1P2*
     & G2P1 - 3*a1(n1)*a2(n2)*G1P1*G2P2 + a1(n1)*a2(n2)*St*G1G2 - 2*
     & a1(n1)*a2(n2)*R*St*G1P1 + 3*a1(n2)*a2(n1)*G1P2*G2P1 + a1(n2)*
     & a2(n1)*G1P1*G2P2 - 3*a1(n2)*a2(n1)*St*G1G2 - 2*a1(n2)*a2(n1)*R*
     & St*G1P1 + DF(n1,n2)*G2P2*A1P1*A2G1 + DF(n1,n2)*G2P2*A1G1*A2P1 + 
     & 3*DF(n1,n2)*G2P1*A1P2*A2G1 + DF(n1,n2)*G2P1*A1G1*A2P2 - DF(n1,n2
     & )*G1P2*A1P1*A2G2 + 3*DF(n1,n2)*G1P2*A1G2*A2P1 - 3*DF(n1,n2)*G1P2
     & *G2P1*A1A2 + DF(n1,n2)*G1P1*A1P2*A2G2 - DF(n1,n2)*G1P1*A1G2*A2P2
     &  - DF(n1,n2)*G1P1*G2P2*A1A2 - 3*DF(n1,n2)*G1G2*A1P2*A2P1 + DF(n1
     & ,n2)*G1G2*A1P1*A2P2 - 3*DF(n1,n2)*St*A1G2*A2G1 - DF(n1,n2)*St*
     & A1G1*A2G2
      SUM45 = SUM45 + 3*DF(n1,n2)*St*G1G2*A1A2 - 2*DF(n1,n2)*Q*G2P2*
     & A1P1*A2P1 - 2*DF(n1,n2)*Q*G2P1*A1P1*A2P2 + 2*DF(n1,n2)*Q*St*A1P1
     & *A2G2 + 2*DF(n1,n2)*R*G1P2*A1P1*A2P1 - 2*DF(n1,n2)*R*G1P1*A1P2*
     & A2P1 - 2*DF(n1,n2)*R*St*A1P1*A2G1 + 2*DF(n1,n2)*R*St*G1P1*A1A2
     &  + 2*DF(n1,n2)*XM2*A1G1*A2G2 + 8*DF(n1,n2)*XM2*RQ*A1P1*A2P2 + 8*
     & DF(n1,n2)*XM2*RQ*A1P1*A2P1 - 4*DF(n1,n2)*XM2*Q*A1P1*A2G2 - 4*DF(
     & n1,n2)*XM2*R*A1G1*A2P2 - 4*DF(n1,n2)*XM2*R*A1G1*A2P1
C
      SUM46 =
     &  - 3*P1(n2)*P2(n1)*A1G2*A2G1 - P1(n2)*P2(n1)*A1G1*A2G2 + 3*
     & P1(n2)*P2(n1)*G1G2*A1A2 + 2*P1(n2)*P2(n1)*Q*A1P2*A2G2 - 2*P1(n2)
     & *P2(n1)*R*A1P2*A2G1 + 2*P1(n2)*P2(n1)*R*G1P2*A1A2 - P1(n2)*
     & G1(n1)*A1P2*A2G2 + 3*P1(n2)*G1(n1)*A1G2*A2P2 - 3*P1(n2)*G1(n1)*
     & G2P2*A1A2 + 2*P1(n2)*G1(n1)*R*A1P2*A2P2 + P1(n2)*G2(n1)*A1P2*
     & A2G1 + P1(n2)*G2(n1)*A1G1*A2P2 - P1(n2)*G2(n1)*G1P2*A1A2 - 2*
     & P1(n2)*G2(n1)*Q*A1P2*A2P2 + 3*P1(n2)*a1(n1)*G2P2*A2G1 + P1(n2)*
     & a1(n1)*G1P2*A2G2 - 3*P1(n2)*a1(n1)*G1G2*A2P2 - 2*P1(n2)*a1(n1)*R
     & *G1P2*A2P2 + P1(n2)*a2(n1)*G2P2*A1G1 - P1(n2)*a2(n1)*G1P2*A1G2
     &  + P1(n2)*a2(n1)*G1G2*A1P2 - 2*P1(n2)*a2(n1)*Q*G2P2*A1P2 + 
     & P2(n1)*G1(n2)*A1P1*A2G2 - P2(n1)*G1(n2)*A1G2*A2P1 - P2(n1)*
     & G1(n2)*G2P1*A1A2 + 2*P2(n1)*G1(n2)*R*A1P2*A2P1 + 3*P2(n1)*G2(n2)
     & *A1P1*A2G1 + P2(n1)*G2(n2)*A1G1*A2P1 - 3*P2(n1)*G2(n2)*G1P1*A1A2
     &  - 2*P2(n1)*G2(n2)*Q*A1P2*A2P1 + P2(n1)*a1(n2)*G2P1*A2G1 - 
     & P2(n1)*a1(n2)*G1P1*A2G2
      SUM46 = SUM46 + P2(n1)*a1(n2)*G1G2*A2P1 - 2*P2(n1)*a1(n2)*R*G1P2*
     & A2P1 + P2(n1)*a2(n2)*G2P1*A1G1 + 3*P2(n1)*a2(n2)*G1P1*A1G2 - 3*
     & P2(n1)*a2(n2)*G1G2*A1P1 - 2*P2(n1)*a2(n2)*Q*G2P1*A1P2 - 2*P2(n1)
     & *a2(n2)*R*G1P2*A1P1 + 2*P2(n1)*a2(n2)*R*G1P1*A1P2 + G1(n1)*
     & G2(n2)*A1P2*A2P1 - 3*G1(n1)*G2(n2)*A1P1*A2P2 + 3*G1(n1)*G2(n2)*
     & St*A1A2 - G1(n1)*a1(n2)*G2P2*A2P1 - G1(n1)*a1(n2)*G2P1*A2P2 + 
     & G1(n1)*a1(n2)*St*A2G2 - 2*G1(n1)*a1(n2)*XM2*A2G2 + 4*G1(n1)*
     & a1(n2)*XM2*R*A2P2 + 4*G1(n1)*a1(n2)*XM2*R*A2P1 + 3*G1(n1)*a2(n2)
     & *G2P2*A1P1 + G1(n1)*a2(n2)*G2P1*A1P2 - 3*G1(n1)*a2(n2)*St*A1G2
     &  - 2*G1(n1)*a2(n2)*R*St*A1P2 + 3*G1(n2)*G2(n1)*A1P2*A2P1 - 
     & G1(n2)*G2(n1)*A1P1*A2P2 + G1(n2)*G2(n1)*St*A1A2 + G1(n2)*a1(n1)*
     & G2P2*A2P1 + G1(n2)*a1(n1)*G2P1*A2P2 - G1(n2)*a1(n1)*St*A2G2 + 2*
     & G1(n2)*a1(n1)*XM2*A2G2 - 4*G1(n2)*a1(n1)*XM2*R*A2P2 - 4*G1(n2)*
     & a1(n1)*XM2*R*A2P1 - G1(n2)*a2(n1)*G2P2*A1P1 - 3*G1(n2)*a2(n1)*
     & G2P1*A1P2 + G1(n2)*a2(n1)*St*A1G2 - 2*G1(n2)*a2(n1)*R*St*A1P2 - 
     & 3*G2(n1)*a1(n2)*G1P2*A2P1
      SUM46 = SUM46 + G2(n1)*a1(n2)*G1P1*A2P2 - G2(n1)*a1(n2)*St*A2G1
     &  + G2(n1)*a2(n2)*G1P2*A1P1 - G2(n1)*a2(n2)*G1P1*A1P2 - G2(n1)*
     & a2(n2)*St*A1G1 + 2*G2(n1)*a2(n2)*Q*St*A1P2 - G2(n2)*a1(n1)*G1P2*
     & A2P1 + 3*G2(n2)*a1(n1)*G1P1*A2P2 - 3*G2(n2)*a1(n1)*St*A2G1 + 
     & G2(n2)*a2(n1)*G1P2*A1P1 - G2(n2)*a2(n1)*G1P1*A1P2 - G2(n2)*
     & a2(n1)*St*A1G1 + 2*G2(n2)*a2(n1)*Q*St*A1P2 - a1(n1)*a2(n2)*G1P2*
     & G2P1 - 3*a1(n1)*a2(n2)*G1P1*G2P2 + 3*a1(n1)*a2(n2)*St*G1G2 + 2*
     & a1(n1)*a2(n2)*R*St*G1P2 + 3*a1(n2)*a2(n1)*G1P2*G2P1 + a1(n2)*
     & a2(n1)*G1P1*G2P2 - a1(n2)*a2(n1)*St*G1G2 + 2*a1(n2)*a2(n1)*R*St*
     & G1P2 - 3*DF(n1,n2)*G2P2*A1P1*A2G1 - DF(n1,n2)*G2P2*A1G1*A2P1 - 
     & DF(n1,n2)*G2P1*A1P2*A2G1 - DF(n1,n2)*G2P1*A1G1*A2P2 - DF(n1,n2)*
     & G1P2*A1P1*A2G2 + DF(n1,n2)*G1P2*A1G2*A2P1 + DF(n1,n2)*G1P2*G2P1*
     & A1A2 + DF(n1,n2)*G1P1*A1P2*A2G2 - 3*DF(n1,n2)*G1P1*A1G2*A2P2 + 3
     & *DF(n1,n2)*G1P1*G2P2*A1A2 - DF(n1,n2)*G1G2*A1P2*A2P1 + 3*DF(n1,
     & n2)*G1G2*A1P1*A2P2 + 3*DF(n1,n2)*St*A1G2*A2G1 + DF(n1,n2)*St*
     & A1G1*A2G2
      SUM46 = SUM46 - 3*DF(n1,n2)*St*G1G2*A1A2 + 2*DF(n1,n2)*Q*G2P2*
     & A1P2*A2P1 + 2*DF(n1,n2)*Q*G2P1*A1P2*A2P2 - 2*DF(n1,n2)*Q*St*A1P2
     & *A2G2 + 2*DF(n1,n2)*R*G1P2*A1P1*A2P2 - 2*DF(n1,n2)*R*G1P1*A1P2*
     & A2P2 + 2*DF(n1,n2)*R*St*A1P2*A2G1 - 2*DF(n1,n2)*R*St*G1P2*A1A2
     &  - 2*DF(n1,n2)*XM2*A1G1*A2G2 - 8*DF(n1,n2)*XM2*RQ*A1P2*A2P2 - 8*
     & DF(n1,n2)*XM2*RQ*A1P2*A2P1 + 4*DF(n1,n2)*XM2*Q*A1P2*A2G2 + 4*DF(
     & n1,n2)*XM2*R*A1G1*A2P2 + 4*DF(n1,n2)*XM2*R*A1G1*A2P1
C
      SUM62 =
     & 3*P1(n2)*P2(n1)*A1G2*A2G1 - 3*P1(n2)*P2(n1)*A1G1*A2G1 - 3*P1(n2)
     & *P2(n1)*G1G2*A1A2 + 2*P1(n2)*P2(n1)*W1*A1A2 - 2*P1(n2)*P2(n1)*R2
     & *St*A1A2 + 2*P1(n2)*P2(n1)*Q*A1P2*A2G1 + 2*P1(n2)*P2(n1)*Q*A1P1*
     & A2G1 - P1(n2)*P2(n1)*Q*A1G1*A2P2 - P1(n2)*P2(n1)*Q*A1G1*A2P1 - 
     & P1(n2)*P2(n1)*Q*G1P2*A1A2 - P1(n2)*P2(n1)*Q*G1P1*A1A2 + 2*P1(n2)
     & *P2(n1)*R*A1P2*A2G1 + 2*P1(n2)*P2(n1)*R*A1P1*A2G1 + P1(n2)*
     & P2(n1)*R*A1G2*A2P2 + P1(n2)*P2(n1)*R*A1G2*A2P1 - P1(n2)*P2(n1)*R
     & *G2P2*A1A2 - P1(n2)*P2(n1)*R*G2P1*A1A2 + 2*P1(n2)*G1(n1)*A1P2*
     & A2G1 - 2*P1(n2)*G1(n1)*A1G2*A2P2 + P1(n2)*G1(n1)*A1G1*A2P2 + 2*
     & P1(n2)*G1(n1)*G2P2*A1A2 - 4*P1(n2)*G1(n1)*G1P2*A1A2 - P1(n2)*
     & G1(n1)*Q*A1P2*A2P2 + P1(n2)*G1(n1)*Q*A1P2*A2P1 - 2*P1(n2)*G1(n1)
     & *Q*A1P1*A2P2 - 2*P1(n2)*G1(n1)*R*A1P2*A2P2 - 2*P1(n2)*G1(n1)*R*
     & A1P1*A2P2 - P1(n2)*G2(n1)*A1P2*A2G1 + 2*P1(n2)*G2(n1)*G1P2*A1A2
     &  - P1(n2)*G2(n1)*R*A1P2*A2P2 - P1(n2)*G2(n1)*R*A1P2*A2P1 - 
     & P1(n2)*a1(n1)*G2P2*A2G1
      SUM62 = SUM62 + 2*P1(n2)*a1(n1)*G1P2*A2G1 + P1(n2)*a1(n1)*G1G2*
     & A2P2 + 2*P1(n2)*a1(n1)*R2*St*A2P2 + 2*P1(n2)*a1(n1)*Q*G1P2*A2P2
     &  + P1(n2)*a1(n1)*Q*G1P2*A2P1 + P1(n2)*a1(n1)*Q*G1P1*A2P2 - 
     & P1(n2)*a1(n1)*R*G2P2*A2P1 + P1(n2)*a1(n1)*R*G2P1*A2P2 - 2*P1(n2)
     & *a2(n1)*G1P2*A1G2 + P1(n2)*a2(n1)*G1P2*A1G1 + P1(n2)*a2(n1)*G1G2
     & *A1P2 + 2*P1(n2)*a2(n1)*R2*St*A1P2 - P1(n2)*a2(n1)*Q*G1P2*A1P2
     &  - 2*P1(n2)*a2(n1)*Q*G1P2*A1P1 + P1(n2)*a2(n1)*Q*G1P1*A1P2 + 
     & P1(n2)*a2(n1)*R*G2P2*A1P2 + P1(n2)*a2(n1)*R*G2P1*A1P2 - 2*P1(n2)
     & *a2(n1)*R*G1P2*A1P2 - 2*P1(n2)*a2(n1)*R*G1P2*A1P1 + 2*P2(n1)*
     & G1(n2)*A1P1*A2G1 - 2*P2(n1)*G1(n2)*A1G2*A2P1 + P2(n1)*G1(n2)*
     & A1G1*A2P1 + 2*P2(n1)*G1(n2)*G2P1*A1A2 - 4*P2(n1)*G1(n2)*G1P1*
     & A1A2 - 2*P2(n1)*G1(n2)*Q*A1P2*A2P1 + P2(n1)*G1(n2)*Q*A1P1*A2P2
     &  - P2(n1)*G1(n2)*Q*A1P1*A2P1 - 2*P2(n1)*G1(n2)*R*A1P2*A2P1 - 2*
     & P2(n1)*G1(n2)*R*A1P1*A2P1 - P2(n1)*G2(n2)*A1P1*A2G1 + 2*P2(n1)*
     & G2(n2)*G1P1*A1A2 - P2(n1)*G2(n2)*R*A1P1*A2P2 - P2(n1)*G2(n2)*R*
     & A1P1*A2P1
      SUM62 = SUM62 - P2(n1)*a1(n2)*G2P1*A2G1 + 2*P2(n1)*a1(n2)*G1P1*
     & A2G1 + P2(n1)*a1(n2)*G1G2*A2P1 + 2*P2(n1)*a1(n2)*R2*St*A2P1 + 
     & P2(n1)*a1(n2)*Q*G1P2*A2P1 + P2(n1)*a1(n2)*Q*G1P1*A2P2 + 2*P2(n1)
     & *a1(n2)*Q*G1P1*A2P1 + P2(n1)*a1(n2)*R*G2P2*A2P1 - P2(n1)*a1(n2)*
     & R*G2P1*A2P2 - 2*P2(n1)*a2(n2)*G1P1*A1G2 + P2(n1)*a2(n2)*G1P1*
     & A1G1 + P2(n1)*a2(n2)*G1G2*A1P1 + 2*P2(n1)*a2(n2)*R2*St*A1P1 + 
     & P2(n1)*a2(n2)*Q*G1P2*A1P1 - 2*P2(n1)*a2(n2)*Q*G1P1*A1P2 - P2(n1)
     & *a2(n2)*Q*G1P1*A1P1 + P2(n1)*a2(n2)*R*G2P2*A1P1 + P2(n1)*a2(n2)*
     & R*G2P1*A1P1 - 2*P2(n1)*a2(n2)*R*G1P1*A1P2 - 2*P2(n1)*a2(n2)*R*
     & G1P1*A1P1 + 4*G1(n1)*G1(n2)*St*A1A2 - 2*G1(n1)*G2(n2)*St*A1A2 - 
     & 2*G1(n1)*a1(n2)*St*A2G1 - G1(n1)*a1(n2)*Q*St*A2P2 - G1(n1)*
     & a1(n2)*Q*St*A2P1 + 2*G1(n1)*a2(n2)*St*A1G2 - G1(n1)*a2(n2)*St*
     & A1G1 + 2*G1(n1)*a2(n2)*Q*St*A1P2 + 2*G1(n1)*a2(n2)*Q*St*A1P1 + 2
     & *G1(n1)*a2(n2)*R*St*A1P2 + 2*G1(n1)*a2(n2)*R*St*A1P1 - 2*G1(n2)*
     & G2(n1)*St*A1A2 - 2*G1(n2)*a1(n1)*St*A2G1 - G1(n2)*a1(n1)*Q*St*
     & A2P2
      SUM62 = SUM62 - G1(n2)*a1(n1)*Q*St*A2P1 + 2*G1(n2)*a2(n1)*St*A1G2
     &  - G1(n2)*a2(n1)*St*A1G1 + 2*G1(n2)*a2(n1)*Q*St*A1P2 + 2*G1(n2)*
     & a2(n1)*Q*St*A1P1 + 2*G1(n2)*a2(n1)*R*St*A1P2 + 2*G1(n2)*a2(n1)*R
     & *St*A1P1 + G2(n1)*a1(n2)*St*A2G1 + G2(n1)*a1(n2)*R*St*A2P2 + 
     & G2(n1)*a1(n2)*R*St*A2P1 + G2(n2)*a1(n1)*St*A2G1 + G2(n2)*a1(n1)*
     & R*St*A2P2 + G2(n2)*a1(n1)*R*St*A2P1 - a1(n1)*a2(n2)*St*G1G2 - 2*
     & a1(n1)*a2(n2)*R2*St**2 - a1(n1)*a2(n2)*Q*St*G1P2 - a1(n1)*a2(n2)
     & *Q*St*G1P1 - a1(n1)*a2(n2)*R*St*G2P2 - a1(n1)*a2(n2)*R*St*G2P1
     &  - a1(n2)*a2(n1)*St*G1G2 - 2*a1(n2)*a2(n1)*R2*St**2 - a1(n2)*
     & a2(n1)*Q*St*G1P2 - a1(n2)*a2(n1)*Q*St*G1P1 - a1(n2)*a2(n1)*R*St*
     & G2P2 - a1(n2)*a2(n1)*R*St*G2P1 + DF(n1,n2)*G2P2*A1P1*A2G1 + DF(
     & n1,n2)*G2P1*A1P2*A2G1 - 2*DF(n1,n2)*G1P2*A1P1*A2G1 + 2*DF(n1,n2)
     & *G1P2*A1G2*A2P1 - DF(n1,n2)*G1P2*A1G1*A2P1 - 2*DF(n1,n2)*G1P2*
     & G2P1*A1A2 - 2*DF(n1,n2)*G1P1*A1P2*A2G1 + 2*DF(n1,n2)*G1P1*A1G2*
     & A2P2 - DF(n1,n2)*G1P1*A1G1*A2P2 - 2*DF(n1,n2)*G1P1*G2P2*A1A2 + 4
     & *DF(n1,n2)*G1P1*G1P2*A1A2
      SUM62 = SUM62 - DF(n1,n2)*G1G2*A1P2*A2P1 - DF(n1,n2)*G1G2*A1P1*
     & A2P2 - 3*DF(n1,n2)*St*A1G2*A2G1 + 3*DF(n1,n2)*St*A1G1*A2G1 + 3*
     & DF(n1,n2)*St*G1G2*A1A2 - 2*DF(n1,n2)*St*W1*A1A2 - 2*DF(n1,n2)*R2
     & *St*A1P2*A2P1 - 2*DF(n1,n2)*R2*St*A1P1*A2P2 + 2*DF(n1,n2)*R2*
     & St**2*A1A2 + DF(n1,n2)*Q*G1P2*A1P2*A2P1 - 2*DF(n1,n2)*Q*G1P2*
     & A1P1*A2P2 + DF(n1,n2)*Q*G1P2*A1P1*A2P1 + DF(n1,n2)*Q*G1P1*A1P2*
     & A2P2 - 2*DF(n1,n2)*Q*G1P1*A1P2*A2P1 + DF(n1,n2)*Q*G1P1*A1P1*A2P2
     &  - 2*DF(n1,n2)*Q*St*A1P2*A2G1 - 2*DF(n1,n2)*Q*St*A1P1*A2G1 + DF(
     & n1,n2)*Q*St*A1G1*A2P2 + DF(n1,n2)*Q*St*A1G1*A2P1 + DF(n1,n2)*Q*
     & St*G1P2*A1A2 + DF(n1,n2)*Q*St*G1P1*A1A2 - DF(n1,n2)*R*G2P2*A1P2*
     & A2P1 + DF(n1,n2)*R*G2P2*A1P1*A2P1 + DF(n1,n2)*R*G2P1*A1P2*A2P2
     &  - DF(n1,n2)*R*G2P1*A1P1*A2P2 + 2*DF(n1,n2)*R*G1P2*A1P2*A2P1 + 2
     & *DF(n1,n2)*R*G1P2*A1P1*A2P1 + 2*DF(n1,n2)*R*G1P1*A1P2*A2P2 + 2*
     & DF(n1,n2)*R*G1P1*A1P1*A2P2 - 2*DF(n1,n2)*R*St*A1P2*A2G1 - 2*DF(
     & n1,n2)*R*St*A1P1*A2G1 - DF(n1,n2)*R*St*A1G2*A2P2 - DF(n1,n2)*R*
     & St*A1G2*A2P1
      SUM62 = SUM62 + DF(n1,n2)*R*St*G2P2*A1A2 + DF(n1,n2)*R*St*G2P1*
     & A1A2 - 2*DF(n1,n2)*XM2*A1G1*A2G1 + 4*DF(n1,n2)*XM2*RQ*A1P2*A2P2
     &  + 4*DF(n1,n2)*XM2*RQ*A1P2*A2P1 + 4*DF(n1,n2)*XM2*RQ*A1P1*A2P2
     &  + 4*DF(n1,n2)*XM2*RQ*A1P1*A2P1 - 4*DF(n1,n2)*XM2*R2*A1P2*A2P2
     &  - 4*DF(n1,n2)*XM2*R2*A1P2*A2P1 - 4*DF(n1,n2)*XM2*R2*A1P1*A2P2
     &  - 4*DF(n1,n2)*XM2*R2*A1P1*A2P1 - 2*DF(n1,n2)*XM2*Q*A1G1*A2P2 - 
     & 2*DF(n1,n2)*XM2*Q*A1G1*A2P1 + 4*DF(n1,n2)*XM2*R*A1P2*A2G1 + 4*
     & DF(n1,n2)*XM2*R*A1P1*A2G1 + 2*DF(n1,n2)*XM2*R*A1G1*A2P2 + 2*DF(
     & n1,n2)*XM2*R*A1G1*A2P1
C
      SUM64B =
     &  - P1(n2)*P2(n1)*A1A2 + P1(n2)*a1(n1)*A2P2 + P1(n2)*a2(n1)*A1P2
     &  + P2(n1)*a1(n2)*A2P1 + P2(n1)*a2(n2)*A1P1 - a1(n1)*a2(n2)*St - 
     & a1(n2)*a2(n1)*St - DF(n1,n2)*A1P2*A2P1 - DF(n1,n2)*A1P1*A2P2 + 
     & DF(n1,n2)*St*A1A2 - 2*DF(n1,n2)*XM2*A1A2
      SUM64C =
     & P1(n2)*P2(n1)*A1A2 - P1(n2)*a1(n1)*A2P2 - P1(n2)*a2(n1)*A1P2 - 
     & P2(n1)*a1(n2)*A2P1 - P2(n1)*a2(n2)*A1P1 + a1(n1)*a2(n2)*St + 
     & a1(n2)*a2(n1)*St + DF(n1,n2)*A1P2*A2P1 + DF(n1,n2)*A1P1*A2P2 - 
     & DF(n1,n2)*St*A1A2 + 2*DF(n1,n2)*XM2*A1A2
C
      AMP(n1,n2) = AMP(n1,n2)
     . + SUM41/Z41 *( 6.)  + SUM42/Z42 *(-6.)
     . + SUM45/Z45 *(-6.)  + SUM46/Z46 *( 6.)
c    . + SUM62/Z62 *(-12.) +(SUM64B *12. -SUM64C *12.)/Z64/2.
     . + SUM62/Z62 *(-12.) +(SUM64B *12.)/Z64
   88 CONTINUE
C
   99 RETURN
      END

      SUBROUTINE SCALAR1(G1,A1, G2,A2, P1,P2)
C*******************************************************************C
C      Dot-products of particle momenta                             C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G1(4),A1(4),G2(4),A2(4),P1(4),P2(4)
      COMMON/DOTPR1/W1,W2,G1G2,G1P1,G1P2,G2P1,G2P2,P1P2
     .       A1A2,A1G1,A1G2,A1P1,A1P2,A2G1,A2G2,A2P1,A2P2
      COMMON/QUARK1/XC,XC2,XB,XB2,XM,XM2,CHC,CHB
      COMMON/OTHER1/XM4,R,Q,R2,Q2,RQ,St,H1,H2
      SAVE
      W1   =DOT_S(G1,G1)
      W2   =DOT_S(G2,G2)
      G1G2 =DOT_S(G1,G2)
      G1P1 =DOT_S(G1,P1)
      G1P2 =DOT_S(G1,P2)
      G2P1 =DOT_S(G2,P1)
      G2P2 =DOT_S(G2,P2)
      P1P2 =DOT_S(P1,P2)
      A1A2 =DOT_S(A1,A2)
      A1G1 =DOT_S(A1,G1)
      A1G2 =DOT_S(A1,G2)
      A1P1 =DOT_S(A1,P1)
      A1P2 =DOT_S(A1,P2)
      A2G1 =DOT_S(A2,G1)
      A2G2 =DOT_S(A2,G2)
      A2P1 =DOT_S(A2,P1)
      A2P2 =DOT_S(A2,P2)
      XM4=XM2*XM2
      R = XC/XM
      Q = XB/XM
      R2=R*R
      Q2=Q*Q
      RQ=R*Q
      St=P1P2+XM2
      H1=W1/2.
      H2=W2/2.
      RETURN
      END
