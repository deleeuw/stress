      SUBROUTINE BSPLINE (N, X, M, Y, NORD, NORM, G, H)
     
      DOUBLE PRECISION X(*), Y(*), G(*), H(*)
      DOUBLE PRECISION V(NORD), Z(NORD)

      INTEGER I, J, K, N, M, NORD, NORM, JINT

      JINT=0

      DO 30 I=1,N
      CALL TG04AD (M, Y, NORD, NORM, X(I), JINT, V, Z)
      DO 20 J=1,NORD
      JJ = JINT - NORD + J
      K = I + (JJ - 1) * N
      G(K) = V(J)
   20 H(K) = Z(J)
   30 CONTINUE
 
      END
     
