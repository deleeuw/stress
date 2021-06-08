* COPYRIGHT (c) 1976 AEA Technology
*######DATE 14 Jan 1993
C       Toolpack tool decs employed.
C       Reference TG04C removed.
C       JINT undefined at label 80, code modified
C       SAVE statements added.
C       ZERO and ONE made PARAMETER.
C
C  EAT 21/6/93 EXTERNAL statement put in for block data so will work on VAXs.
C
C
      SUBROUTINE TG04A(N,X,K,NORM,XVALUE,JINT,V,VINT)
C      ********************************************************
C                            PURPOSE
C      ********************************************************
C
C      GIVEN THE DATA POINTS X(1)<=X(2)<=,...,<=X(N-1)<=X(N)
C      AND A VALUE OF X, X=XVALUE, WHERE X(1)<=XVALUE<=X(N),
C      THIS SUBROUTINE COMPUTES THE VALUES OF THE B-SPLINES OF
C      DEGREE K-1, (MKI(X)) 1<=I<=N-K, AT X=XVALUE, AND THE
C      VALUES OF THE CORRESPONDING INTEGRALS OVER THE RANGE
C      X(I)<=X<=XVALUE.
C      N.B. THE FUNCTION MKI(X) IS A SPLINE OF DEGREE K-1 WITH
C           KNOTS AT THE POINTS X(I),X(I+1),...,X(I+K), AND IS
C           NON-ZERO ONLY OVER THE RANGE X(I)<X<X(I+K).
C
C      ********************************************************
C     .. Parameters ..
      REAL ZERO,ONE
      PARAMETER (ZERO=0.0,ONE=1.0)
C     ..
C     .. Scalar Arguments ..
      REAL XVALUE
      INTEGER JINT,K,N,NORM
C     ..
C     .. Array Arguments ..
      REAL V(*),VINT(*),X(*)
C     ..
C     .. Local Scalars ..
      REAL E1,E2,E3,XLEFT,XRITE
      INTEGER I,I1,I2,IL,IPK,IR,IRMI,IRP1,J,JINTMJ,JINTMK,JINTPJ,JINTPL,
     +        JJ1,JTPJMK,KP1,L,MIDDLE,MJ,NK,NKJ
C     ..
C     .. Data block external statement
      EXTERNAL TG04C
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC FLOAT,MAX0,MIN0
C     ..
C     .. Common blocks ..
      COMMON /TG04B/LP,IFAIL
      INTEGER IFAIL,LP
C     ..
C     .. Save statement ..
      SAVE /TG04B/
C     ..
C     .. Executable Statements ..
      IFAIL = 0
C      ********************************************************
C      CHECK THAT N,K, AND NORM HAVE SENSIBLE VALUES.
C      ********************************************************
      IF (N.GE.K+1) GO TO 10
      IFAIL = 1
      GO TO 330

   10 IF (K.GE.1) GO TO 20
      IFAIL = 2
      GO TO 330

   20 IF (NORM.GE.1 .AND. NORM.LE.2) GO TO 30
      NORM = 1
      IFAIL = 5
      IF (LP.EQ.0) GO TO 30
      WRITE (LP,FMT=380)
C      N,K, AND NORM HAVE SENSIBLE VALUES.
C      ********************************************************
   30 XLEFT = X(1)
      XRITE = X(N)
C      ********************************************************
C      INITIALISING V AND VINT TO ZERO.
C      ********************************************************
      DO 40 I = 1,K
        VINT(I) = ZERO
   40 V(I) = ZERO
C      END OF INITIALISATION.
C      ********************************************************
C
      IF (XVALUE.GE.XLEFT) GO TO 50
      JINT = 1
      RETURN
C
C      **********************************************************
C      CHECK, AND IF NECESSARY COMPUTE, THE VALUE OF JINT.
C      **********************************************************
   50 IF (XVALUE.LT.XRITE) GO TO 80
C      CHECK TO SEE IF XVALUE(=XRITE) IS A MULTIPLE KNOT
      I1 = N
      I2 = N - 1
      DO 60 I = 1,I2
        IRMI = N - I
        IF (X(IRMI).LT.XVALUE) GO TO 70
        I1 = IRMI - 1
   60 CONTINUE
   70 JINT = I1
      IF (JINT.GT.0) GO TO 130
      GO TO 180
C
C   80 IF (JINT.LT.1 .OR. JINT.GE.N) GO TO 90
C      IF (XVALUE.LT.X(JINT)) GO TO 90
C      IF (XVALUE.LT.X(JINT+1)) GO TO 130
C   90 IL = 1
   80 IL = 1
      IR = N
  100 IF (IR-IL.LE.1) GO TO 120
      MIDDLE = (IL+IR)/2
      IF (XVALUE.LT.X(MIDDLE)) GO TO 110
      IL = MIDDLE
      GO TO 100

  110 IR = MIDDLE
      GO TO 100

  120 JINT = IL
C      **********************************************************
C
C      **********************************************************
C      CHECK THAT THE DATA POINTS USED IN THE CALCULATION OF THE
C      B-SPLINES ARE IN ASCENDING ORDER
C      **********************************************************
  130 JINTMK = JINT - K
      IR = MIN0(N,JINT+K) - 1
      IF (JINT.EQ.N) IL = MIN0(IR,MAX0(1,JINTMK))
      IF (JINT.LT.N) IL = MIN0(IR,MAX0(1,JINTMK+1))
      DO 140 I = IL,IR
        IF (X(I).GT.X(I+1)) GO TO 150
  140 CONTINUE
      GO TO 160

  150 IFAIL = 3
      GO TO 330
C      **********************************************************
C
C      **********************************************************
C      CHECK THAT NO MORE THAN K OF THESE DATA POINTS COALESCE
C      **********************************************************
  160 I = IL
      IRP1 = IR + 1
  170 IPK = I + K
      IF (IPK.GT.IRP1) GO TO 190
      IF (X(I).GE.X(IPK)) GO TO 180
      I = I + 1
      GO TO 170

  180 IFAIL = 4
      GO TO 330
C      ********************************************************
  190 JJ1 = N - JINT
      IF (JJ1.GT.0) GO TO 200
      IF (NORM.EQ.1) VINT(1) = ONE/FLOAT(K)
      IF (NORM.EQ.2) VINT(1) = (X(JINT)-X(JINTMK))/FLOAT(K)
      RETURN
C
C      ********************************************************
C      COMPUTE M1JINT(XVALUE) AND (XVALUE-X(JINT))*M1JINT(XVALUE)
C      AND TEST FOR K=1
C      ********************************************************
  200 E1 = XVALUE - X(JINT)
      E2 = X(JINT+1) - XVALUE
      V(1) = ONE/ (X(JINT+1)-X(JINT))
      VINT(1) = E1*V(1)
      IF (K.GT.1) GO TO 210
      IF (NORM.EQ.1) RETURN
      V(1) = ONE
      VINT(1) = E1
      RETURN
C      ********************************************************
C
C      ********************************************************
C      COMPUTE AND STORE IN V(J) THE VALUE OF MJJINT(XVALUE)
C      AND STORE IN VINT(J) THE CORRESPONDING VALUE OF
C      (XVALUE-X(JINT))*MJJINT(XVALUE) FOR J=2,...,NK
C      ********************************************************
  210 NK = MIN0(K,JJ1)
      IF (NK.EQ.1) GO TO 240
      DO 230 J = 2,NK
        JINTPJ = JINT + J
        IF (J.EQ.K .AND. NORM.EQ.2) GO TO 220
        V(J) = E1*V(J-1)/ (X(JINTPJ)-X(JINT))
        VINT(J) = E1*V(J)
        GO TO 230

  220   V(J) = E1*V(J-1)
        VINT(J) = E1*V(J)/ (X(JINTPJ)-X(JINT))
  230 CONTINUE
C      ********************************************************
C
  240 IF (JINT.EQ.1) GO TO 300
C
C      ********************************************************
C      COMPUTE AND STORE IN V(1) THE VALUE OF
C      MJ+1JINT-J(XVALUE) AND STORE IN VINT(1) THE APPROPRIATE
C      SUM FOR THE INTEGRAL
C      ********************************************************
      MJ = MIN0(JINT,K) - 1
      DO 290 J = 1,MJ
        JINTMJ = JINT - J
        E3 = XVALUE - X(JINTMJ)
        IF (J+1.EQ.K .AND. NORM.EQ.2) GO TO 250
        V(1) = E2*V(1)/ (X(JINT+1)-X(JINTMJ))
        VINT(1) = VINT(1) + E3*V(1)
        GO TO 260

  250   V(1) = E2*V(1)
        VINT(1) = VINT(1) + E3*V(1)/ (X(JINT+1)-X(JINTMJ))
C      ********************************************************
C
  260   NKJ = MIN0(K-J,JJ1)
        IF (NKJ.LE.1) GO TO 290
C
C      ********************************************************
C      COMPUTE AND STORE IN V(L) THE VALUE MJ+LJINT-J(XVALUE)
C      AND STORE IN VINT(L) THE ACCUMULATED SUM FOR THE
C      INTEGRALS, FOR L=2,...,NKJ=MIN(K-J,JJ1)
C      *********************************************************
        DO 280 L = 2,NKJ
          JINTPL = JINT + L
          IF (J+L.EQ.K .AND. NORM.EQ.2) GO TO 270
          V(L) = (E3*V(L-1)+ (X(JINTPL)-XVALUE)*V(L))/
     +           (X(JINTPL)-X(JINTMJ))
          VINT(L) = VINT(L) + E3*V(L)
          GO TO 280

  270     V(L) = E3*V(L-1) + (X(JINTPL)-XVALUE)*V(L)
          VINT(L) = VINT(L) + E3*V(L)/ (X(JINTPL)-X(JINTMJ))
  280   CONTINUE
  290 CONTINUE
C      *********************************************************
C
C      ********************************************************
C      FORMING THE INTEGRALS AND RESETTING THE APPROPRIATE
C      LOCATIONS OF V AND VINT TO ZERO
C      ********************************************************
  300 DO 320 J = 1,K
        JINTPJ = JINT + J
        JTPJMK = JINTPJ - K
        IF (JINTPJ.GT.N) RETURN
        IF (JTPJMK.LT.1) GO TO 310
        IF (NORM.EQ.1) VINT(J) = VINT(J)/FLOAT(K)
        IF (NORM.EQ.2) VINT(J) = (X(JINTPJ)-X(JTPJMK))*VINT(J)/FLOAT(K)
        GO TO 320

  310   V(J) = ZERO
        VINT(J) = ZERO
  320 CONTINUE
      RETURN
C      *********************************************************
C
C      DIAGNOSTIC PRINTING
  330 IF (LP.EQ.0) RETURN
      GO TO (340,350,360,370) IFAIL

  340 KP1 = K + 1
      WRITE (LP,FMT=400) N,KP1
      RETURN

  350 WRITE (LP,FMT=390) K
      RETURN

  360 WRITE (LP,FMT=410)
      RETURN

  370 WRITE (LP,FMT=420) K
      RETURN

  380 FORMAT (1X,18HMESSAGE FROM TG04A,3X,
     +       41HNORM WAS .NE.1 OR 2 AND HAS BEEN SET TO 1)
  390 FORMAT (1X,18HMESSAGE FROM TG04A,3X,2HK=,I4,18H AND IS NOT .GT. 1)
  400 FORMAT (1X,18HMESSAGE FROM TG04A,3X,2HN=,I4,
     +       20H AND IS NOT .GE.K+1=,I4)
  410 FORMAT (1X,18HMESSAGE FROM TG04A,3X,
     +       43HTHE DATA POINTS ARE NOT IN ASCENDING ORDER.)
  420 FORMAT (1X,18HMESSAGE FROM TG04A,3X,10HMORE THAN ,I4,
     +       22H DATA POINTS COALESCE.)

      END
      BLOCK DATA TG04C
C     .. Common blocks ..
      COMMON /TG04B/LP,IFAIL
      INTEGER IFAIL,LP
C     ..
C     .. Save statement ..
      SAVE /TG04B/
C     ..
C     .. Data statements ..
      DATA LP/6/
C     ..
C     .. Executable Statements ..
      END
