#include "odepack_cpp/odepack.hpp"

namespace odepack_cpp {
/**
 * @fn DGEFA
 * 
C***BEGIN PROLOGUE  DGEFA
C***PURPOSE  Factor a matrix using Gaussian elimination.
C***CATEGORY  D2A1
C***TYPE      odepack_cpp_real PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
C***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
C             MATRIX FACTORIZATION
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGEFA factors a odepack_cpp_real precision matrix by Gaussian elimination.
C
C     DGEFA is usually called by DGECO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
C
C     On Entry
C
C        A       odepack_cpp_real PRECISION(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGESL or DGEDI will divide by zero
C                     if called.  Use  RCOND  in DGECO for a reliable
C                     indication of singularity.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGEFA
 */
void DGEFA(odepack_cpp_real *a, int lda, int n, int *ipvt, int &info) {
#ifndef MATA
#define MATA(i, j) ARRAY2D(a, lda, i, j)
#endif
//
    odepack_cpp_real t;
    int j, k, kp1, l, nm1;
//
// GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
//
// FIRST EXECUTABLE STATEMENT DGEFA
    info = 0;
    nm1 = n - 1;
    if (nm1 < 1) goto LABEL_70;
    for (k = 1; k <= nm1; ++k) {
        kp1 = k + 1;
//
//      FIND L = PIVOT INDEX
//
        l = IDAMAX(n-k+1, &MATA(k, k), 1) + k - 1;
        ARRAY1D(ipvt, k) = l;
//
//      ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
// 
        if (MATA(l, k) == 0.0) goto LABEL_40;
//
//      intERCHANGE IF NECESSARY
// 
        if (l == k) goto LABEL_10;
        t = MATA(l, k);
        MATA(l, k) = MATA(k, k);
        MATA(k, k) = t;
LABEL_10:
//
//      COMPUTE MULTIPLIERS
//
        t = - 1.0 / MATA(k, k);
        DSCAL(n-k, t, &MATA(k+1, k), 1);
//
//      ROW ELIMINATION WITH COLUMN INDEXING
//
        for (j = kp1; j <= n; ++j) {
            t = MATA(l, j);
            if (l == k) goto LABEL_20;
            MATA(l, j) = MATA(k, j);
            MATA(k, j) = t;
LABEL_20:
            DAXPY(n-k, t, &MATA(k+1, k), 1, &MATA(k+1, j), 1);
        }
        goto LABEL_50;
LABEL_40:
        info = k;
LABEL_50:
        continue;
    }
LABEL_70:
    ARRAY1D(ipvt, n) = n;
    if (MATA(n, n) == 0.0) info = n;
    return;
//
#ifdef MATA
#undef MATA
#endif
}


/**
 * @fn DGESL
 * 
C***BEGIN PROLOGUE  DGESL
C***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the
C            factors computed by DGECO or DGEFA.
C***CATEGORY  D2A1
C***TYPE      odepack_cpp_real PRECISION (SGESL-S, DGESL-D, CGESL-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGESL solves the odepack_cpp_real precision system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGECO or DGEFA.
C
C     On Entry
C
C        A       odepack_cpp_real PRECISION(LDA, N)
C                the output from DGECO or DGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGECO or DGEFA.
C
C        B       odepack_cpp_real PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B  where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGECO has set RCOND .GT. 0.0
C        or DGEFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DDOT
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGESL
*/
void DGESL(odepack_cpp_real *a, const int lda, int n, int *ipvt, odepack_cpp_real *b, int job) {
#ifndef MATA
#define MATA(i, j) ARRAY2D(a, lda, i, j)
#endif
//
    int i, j, k, kb, l, nm1;
    odepack_cpp_real t;
//***FIRST EXECUTABLE STATEMENT  DGESL
    nm1 = n - 1;
    if (job != 0) goto LABEL_50;
//
//  JOB = 0 , SOLVE  A * X = B
//  FIRST SOLVE  L*Y = B
//
    if (nm1 < 1) goto LABEL_30;
    for (k = 1; k <= nm1; ++k) {
        l = ARRAY1D(ipvt, k);
        t = ARRAY1D(b, l);
        if (l == k) goto LABEL_10;
        ARRAY1D(b, l) = ARRAY1D(b, k);
        ARRAY1D(b, k) = t;
LABEL_10:
        DAXPY(n-k, t, &MATA(k+1, k), 1, &ARRAY1D(b, k+1), 1);
    }
LABEL_30:
//
//  NOW SOLVE TRANS(L)*X = Y
//
    for (kb = 1; kb <= n; ++kb) {
        k = n + 1 - kb;
        ARRAY1D(b, k) /= MATA(k, k);
        t = - ARRAY1D(b, k);
        DAXPY(k-1, t, &MATA(1, k), 1, &ARRAY1D(b, 1), 1);
    }
    goto LABEL_100;
LABEL_50:
//
//  JOB = NONZERO, SOLVE  TRANS(A) * X = B
//  FIRST SOLVE  TRANS(U)*Y = B
//
    for (k = 1; k <= n; ++k) {
        t = DDOT(k-1, &MATA(1, k), 1, &ARRAY1D(b, 1), 1);
        ARRAY1D(b, k) = (ARRAY1D(b, k) - t) / MATA(k, k);
    }
//
//  NOW SOLVE  U*X = Y
//
    if (nm1 < 1) goto LABEL_90;
    for (kb = 1; kb <= nm1; ++kb) {
        k  = n - kb;
        ARRAY1D(b, k) += DDOT(n-k, &MATA(k+1, k), 1, &ARRAY1D(b, k+1), 1);
        l = ARRAY1D(ipvt, k);
        if (l == k) goto LABEL_70;
        t = ARRAY1D(b, l);
        ARRAY1D(b, l) = ARRAY1D(b, k);
        ARRAY1D(b, k) = t;
LABEL_70:
        continue;
    }
LABEL_90:
LABEL_100:
    return;
//
#ifdef MATA
#undef MATA
#endif
}


/**
 * @fn DGBFA
 * 
C***BEGIN PROLOGUE  DGBFA
C***PURPOSE  Factor a band matrix using Gaussian elimination.
C***CATEGORY  D2A2
C***TYPE      odepack_cpp_real PRECISION (SGBFA-S, DGBFA-D, CGBFA-C)
C***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGBFA factors a odepack_cpp_real precision band matrix by elimination.
C
C     DGBFA is usually called by DGBCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C
C     On Entry
C
C        ABD     odepack_cpp_real PRECISION(LDA, N)
C                contains the matrix in band storage.  The columns
C                of the matrix are stored in the columns of  ABD  and
C                the diagonals of the matrix are stored in rows
C                ML+1 through 2*ML+MU+1 of  ABD .
C                See the comments below for details.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C                LDA must be .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C                0 .LE. ML .LT.  N .
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C                0 .LE. MU .LT.  N .
C                More efficient if  ML .LE. MU .
C     On Return
C
C        ABD     an upper triangular matrix in band storage and
C                the multipliers which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGBSL will divide by zero if
C                     called.  Use  RCOND  in DGBCO for a reliable
C                     indication of singularity.
C
C     Band Storage
C
C           If  A  is a band matrix, the following program segment
C           will set up the input.
C
C                   ML = (band width below the diagonal)
C                   MU = (band width above the diagonal)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX(1, J-MU)
C                      I2 = MIN(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
C           In addition, the first  ML  rows in  ABD  are used for
C           elements generated during the triangularization.
C           The total number of rows needed in  ABD  is  2*ML+MU+1 .
C           The  ML+MU by ML+MU  upper left triangle and the
C           ML by ML  lower right triangle are not referenced.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGBFA
 */
void DGBFA(odepack_cpp_real *abd, int lda, int n, int ml, int mu, int *ipvt, int &info) {
#ifndef ABD
#define ABD(i, j) ARRAY2D(abd, lda, i, j)
#endif
//
    odepack_cpp_real t;
    int i, i0, j, ju, jz, j0, j1, k, kp1, l, lm, m, mm, nm1;
//
//***FIRST EXECUTABLE STATEMENT  DGBFA
    m = ml + mu + 1;
    info = 0;
//
//  ZERO INITIAL FILL-IN COLUMNS
//
    j0 = mu + 2;
    j1 = std::min(n, m) - 1;
    if (j1 < j0) goto LABEL_30;
    for (jz = j0; jz <= j1; ++jz) {
        i0 = m + 1 - jz;
        for (i = i0; i <= ml; ++i) {
            ABD(i, jz) = 0.0;
        }
    }
LABEL_30:
    jz = j1;
    ju = 0;
//
//  GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
//
    nm1 = n - 1;
    if (nm1 < 1) goto LABEL_130;
    for (k = 1; k <= nm1; ++k) {
        kp1 = k + 1;
//
//      ZERO NEXT FILL-IN COLUMN
//
        jz++;
        if (jz > n) goto LABEL_50;
        if (ml < 1) goto LABEL_50;
        for (i = 1; i <= ml; ++i) {
            ABD(i, jz) = 0.0;
        }
LABEL_50:
//
//      FIND L = PIVOT INDEX
//
        lm = std::min(ml, n-k);
        l = IDAMAX(lm+1, &ABD(m, k), 1) + m - 1;
        ARRAY1D(ipvt, k) = l + k - m;
//
//      ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
//
        if (ABD(l, k) == 0.0) goto LABEL_100;
//
//      intERCHANGE IF NECESSARY
//
        if (l == m) goto LABEL_60;
        t = ABD(l, k);
        ABD(l, k) = ABD(m, k);
        ABD(m, k) = t;
LABEL_60:
//
//      COMPUTE MULTIPLIERS
//
        t = -1.0 / ABD(m, k);
        DSCAL(lm, t, &ABD(m+1, k), 1);
//
//      ROW ELIMINATION WITH COLUMN INDEXING
//
        ju = std::min(std::max(ju, mu + ARRAY1D(ipvt, k)), n);
        mm = m;
        if (ju < kp1) goto LABEL_90;
        for (j = kp1; j <= ju; ++j) {
            l--;
            mm--;
            t = ABD(l, j);
            if (l == mm) goto LABEL_70;
            ABD(l, j) = ABD(mm, j);
            ABD(mm, j) = t;
LABEL_70:
            DAXPY(lm, t, &ABD(m+1, k), 1, &ABD(mm+1, k), 1);
        }
LABEL_90:
        goto LABEL_110;
LABEL_100:
        info = k;
LABEL_110:
LABEL_120:
        continue;
    }
LABEL_130:
    ARRAY1D(ipvt, n) = n;
    if (ABD(m, n) == 0.0) info = n;
    return;
//
#ifdef ABD
#undef ABD
#endif
}


/**
 * @fn DGBSL
 * 
C***BEGIN PROLOGUE  DGBSL
C***PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using
C            the factors computed by DGBCO or DGBFA.
C***CATEGORY  D2A2
C***TYPE      odepack_cpp_real PRECISION (SGBSL-S, DGBSL-D, CGBSL-C)
C***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGBSL solves the odepack_cpp_real precision band system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGBCO or DGBFA.
C
C     On Entry
C
C        ABD     odepack_cpp_real PRECISION(LDA, N)
C                the output from DGBCO or DGBFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGBCO or DGBFA.
C
C        B       odepack_cpp_real PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B , where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGBCO has set RCOND .GT. 0.0
C        or DGBFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DDOT
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGBSL
 */
void DGBSL(odepack_cpp_real *abd, int lda, int n, int ml, int mu, int *ipvt, odepack_cpp_real *b, int job) {
#ifndef ABD
#define ABD(i, j) (ARRAY2D(abd, lda, i, j))
#endif
//
    int m, nm1, k, lm, l, kb, la, lb;
    odepack_cpp_real t;
//***FIRST EXECUTABLE STATEMENT  DGBSL
    m   = mu + ml + 1;
    nm1 = n - 1;
    if (job != 0) goto LABEL_50;
//
//  JOB = 0 , SOLVE  A * X = B
//  FIRST SOLVE L*Y = B
//
    if (ml == 0) goto LABEL_30;
    if (nm1 < 1) goto LABEL_30;
    for (k = 1; k <= nm1; ++k) {
        lm = std::min(ml, n-k);
        l = ARRAY1D(ipvt, k);
        t = ARRAY1D(b, l);
        if (l == k) goto LABEL_10;
        ARRAY1D(b, l) = ARRAY1D(b, k);
        ARRAY1D(b, k) = t;
LABEL_10:
        DAXPY(lm, t, &ABD(m+1, k), 1, &ARRAY1D(b, k+1), 1);
    }
LABEL_30:
//
//  NOW SOLVE  U*X = Y
//
    for (kb = 1; kb <= n; ++kb) {
        k = n + 1 - kb;
        ARRAY1D(b, k) /= ABD(m, k);
        lm = std::min(k, m) - 1;
        la = m - lm;
        lb = k - lm;
        t = - ARRAY1D(b, k);
        DAXPY(lm, t, &ABD(la, k), 1, &ARRAY1D(b, lb), 1);
    }
    goto LABEL_100;
LABEL_50:
//
//  JOB = NONZERO, SOLVE  TRANS(A) * X = B
//  FIRST SOLVE  TRANS(U)*Y = B
//
    for (k = 1; k <= n; ++k) {
        lm = std::min(k, m) - 1;
        la = m - lm;
        lb = k - lm;
        t = DDOT(lm, &ABD(la, k), 1, &ARRAY1D(b, lb), 1);
        ARRAY1D(b, k) = (ARRAY1D(b, k) - t) / ABD(m, k);
    }
//
//  NOW SOLVE  U*X = Y
//
    if (ml == 0) goto LABEL_90;
    if (nm1 < 1) goto LABEL_90;
    for (kb = 1; kb <= nm1; ++kb) {
        k = n - kb;
        lm = std::min(ml, n-k);
        ARRAY1D(b, k) += DDOT(lm, &ABD(m+1, k), 1, &ARRAY1D(b, k+1), 1);
        l = ARRAY1D(ipvt, k);
        if (l == k) goto LABEL_70;
        t = ARRAY1D(b, l);
        ARRAY1D(b, l) = ARRAY1D(   b, k);
        ARRAY1D(b, k) = t;
LABEL_70:
        continue;
    }
LABEL_90:
LABEL_100:
    return;
//
#ifdef ABD
#undef ABD
#endif
}


/**
 * @fn DAXPY
 * 
C***BEGIN PROLOGUE  DAXPY
C***PURPOSE  Compute a constant times a vector plus a vector.
C***CATEGORY  D1A7
C***TYPE      odepack_cpp_real PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  odepack_cpp_real precision scalar multiplier
C       DX  odepack_cpp_real precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  odepack_cpp_real precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  odepack_cpp_real precision result (unchanged if N .LE. 0)
C
C     Overwrite odepack_cpp_real precision DY with odepack_cpp_real precision DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DAXPY
 */
void DAXPY(int n, odepack_cpp_real da, odepack_cpp_real *dx, int incx, odepack_cpp_real *dy, int incy) {
    int i, ix, iy, m, ns, mp1;
//***FIRST EXECUTABLE STATEMENT  DAXPY
    if (n <= 0 || da == 0.0) return;
    if (incx == incy) {
        if (incx < 1) {
            goto LABEL_5;
        } else if (incx == 1) {
            goto LABEL_20;
        } else if (incx > 1) {
            goto LABEL_60;
        }
    }
//
//  Code for unequal or nonpositive increments.
//
LABEL_5:
    ix = 1;
    iy = 1;
    if (incx < 0) ix = (-n + 1) * incx + 1;
    if (incy < 0) iy = (-n + 1) * incy + 1;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(dy, iy) += da * ARRAY1D(dx, ix);
        ix += incx;
        iy += incy;
    }
    return;
//
//  Code for both increments equal to 1.
//
//  Clean-up loop so remaining vector length is a multiple of 4.
//
LABEL_20:
    m = n % 4;
    if (m == 0) goto LABEL_40;
    for (i = 1; i <= m; ++i) {
        ARRAY1D(dy, i) += da * ARRAY1D(dx, i);
    }
    if (n < 4) return;
LABEL_40:
    mp1 = m + 1;
    for (i = mp1; i <= n; i += 4) {
        ARRAY1D(dy,   i) += da * ARRAY1D(dx,   i);
        ARRAY1D(dy, i+1) += da * ARRAY1D(dx, i+1);
        ARRAY1D(dy, i+2) += da * ARRAY1D(dx, i+2);
        ARRAY1D(dy, i+3) += da * ARRAY1D(dx, i+3);
    }
    return;
//
//  Code for equal, positive, non-unit increments.
//
LABEL_60:
    ns = n * incx;
    for (i = 1; i <= ns; i += incx) {
        ARRAY1D(dy, i) += da * ARRAY1D(dx, i);
    }
    return;
}


/**
 * @fn DCOPY
 * 
C***BEGIN PROLOGUE  DCOPY
C***PURPOSE  Copy a vector.
C***CATEGORY  D1A5
C***TYPE      odepack_cpp_real PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
C***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  odepack_cpp_real precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  odepack_cpp_real precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  copy of vector DX (unchanged if N .LE. 0)
C
C     Copy odepack_cpp_real precision DX to odepack_cpp_real precision DY.
C     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DCOPY
 */
void DCOPY(int n, odepack_cpp_real *dx, int incx, odepack_cpp_real *dy, int incy) {
    int i, ix, iy, m, ns, mp1;
//***FIRST EXECUTABLE STATEMENT  DCOPY
    if (n <= 0) return;
    if (incx == incy) {
        if (incx < 1) {
            goto LABEL_5;
        } else if (incx == 1) {
            goto LABEL_20;
        } else {
            goto LABEL_60;
        }
    }
//
//  Code for unequal or nonpositive increments.
//
LABEL_5:
    ix = 1;
    iy = 1;
    if (incx < 0) ix = (-n + 1) * incx + 1;
    if (incy < 0) iy = (-n + 1) * incy + 1;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(dy, iy) = ARRAY1D(dx, ix);
        ix += incx;
        iy += incy;
    }
    return;
//
//  Code for both increments equal to 1.
//
//  Clean-up loop so remaining vector length is a multiple of 7.
//
LABEL_20:
    m = n % 7;
    if (m == 0) goto LABEL_40;
    for (i = 1; i <= m; ++i) {
        ARRAY1D(dy, i) = ARRAY1D(dx, i);
    }
    if (m < 7) return;
LABEL_40:
    mp1 = m + 1;
    for (i = mp1; i <= n; i += 7) {
        ARRAY1D(dy, i  ) = ARRAY1D(dx, i  );
        ARRAY1D(dy, i+1) = ARRAY1D(dx, i+1);
        ARRAY1D(dy, i+2) = ARRAY1D(dx, i+2);
        ARRAY1D(dy, i+3) = ARRAY1D(dx, i+3);
        ARRAY1D(dy, i+4) = ARRAY1D(dx, i+4);
        ARRAY1D(dy, i+5) = ARRAY1D(dx, i+5);
        ARRAY1D(dy, i+6) = ARRAY1D(dx, i+6);
    }
    return;
//
// Code for equal, positive, non-unit increments.
//
LABEL_60:
    ns = n * incx;
    for (i = 1; i <= ns; i += incx) {
        ARRAY1D(dy, i) = ARRAY1D(dx, i);
    }
    return;
}


/**
 * @fn DDOT
 * 
C***BEGIN PROLOGUE  DDOT
C***PURPOSE  Compute the inner product of two vectors.
C***CATEGORY  D1A4
C***TYPE      odepack_cpp_real PRECISION (SDOT-S, DDOT-D, CDOTU-C)
C***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  odepack_cpp_real precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  odepack_cpp_real precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C     DDOT  odepack_cpp_real precision dot product (zero if N .LE. 0)
C
C     Returns the dot product of odepack_cpp_real precision DX and DY.
C     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DDOT
 */
odepack_cpp_real DDOT(int n, odepack_cpp_real *dx, int incx, odepack_cpp_real *dy, int incy) {
    int i, ix, iy, m, ns, mp1;
    odepack_cpp_real ddot;
//***FIRST EXECUTABLE STATEMENT DDOT
    ddot = 0.0;
    if (n <= 0) return ddot;
    if (incx == incy) {
        if (incx < 1) {
            goto LABEL_5;
        } else if (incx == 1) {
            goto LABEL_20;
        } else {
            goto LABEL_60;
        }
    }
//
//  Code for unequal or nonpositive increments.
//
LABEL_5:
    ix = 1;
    iy = 1;
    if (incx < 0) ix = (-n + 1) * incx + 1;
    if (incy < 0) iy = (-n + 1) * incy + 1;
    for (i = 0; i < n; ++i) {
        ddot += ARRAY1D(dx, ix) * ARRAY1D(dy, iy);
        ix += incx;
        iy += incy;
    }
    return ddot;
//
//  Code for both increments equal to 1.
//
//  Clean-up loop so remaining vector length is a multiple of 5.
//
LABEL_20:
    m = n % 5;
    if (m == 0) goto LABEL_40;
    for (i = 1; i <= m; ++i) {
        ddot += ARRAY1D(dx, i) * ARRAY1D(dy, i);
    }
    if (n < 5) return ddot;
LABEL_40:
    mp1 = m + 1;
    for (i = mp1; i <= n; i += 5) {
        ddot += ARRAY1D(dx,   i) * ARRAY1D(dy,   i)
            + ARRAY1D(dx, i+1) * ARRAY1D(dy, i+1)
            + ARRAY1D(dx, i+2) * ARRAY1D(dy, i+2)
            + ARRAY1D(dx, i+3) * ARRAY1D(dy, i+3)
            + ARRAY1D(dx, i+4) * ARRAY1D(dy, i+4);
    }
    return ddot;
//
// Code for equal, positive, non-unit increments.
//
LABEL_60:
    ns = n * incx;
    for (i = 1; i <= ns; i += incx) {
        ddot += ARRAY1D(dx, i) * ARRAY1D(dy, i);
    }
    return ddot;
}


/**
 * @fn DNRM2
 * 
C***BEGIN PROLOGUE  DNRM2
C***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.
C***CATEGORY  D1A3B
C***TYPE      odepack_cpp_real PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)
C***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
C             LINEAR ALGEBRA, UNITARY, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  odepack_cpp_real precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C    DNRM2  odepack_cpp_real precision result (zero if N .LE. 0)
C
C     Euclidean norm of the N-vector stored in DX with storage
C     increment INCX.
C     If N .LE. 0, return with result = 0.
C     If N .GE. 1, then INCX must be .GE. 1
C
C     Four phase method using two built-in constants that are
C     hopefully applicable to all machines.
C         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
C         CUTHI = minimum of  SQRT(V)      over all known machines.
C     where
C         EPS = smallest no. such that EPS + 1. .GT. 1.
C         U   = smallest positive no.   (underflow limit)
C         V   = largest  no.            (overflow  limit)
C
C     Brief outline of algorithm.
C
C     Phase 1 scans zero components.
C     move to phase 2 when a component is nonzero and .LE. CUTLO
C     move to phase 3 when a component is .GT. CUTLO
C     move to phase 4 when a component is .GE. CUTHI/M
C     where M = N for X() real and M = 2*N for complex.
C
C     Values for CUTLO and CUTHI.
C     From the environmental parameters listed in the IMSL converter
C     document the limiting values are as follows:
C     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
C                   Univac and DEC at 2**(-103)
C                   Thus CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
C                   Thus CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
C                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
C     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DNRM2
 */
odepack_cpp_real DNRM2(int n, odepack_cpp_real *dx, int incx) {
    odepack_cpp_real dnrm2, sum, xmax, hitest;
    int i, j, nn, next;
    odepack_cpp_real zero  = 0.0;
    odepack_cpp_real one   = 1.0;
//
    odepack_cpp_real cutlo = 8.232e-11;
    odepack_cpp_real cuthi = 1.304e19;
//***FIRST EXECUTABLE STATEMENT  DNRM2
    if (n > 0) goto LABEL_10;
    dnrm2 = zero;
    goto LABEL_300;
//
LABEL_10:
    next = 30;
    sum  = zero;
    nn   = n * incx;
//
// BEGIN MAIN LOOP
//
    i = 1;
LABEL_20:
    if (next == 30) {
        goto LABEL_30;
    } else if (next == 50) {
        goto LABEL_50;
    } else if (next == 70) {
        goto LABEL_70;
    } else if (next == 110) {
        goto LABEL_110;
    }
LABEL_30:
    if (std::abs(ARRAY1D(dx, i)) > cutlo) goto LABEL_85;
    next = 50;
    xmax = zero;
//
//  PHASE 1.  SUM IS ZERO
//
LABEL_50:
    if (ARRAY1D(dx, i) == zero) goto LABEL_200;
    if (std::abs(ARRAY1D(dx, i)) > cutlo) goto LABEL_85;
//
//  PREPARE FOR PHASE 2.
//
    next = 70;
    goto LABEL_105;
//
//  PREPARE FOR PHASE 4.
//
LABEL_100:
    i = j;
    next = 110;
    sum = (sum / ARRAY1D(dx, i)) / ARRAY1D(dx, i);
LABEL_105:
    xmax = std::abs(ARRAY1D(dx, i));
    goto LABEL_115;
//
//  PHASE 2.  SUM IS SMALL.
//  SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
//
LABEL_70:
    if (std::abs(ARRAY1D(dx, i)) > cutlo) goto LABEL_75;
//
//  COMMON CODE FOR PHASES 2 AND 4.
//  IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
//
LABEL_110:
    if (std::abs(ARRAY1D(dx, i)) <= xmax) goto LABEL_115;
    sum  = one + sum * (xmax / ARRAY1D(dx, i)) * (xmax / ARRAY1D(dx, i));
    xmax = std::abs(ARRAY1D(dx, i));
    goto LABEL_200;
//
LABEL_115:
    sum = sum + (ARRAY1D(dx, i) / xmax) * (ARRAY1D(dx, i) / xmax);
    goto LABEL_200;
//
//  PREPARE FOR PHASE 3.
//
LABEL_75:
    sum = (sum * xmax) * xmax;
//
//  FOR odepack_cpp_real OR D.P. SET HITEST = CUTHI/N
//  FOR COMPLEX      SET HITEST = CUTHI/(2*N)
//
LABEL_85:
    hitest = cuthi / static_cast<odepack_cpp_real>(n);
//
//  PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
//
    for (j = i; i <= nn; j += incx) {
        if (std::abs(ARRAY1D(dx, j)) >= hitest) goto LABEL_100;
        sum = sum + ARRAY1D(dx, j) * ARRAY1D(dx, j);
    }
    dnrm2 = std::sqrt(sum);
    goto LABEL_300;
//
LABEL_200:
    i = i + incx;
    if (i <= nn) goto LABEL_20;
//
//  END OF MAIN LOOP.
//
//  COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
//
    dnrm2 = xmax * std::sqrt(sum);
LABEL_300:
    return dnrm2;
}


/**
 * @fn DSCAL
 * 
C***BEGIN PROLOGUE  DSCAL
C***PURPOSE  Multiply a vector by a constant.
C***CATEGORY  D1A6
C***TYPE      odepack_cpp_real PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  odepack_cpp_real precision scale factor
C       DX  odepack_cpp_real precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C       DX  odepack_cpp_real precision result (unchanged if N.LE.0)
C
C     Replace odepack_cpp_real precision DX by odepack_cpp_real precision DA*DX.
C     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DSCAL
 */
void DSCAL(int n, odepack_cpp_real da, odepack_cpp_real *dx, int incx) {
    int i, ix, m, mp1;
//***FIRST EXECUTABLE STATEMENT  DSCAL
    if (n <= 0) return;
    if (incx == 1) goto LABEL_20;
//
// Code for increment not equal to 1.
//
    ix = 1;
    if (incx < 0) ix = (-n + 1) * incx + 1;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(dx, ix) *= da;
        ix += incx;
    }
    return;
//
// Code for increment equal to 1.
// 
// Clean-up loop so remaining vector length is a multiple of 5.
//
LABEL_20:
    m = n % 5;
    if (m == 0) goto LABEL_40;
    for (i = 1; i <= m; ++i) {
        ARRAY1D(dx, i) *= da;
    }
    if (n < 5) return;
LABEL_40:
    mp1 = m + 1;
    for (i = mp1; i <= n; i += 5) {
        ARRAY1D(dx,   i) *= da;
        ARRAY1D(dx, i+1) *= da;
        ARRAY1D(dx, i+2) *= da;
        ARRAY1D(dx, i+3) *= da;
        ARRAY1D(dx, i+4) *= da;
    }
    return;
}


/**
 * @fn IDAMAX
 * 
C***BEGIN PROLOGUE  IDAMAX
C***PURPOSE  Find the smallest index of that component of a vector
C            having the maximum magnitude.
C***CATEGORY  D1A2
C***TYPE      odepack_cpp_real PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  odepack_cpp_real precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C   IDAMAX  smallest index (zero if N .LE. 0)
C
C     Find smallest index of maximum magnitude of odepack_cpp_real precision DX.
C     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  IDAMAX
 */
int IDAMAX(int n, odepack_cpp_real *dx, int incx) {
    int i, ix, idamax;
    odepack_cpp_real dmax, xmag;
//***FIRST EXECUTABLE STATEMENT  IDAMAX
    idamax = 0;
    if (n <= 0) return idamax;
    idamax = 1;
    if (n == 1) return idamax;
//
    if (incx == 1) goto LABEL_20;
//
//  Code for increments not equal to 1.
//
    ix = 1;
    if (incx < 0) ix = (-n + 1) * incx + 1;
    dmax = std::abs(ARRAY1D(dx, ix));
    ix += incx;
    for (i = 2; i <= n; ++i) {
        xmag = std::abs(ARRAY1D(dx, ix));
        if (xmag > dmax) {
            idamax = i;
            dmax = xmag;
        }
        ix += incx;
    }
    return idamax;
//
//  Code for increments equal to 1.
//
LABEL_20:
    dmax = std::abs(ARRAY1D(dx, 1));
    for (i = 2; i <= n; ++i) {
        xmag = std::abs(ARRAY1D(dx, i));
        if (xmag > dmax) {
            idamax = i;
            dmax = xmag;
        }
    }
    return idamax;
}


/**
 * @fn XERRWD
 * 
C***BEGIN PROLOGUE  XERRWD
C***SUBSIDIARY
C***PURPOSE  Write error message with values.
C***CATEGORY  R3C
C***TYPE      odepack_cpp_real PRECISION (XERRWV-S, XERRWD-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  Subroutines XERRWD, XSETF, XSETUN, and the function routine IXSAV,
C  as given here, constitute a simplified version of the SLATEC error
C  handling package.
C
C  All arguments are input arguments.
C
C  MSG    = The message (character array).
C  NMES   = The length of MSG (number of characters).
C  NERR   = The error number (not used).
C  LEVEL  = The error level..
C           0 or 1 means recoverable (control returns to caller).
C           2 means fatal (run is aborted--see note below).
C  NI     = Number of integers (0, 1, or 2) to be printed with message.
C  I1,I2  = Integers to be printed, depending on NI.
C  NR     = Number of reals (0, 1, or 2) to be printed with message.
C  R1,R2  = Reals to be printed, depending on NR.
C
C  Note..  this routine is machine-dependent and specialized for use
C  in limited context, in the following ways..
C  1. The argument MSG is assumed to be of type CHARACTER, and
C     the message is printed with a format of (1X,A).
C  2. The message is assumed to take only one line.
C     Multi-line messages are generated by repeated calls.
C  3. If LEVEL = 2, control passes to the statement   STOP
C     to abort the run.  This statement may be machine-dependent.
C  4. R1 and R2 are assumed to be in odepack_cpp_real precision and are printed
C     in D21.13 format.
C
C***ROUTINES CALLED  IXSAV
C***REVISION HISTORY  (YYMMDD)
C   920831  DATE WRITTEN
C   921118  Replaced MFLGSV/LUNSAV by IXSAV. (ACH)
C   930329  Modified prologue to SLATEC format. (FNF)
C   930407  Changed MSG from CHARACTER*1 array to variable. (FNF)
C   930922  Minor cosmetic change. (FNF)
C***END PROLOGUE  XERRWD
C
C*Internal Notes:
C
C For a different default logical unit number, IXSAV (or a subsidiary
C routine that it calls) will need to be modified.
C For a different run-abort command, change the statement following
C statement 100 at the end.
C-----------------------------------------------------------------------
C Subroutines called by XERRWD.. None
C Function routine called by XERRWD.. IXSAV
C-----------------------------------------------------------------------
C**End
 */
void XERRWD(std::string msg, int nmes, int nerr, int level, int ni, int i1, int i2, int nr, odepack_cpp_real r1, odepack_cpp_real r2) {
    int lunit, mesflg;
//***FIRST EXECUTABLE STATEMENT  XERRWD
    lunit = IXSAV(1, 0, false);
    mesflg = IXSAV(2, 0, false);
    if (mesflg == 0) goto LABEL_100;
//
//  Write the message.
//
    std::cout << msg << std::endl;
    if (ni == 1) {
        std::cout << "In above message, i1 = " << i1 << std::endl;
    }
    if (ni == 2) {
        std::cout << "In above message, i1 = " << i1 << ", i2 = " << i2 << std::endl;
    }
    if (nr == 1) {
        std::cout << "In above message, r1 = " << r1 << std::endl;
    }
    if (nr == 2) {
        std::cout << "In above message, r1 = " << r1 << ", r2 = " << r2 << std::endl;
    }
//
//  Abort the run if LEVEL = 2.
//
LABEL_100:
    if (level != 2) return;
    std::exit(1);
}


int IXSAV(int ipar, int ivalue, bool iset) {
    int ixsav;
    int lunit = -1;
    int mesflg = 1;
//***FIRST EXECUTABLE STATEMENT  IXSAV
    if (ipar == 1) {
        if (lunit == -1) lunit = IUMACH();
        ixsav = lunit;
        if (iset) lunit = ivalue;
    }
//
    if (ipar == 2) {
        ixsav = mesflg;
        if (iset) mesflg = ivalue;
    }
// 
    return ixsav;
}


int IUMACH() {
    return 6;
}
} // end namespace odepack_cpp
