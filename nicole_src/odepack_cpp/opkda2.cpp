/**
 * @fn opkda2.cpp
 * @note
 * 2025/01/21 kawasaki
*/

#include "odepack.hpp"

namespace odepack_cpp
{

void Odepack::DGEFA(double *a, const int lda, const int n, int *ipvt, int &info)
{
#ifndef MATA
#define MATA(i, j) MATF(a, lda, i, j)
#endif
//
    double t;
    int j, k, kp1, l, nm1;
//
// GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
//
// FIRST EXECUTABLE STATEMENT DGEFA
    info = 0;
    nm1  = n - 1;
    if (nm1 < 1) goto LABEL_70;
    for (k = 1; k <= nm1; ++k) {
        kp1 = k + 1;
//
//      FIND L = PIVOT INDEX
//
        l = IDAMAX(n-k+1, &MATA(k, k), 1) + k - 1;
        ARRAYF(ipvt, k) = l;
//
//      ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
// 
        if (MATA(l, k) == 0.0) goto LABEL_40;
//
//      INTERCHANGE IF NECESSARY
// 
        if (l == k) goto LABEL_10;
        t          = MATA(l, k);
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
            t          = MATA(l, j);
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
    ARRAYF(ipvt, n) = n;
    if (MATA(n, n) == 0.0) info = n;
    return;
//
#ifdef MATA
#undef MATA
#endif
}

/*
    A : n*n array
    n : size
    ipvt :
    b : n array
    job
*/
void Odepack::DGESL(const double *a, const int lda, const int n, int *ipvt, double *b, const int job)
{
#ifndef MATA
#define MATA(i, j) MATF(a, lda, i, j)
#endif
//
    int i, j, k, kb, l, nm1;
    double t;
//***FIRST EXECUTABLE STATEMENT  DGESL
    nm1 = n - 1;
    if (job != 0) goto LABEL_50;
//
//  JOB = 0 , SOLVE  A * X = B
//  FIRST SOLVE  L*Y = B
//
    if (nm1 < 1) goto LABEL_30;
    for (k = 1; k <= nm1; ++k) {
        l            = ARRAYF(ipvt, k);
        t            = ARRAYF(   b, l);
        if (l == k) goto LABEL_10;
        ARRAYF(b, l) = ARRAYF(b, k);
        ARRAYF(b, k) = t;
LABEL_10:
        DAXPY(n-k, t, &MATA(k+1, k), 1, &ARRAYF(b, k+1), 1);
    }
LABEL_30:
//
//  NOW SOLVE TRANS(L)*X = Y
//
    for (kb = 1; kb <= n; ++kb) {
        k             = n + 1 - kb;
        ARRAYF(b, k) /= MATA(k, k);
        t             = - ARRAYF(b, k);
        DAXPY(k-1, t, &MATA(1, k), 1, &ARRAYF(b, 1), 1);
    }
    goto LABEL_100;
LABEL_50:
//
//  JOB = NONZERO, SOLVE  TRANS(A) * X = B
//  FIRST SOLVE  TRANS(U)*Y = B
//
    for (k = 1; k <= n; ++k) {
        t            = DDOT(k-1, &MATA(1, k), 1, &ARRAYF(b, 1), 1);
        ARRAYF(b, k) = (ARRAYF(b, k) - t) / MATA(k, k);
    }
//
//  NOW SOLVE  U*X = Y
//
    if (nm1 < 1) goto LABEL_90;
    for (kb = 1; kb <= nm1; ++kb) {
        k             = n - kb;
        ARRAYF(b, k) += DDOT(n-k, &MATA(k+1, k), 1, &ARRAYF(b, k+1), 1);
        l             = ARRAYF(ipvt, k);
        if (l == k) goto LABEL_70;
        t             = ARRAYF(b, l);
        ARRAYF(b, l)  = ARRAYF(b, k);
        ARRAYF(b, k)  = t;
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


void Odepack::DGBFA(double *abd, const int lda, const int n, const int ml, const int mu, int *ipvt, int &info)
{
#ifndef ABD
#define ABD(i, j) MATF(abd, lda, i, j)
#endif
//
    double t;
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
        jz += 1;
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
        l  = IDAMAX(lm+1, &ABD(m, k), 1) + m - 1;
        ARRAYF(ipvt, k) = l + k - m;
//
//      ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
//
        if (ABD(l, k) == 0.0) goto LABEL_100;
//
//      INTERCHANGE IF NECESSARY
//
        if (l == m) goto LABEL_60;
        t         = ABD(l, k);
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
        ju = std::min(std::max(ju, mu + ARRAYF(ipvt, k)), n);
        mm = m;
        if (ju < kp1) goto LABEL_90;
        for (j = kp1; j <= ju; ++j) {
            l--;
            mm--;
            t = ABD(l, j);
            if (l == mm) goto LABEL_70;
            ABD( l, j) = ABD(mm, j);
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
    ARRAYF(ipvt, n) = n;
    if (ABD(m, n) == 0.0) info = n;
    return;
//
#ifdef ABD
#undef ABD
#endif
}

void Odepack::DGBSL(double *abd, const int lda, const int n, const int ml, const int mu, int *ipvt, double *b, const int job)
{
#ifndef ABD
#define ABD(i, j) (MATF(abd, lda, i, j))
#endif
//
    int m, nm1, k, lm, l, kb, la, lb;
    double t;
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
        l  = ARRAYF(ipvt, k);
        t  = ARRAYF(   b, l);
        if (l == k) goto LABEL_10;
        ARRAYF(b, l) = ARRAYF(b, k);
        ARRAYF(b, k) = t;
LABEL_10:
        DAXPY(lm, t, &ABD(m+1, k), 1, &ARRAYF(b, k+1), 1);
    }
LABEL_30:
//
//  NOW SOLVE  U*X = Y
//
    for (kb = 1; kb <= n; ++kb) {
        k  = n + 1 - kb;
        ARRAYF(b, k) /= ABD(m, k);
        lm = std::min(k, m) - 1;
        la = m - lm;
        lb = k - lm;
        t  = - ARRAYF(b, k);
        DAXPY(lm, t, &ABD(la, k), 1, &ARRAYF(b, lb), 1);
    }
    goto LABEL_100;
LABEL_50:
//
//  JOB = NONZERO, SOLVE  TRANS(A) * X = B
//  FIRST SOLVE  TRANS(U)*Y = B
//
    for (k = 1; k <= n; ++k) {
        lm           = std::min(k, m) - 1;
        la           = m - lm;
        lb           = k - lm;
        t            = DDOT(lm, &ABD(la, k), 1, &ARRAYF(b, lb), 1);
        ARRAYF(b, k) = (ARRAYF(b, k) - t) / ABD(m, k);
    }
//
//  NOW SOLVE  U*X = Y
//
    if (ml == 0) goto LABEL_90;
    if (nm1 < 1) goto LABEL_90;
    for (kb = 1; kb <= nm1; ++kb) {
        k            = n - kb;
        lm           = std::min(ml, n-k);
        ARRAYF(b, k) = ARRAYF(   b, k) + DDOT(lm, &ABD(m+1, k), 1, &ARRAYF(b, k+1), 1);
        l            = ARRAYF(ipvt, k);
        if (l == k) goto LABEL_70;
        t            = ARRAYF(   b, l);
        ARRAYF(b, l) = ARRAYF(   b, k);
        ARRAYF(b, k) = t;
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

void Odepack::DAXPY(const int n, const double da, const double *dx, const int incx, double *dy, const int incy)
{
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
        ARRAYF(dy, iy) += da * ARRAYF(dx, ix);
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
        ARRAYF(dy, i) += da * ARRAYF(dx, i);
    }
    if (n < 4) return;
LABEL_40:
    mp1 = m + 1;
    for (i = mp1; i <= n; i += 4) {
        ARRAYF(dy,   i) += da * ARRAYF(dx,   i);
        ARRAYF(dy, i+1) += da * ARRAYF(dx, i+1);
        ARRAYF(dy, i+2) += da * ARRAYF(dx, i+2);
        ARRAYF(dy, i+3) += da * ARRAYF(dx, i+3);
    }
    return;
//
//  Code for equal, positive, non-unit increments.
//
LABEL_60:
    ns = n * incx;
    for (i = 1; i <= ns; i += incx) {
        ARRAYF(dy, i) += da * ARRAYF(dx, i);
    }
    return;
}

void Odepack::DCOPY(const int n, const double *dx, const int incx, double *dy, const int incy)
{
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
        ARRAYF(dy, iy) = ARRAYF(dx, ix);
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
        ARRAYF(dy, i) = ARRAYF(dx, i);
    }
    if (m < 7) return;
LABEL_40:
    mp1 = m + 1;
    for (i = mp1; i <= n; i += 7) {
        ARRAYF(dy,   i) = ARRAYF(dx,   i);
        ARRAYF(dy, i+1) = ARRAYF(dx, i+1);
        ARRAYF(dy, i+2) = ARRAYF(dx, i+2);
        ARRAYF(dy, i+3) = ARRAYF(dx, i+3);
        ARRAYF(dy, i+4) = ARRAYF(dx, i+4);
        ARRAYF(dy, i+5) = ARRAYF(dx, i+5);
        ARRAYF(dy, i+6) = ARRAYF(dx, i+6);
    }
    return;
//
// Code for equal, positive, non-unit increments.
//
LABEL_60:
    ns = n * incx;
    for (i = 1; i <= ns; i += incx) {
        ARRAYF(dy, i) = ARRAYF(dx, i);
    }
    return;
}

double Odepack::DDOT(const int n, const double *dx, const int incx, const double *dy, const int incy)
{
    int i, ix, iy, m, ns, mp1;
    double ddot;
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
        ddot += ARRAYF(dx, ix) * ARRAYF(dy, iy);
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
        ddot += ARRAYF(dx, i) * ARRAYF(dy, i);
    }
    if (n < 5) return ddot;
LABEL_40:
    mp1 = m + 1;
    for (i = mp1; i <= n; i += 5) {
        ddot += ARRAYF(dx,   i) * ARRAYF(dy,   i)
              + ARRAYF(dx, i+1) * ARRAYF(dy, i+1)
              + ARRAYF(dx, i+2) * ARRAYF(dy, i+2)
              + ARRAYF(dx, i+3) * ARRAYF(dy, i+3)
              + ARRAYF(dx, i+4) * ARRAYF(dy, i+4);
    }
    return ddot;
//
// Code for equal, positive, non-unit increments.
//
LABEL_60:
    ns = n * incx;
    for (i = 1; i <= ns; i += incx) {
        ddot += ARRAYF(dx, i) * ARRAYF(dy, i);
    }
    return ddot;
}

double Odepack::DNRM2(const int n, const double *dx, const int incx)
{
    double dnrm2, sum, xmax, hitest;
    int i, j, nn, next;
    const double zero  = 0.0;
    const double one   = 1.0;
//
    const double cutlo = 8.232e-11;
    const double cuthi = 1.304e19;
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
    if (std::abs(ARRAYF(dx, i)) > cutlo) goto LABEL_85;
    next = 50;
    xmax = zero;
//
//  PHASE 1.  SUM IS ZERO
//
LABEL_50:
    if (ARRAYF(dx, i) == zero) goto LABEL_200;
    if (std::abs(ARRAYF(dx, i)) > cutlo) goto LABEL_85;
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
    sum = (sum / ARRAYF(dx, i)) / ARRAYF(dx, i);
LABEL_105:
    xmax = std::abs(ARRAYF(dx, i));
    goto LABEL_115;
//
//  PHASE 2.  SUM IS SMALL.
//  SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
//
LABEL_70:
    if (std::abs(ARRAYF(dx, i)) > cutlo) goto LABEL_75;
//
//  COMMON CODE FOR PHASES 2 AND 4.
//  IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
//
LABEL_110:
    if (std::abs(ARRAYF(dx, i)) <= xmax) goto LABEL_115;
    sum  = one + sum * (xmax / ARRAYF(dx, i)) * (xmax / ARRAYF(dx, i));
    xmax = std::abs(ARRAYF(dx, i));
    goto LABEL_200;
//
LABEL_115:
    sum = sum + (ARRAYF(dx, i) / xmax) * (ARRAYF(dx, i) / xmax);
    goto LABEL_200;
//
//  PREPARE FOR PHASE 3.
//
LABEL_75:
    sum = (sum * xmax) * xmax;
//
//  FOR REAL OR D.P. SET HITEST = CUTHI/N
//  FOR COMPLEX      SET HITEST = CUTHI/(2*N)
//
LABEL_85:
    hitest = cuthi / static_cast<double>(n);
//
//  PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
//
    for (j = i; i <= nn; j += incx) {
        if (std::abs(ARRAYF(dx, j)) >= hitest) goto LABEL_100;
        sum = sum + ARRAYF(dx, j) * ARRAYF(dx, j);
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

void Odepack::DSCAL(const int n, const double da, double *dx, const int incx)
{
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
        ARRAYF(dx, ix) *= da;
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
        ARRAYF(dx, i) *= da;
    }
    if (n < 5) return;
LABEL_40:
    mp1 = m + 1;
    for (i = mp1; i <= n; i += 5) {
        ARRAYF(dx,   i) *= da;
        ARRAYF(dx, i+1) *= da;
        ARRAYF(dx, i+2) *= da;
        ARRAYF(dx, i+3) *= da;
        ARRAYF(dx, i+4) *= da;
    }
    return;
}

int Odepack::IDAMAX(const int n, double *dx, const int incx)
{
    int i, ix, idamax;
    double dmax, xmag;
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
    dmax = std::abs(ARRAYF(dx, ix));
    ix += incx;
    for (i = 2; i <= n; ++i) {
        xmag = std::abs(ARRAYF(dx, ix));
        if (xmag > dmax) {
            idamax = i;
            dmax   = xmag;
        }
        ix += incx;
    }
    return idamax;
//
//  Code for increments equal to 1.
//
LABEL_20:
    dmax = std::abs(ARRAYF(dx, 1));
    for (i = 2; i <= n; ++i) {
        xmag = std::abs(ARRAYF(dx, i));
        if (xmag > dmax) {
            idamax = i;
            dmax   = xmag;
        }
    }
    return idamax;
}

void Odepack::XERRWD(const std::string msg, const int nmes, const int nerr, const int level, 
        const int ni, const int i1, const int i2, const int nr, const double r1, const double r2)
{
    int lunit, mesflg;
//***FIRST EXECUTABLE STATEMENT  XERRWD
    lunit  = IXSAV(1, 0, false);
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

int Odepack::IXSAV(const int ipar, const int ivalue, const bool iset)
{
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

int Odepack::IUMACH()
{
    return 6;
}

} // end namespace odepack_cpp