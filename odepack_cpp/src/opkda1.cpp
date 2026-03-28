#include "odepack_cpp/odepack.hpp"

namespace odepack_cpp {
/**
 * PURPOSE  Compute the unit roundoff of the machine.
 * AUTHOR  Hindmarsh, Alan C., (LLNL)
 * DESCRIPTION
 * 
 * The unit roundoff is defined as the smallest positive machine
 * number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
 * in a machine-independent manner.
 * 
 * @return the unit roundoff of the machine.
 */
odepack_cpp_real Odepack::DUMACH() {
    odepack_cpp_real u, comp;
    u = 1.0;
LABEL_10:
    u *= 0.5;
    DUMSUM(1.0, u, comp);
    if (comp != 1.0) goto LABEL_10;
    return u * 2.0;
}


/**
 * Routine to force normal storing of A + B, for DUMACH.
 */
void Odepack::DUMSUM(odepack_cpp_real a, odepack_cpp_real b, odepack_cpp_real &c) {
    c = a + b;
    return;
}


/**
 * PURPOSE  Set ODE integrator coefficients.
 * AUTHOR  Hindmarsh, Alan C., (LLNL)
 * DESCRIPTION
 * 
 * DCFODE is called by the integrator routine to set coefficients
 * needed there. The coefficients for the current method, as
 * given by the value of METH, are set for all orders and saved.
 * The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2.
 * (A smaller value of the maximum order is also allowed.)
 * DCFODE is called once at the beginning of the problem,
 * and is not called again unless and until METH is changed.
 * 
 * The ELCO array contains the basic method coefficients.
 * The coefficients el(i), 1 .le. i .le. nq+1, for the method of
 * order nq are stored in ELCO(i,nq).  They are given by a genetrating
 * polynomial, i.e.,
 *     l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
 * For the implicit Adams methods, l(x) is given by
 *     dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
 * For the BDF methods, l(x) is given by
 *     l(x) = (x+1)*(x+2)* ... *(x+nq)/K,
 * where         K = factorial(nq)*(1 + 1/2 + ... + 1/nq).
 * 
 * The TESCO array contains test constants used for the
 * local error test and the selection of step size and/or order.
 * At order nq, TESCO(k,nq) is used for the selection of step
 * size at order nq - 1 if k = 1, at order nq if k = 2, and at order
 * nq + 1 if k = 3.
 */
void Odepack::DCFODE(int meth, odepack_cpp_real *elco, odepack_cpp_real *tesco) {
#ifndef ELCO
#define ELCO(i, j) ARRAY2D(elco, 13, i, j)
#endif
#ifndef TESCO
#define TESCO(i, j) ARRAY2D(tesco, 3, i, j)
#endif
//
    int i, ib, nq, nqm1, nqp1;
    odepack_cpp_real agamq, fnq, fnqm1, pc[12], pint, ragq, 
        rqfac, rq1fac, tsign, xpin;
//***FIRST EXECUTABLE STATEMENT  DCFODE
    if (meth == 1) {
        goto LABEL_100;
    } else if (meth == 2) {
        goto LABEL_200;
    }
//
LABEL_100:
    ELCO(1, 1) = 1.0;
    ELCO(2, 1) = 1.0;
    TESCO(1, 1) = 0.0;
    TESCO(2, 1) = 2.0;
    TESCO(1, 2) = 1.0;
    TESCO(3, 12) = 0.0;
    ARRAY1D(pc, 1) = 1.0;
    rqfac = 1.0;
    for (nq = 2; nq <= 12; ++nq) {
//-----------------------------------------------------------------------
// The PC array will contain the coefficients of the polynomial
// p(x) = (x+1)*(x+2)*...*(x+nq-1).
// Initially, p(x) = 1.
//-----------------------------------------------------------------------
        rq1fac = rqfac;
        rqfac = rqfac / static_cast<odepack_cpp_real>(nq);
        nqm1 = nq - 1;
        fnqm1 = static_cast<odepack_cpp_real>(nqm1);
        nqp1 = nq + 1;
// Form coefficients of p(x)*(x+nq-1). ----------------------------------
        ARRAY1D(pc, nq) = 0.0;
        for (ib = 1; ib <= nqm1; ++ib) {
            i = nqp1 - ib;
            ARRAY1D(pc, i) = ARRAY1D(pc, i - 1) + fnqm1 * ARRAY1D(pc, i);
        }
        ARRAY1D(pc, 1) = fnqm1 * ARRAY1D(pc, 1);
// Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        pint = ARRAY1D(pc, 1);
        xpin = ARRAY1D(pc, 1) / 2.0;
        tsign = 1.0;
        for (i = 2; i <= nq; ++i) {
            tsign = - tsign;
            pint = pint + tsign * ARRAY1D(pc, i) / static_cast<odepack_cpp_real>(i);
            xpin = xpin + tsign * ARRAY1D(pc, i) / static_cast<odepack_cpp_real>(i + 1);
        }
// Store coefficients in elco and tesco. --------------------------------
        ELCO(1, nq) = pint * rq1fac;
        ELCO(2, nq) = 1.0;
        for (i = 2; i <= nq; ++i) {
            ELCO(i + 1, nq) = rq1fac * ARRAY1D(pc, i) / static_cast<odepack_cpp_real>(i);
        }
        agamq = rqfac * xpin;
        ragq = 1.0 / agamq;
        TESCO(2, nq) = ragq;
        if (nq < 12) TESCO(1, nqp1) = ragq * rqfac / static_cast<odepack_cpp_real>(nqp1);
        TESCO(3, nqm1) = ragq;
    }
    return;
//
LABEL_200:
    ARRAY1D(pc, 1) = 1.0;
    rq1fac = 1.0;
    for (nq = 1; nq <= 5; ++nq) {
//-----------------------------------------------------------------------
// The PC array will contain the coefficients of the polynomial
//     p(x) = (x+1)*(x+2)*...*(x+nq).
// Initially, p(x) = 1.
//-----------------------------------------------------------------------
        fnq = static_cast<odepack_cpp_real>(nq);
        nqp1 = nq + 1;
// Form coefficients of p(x)*(x+nq). ------------------------------------
        ARRAY1D(pc, nqp1) = 0.0;
        for (ib = 1; ib <= nq; ++ib) {
            i = nq + 2 - ib;
            ARRAY1D(pc, i) = ARRAY1D(pc, i - 1) + fnq * ARRAY1D(pc, i);
        }
        ARRAY1D(pc, 1) = fnq * ARRAY1D(pc, 1);
// Store coefficients in elco and tesco. --------------------------------
        for (i = 1; i <= nqp1; ++i) {
            ELCO(i, nq) = ARRAY1D(pc, i) / ARRAY1D(pc, 2);
        }
        ELCO(2, nq) = 1.0;
        TESCO(1, nq) = rq1fac;
        TESCO(2, nq) = (static_cast<odepack_cpp_real>(nqp1)) / ELCO(1, nq);
        TESCO(3, nq) = (static_cast<odepack_cpp_real>(nq + 2)) / ELCO(1, nq);
        rq1fac /= fnq;
    } // end for nq
    return;
//
#ifdef ELCO
#undef ELCO
#endif
#ifdef TESCO
#undef TESCO
#endif
}


/**
 * PURPOSE  Interpolate solution derivatives.
 * AUTHOR  Hindmarsh, Alan C., (LLNL)
 * DESCRIPTION
 * 
 * DINTDY computes interpolated values of the K-th derivative of the
 * dependent variable vector y, and stores it in DKY.  This routine
 * is called within the package with K = 0 and T = TOUT, but may
 * also be called by the user for any K up to the current order.
 * (See detailed instructions in the usage documentation.)
 * 
 * The computed values in DKY are gotten by interpolation using the
 * Nordsieck history array YH.  This array corresponds uniquely to a
 * vector-valued polynomial of degree NQCUR or less, and DKY is set
 * to the K-th derivative of this polynomial at T.
 * The formula for DKY is:
 *              q
 *  DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
 *             j=K
 * where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
 * The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are
 * communicated by COMMON.  The above sum is done in reverse order.
 * IFLAG is returned negative if either K or T is out of bounds.
 */
void Odepack::DINTDY(odepack_cpp_real t, int k, odepack_cpp_real *yh, int nyh, odepack_cpp_real *dky, int &iflag) {
#ifndef YH
#define YH(i, j) ARRAY2D(yh, nyh, i, j)
#endif
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
//
    int i, ic, j, jb, jb2, jj, jj1, jp1;
    odepack_cpp_real c, r, s, tp;
    std::string msg;
//
//***FIRST EXECUTABLE STATEMENT  DINTDY
    iflag = 0;
    if (k < 0 || k > nq) goto LABEL_80;
    tp = tn - hu - 100.0 * uround * (std::abs(tn) + std::abs(hu)) * (hu >= 0.0 ? 1.0 : -1.0);
    if ((t - tp) * (t - tn) > 0.0) goto LABEL_90;
//
    s = (t - tn) / h;
    ic = 1;
    if (k == 0) goto LABEL_15;
    jj1 = l - k;
    for (jj = jj1; jj <= nq; ++jj) {
        ic *= jj;
    }
LABEL_15:
    c = static_cast<odepack_cpp_real>(ic);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(dky, i) = c * YH(i, l);
    }
    if (k == nq) goto LABEL_55;
    jb2 = nq - k;
    for (jb = 1; jb <= jb2; ++jb) {
        j = nq - jb;
        jp1 = j + 1;
        ic = 1;
        if (k == 0) goto LABEL_35;
        jj1 = jp1 - k;
        for (jj = jj1; jj <= j; ++jj) {
            ic *= jj;
        }
LABEL_35:
        c = static_cast<odepack_cpp_real>(ic);
        for (i = 1; i <= n; ++i) {
            ARRAY1D(dky, i) = c * YH(i, jp1) + s * ARRAY1D(dky, i);
        }
    }
    if (k == 0) return;
LABEL_55:
    r = std::pow(h, -static_cast<odepack_cpp_real>(k));
    for (i = 1; i <= n; ++i) {
        ARRAY1D(dky, i) *= r;
    }
    return;
//
LABEL_80:
    msg = "DINTDY-  K (=I1) illegal      ";
    XERRWD(msg, 30, 51, 0, 1, k, 0, 0, 0.0, 0.0);
    iflag = -1;
    return;
LABEL_90:
    msg = "DINTDY-  T (=R1) illegal      ";
    XERRWD(msg, 30, 51, 0, 1, k, 0, 0, 0.0, 0.0);
    msg = "      T not in interval TCUR - HU (= R1) to TCUR (=R2)      ";
    XERRWD(msg, 60, 52, 0, 0, 0, 0, 2, tp, tn);
    iflag = -2;
    return;
//
#ifdef YH
#undef YH
#endif
}


/**
 * PURPOSE  Compute and process Newton iteration matrix.
 * AUTHOR  Hindmarsh, Alan C., (LLNL)
 * DESCRIPTION
 * 
 * DPREPJ is called by DSTODE to compute and process the matrix
 * P = I - h*el(1)*J , where J is an approximation to the Jacobian.
 * Here J is computed by the user-supplied routine JAC if
 * MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
 * If MITER = 3, a diagonal approximation to J is used.
 * J is stored in WM and replaced by P.  If MITER .ne. 3, P is then
 * subjected to LU decomposition in preparation for later solution
 * of linear systems with P as coefficient matrix.  This is done
 * by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
 * 
 * In addition to variables described in DSTODE and DLSODE prologues,
 * communication with DPREPJ uses the following:
 * Y     = array containing predicted values on entry.
 * FTEM  = work array of length N (ACOR in DSTODE).
 * SAVF  = array containing f evaluated at predicted y.
 * WM    = real work space for matrices.  On output it contains the
 *         inverse diagonal matrix if MITER = 3 and the LU decomposition
 *         of P if MITER is 1, 2 , 4, or 5.
 *         Storage of matrix elements starts at WM(3).
 *         WM also contains the following matrix-related data:
 *         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
 *         WM(2) = H*EL0, saved for later use if MITER = 3.
 * IWM   = integer work space containing pivot information, starting at
 *         IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
 *         parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
 * EL0   = EL(1) (input).
 * IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
 *         P matrix found to be singular.
 * JCUR  = output flag = 1 to indicate that the Jacobian matrix
 *         (or approximation) is now current.
 * This routine also uses the COMMON variables EL0, H, TN, UROUND,
 * MITER, N, NFE, and NJE.
 */
void Odepack::DPREPJ(int neq,  odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *ewt, 
    odepack_cpp_real *ftem, odepack_cpp_real *savf, odepack_cpp_real *wm, int *iwm, ODEPACK_FUNCTION f, 
    ODEPACK_JACOBIAN1 jac, void *user_data) {
#ifndef YH
#define YH(i, j) ARRAY2D(yh, nyh, i, j)
#endif
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
//
    int i, i1, i2, ier, ii, j, j1, jj, lenp,
        mba, mband, meb1, meband, ml, ml3, mu, np1;
    odepack_cpp_real con, di, fac, hl0, r, r0, srur, yi, yj, yjj;
//
//***FIRST EXECUTABLE STATEMENT  DPREPJ
    nje++;
    ierpj = 0;
    jcur = 1;
    hl0 = h * el0;
    if (miter == 1) {
        goto LABEL_100;
    } else if (miter == 2) {
        goto LABEL_200;
    } else if (miter == 3) {
        goto LABEL_300;
    } else if (miter == 4) {
        goto LABEL_400;
    } else if (miter == 5) {
        goto LABEL_500;
    }
// If MITER = 1, call JAC and multiply by scalar. -----------------------
LABEL_100:
    lenp = n * n;
    for (i = 1; i <= lenp; ++i) {
        ARRAY1D(wm, i + 2) = 0.0;
    }
    (*jac)(neq, tn, y, 0, 0, &ARRAY1D(wm, 3), n, user_data);
    con = -hl0;
    for (i = 1; i <= lenp; ++i) {
        ARRAY1D(wm, i + 2) *= con;
    }
    goto LABEL_240;
// If MITER = 2, make N calls to F to approximate J. --------------------
LABEL_200:
    fac = DVNORM(n, savf, ewt);
    r0 = 1000.0 * std::abs(h) * uround * static_cast<odepack_cpp_real>(n) * fac;
    if (r0 == 0.0) r0 = 1.0;
    srur = ARRAY1D(wm, 1);
    j1 = 2;
    for (j = 1; j <= n; ++j) {
        yj = ARRAY1D(y, j);
        r = std::max(srur * std::abs(yj), r0 / ARRAY1D(ewt, j));
        ARRAY1D(y, j) += r;
        fac = -hl0 / r;
        (*f)(neq, tn, y, ftem, user_data);
        for (i = 1; i <= n; ++i) {
            ARRAY1D(wm, i + j1) = (ARRAY1D(ftem, i) - ARRAY1D(savf, i)) * fac;
        }
        ARRAY1D(y, j) = yj;
        j1 += n;
    }
    nfe += n;
// Add identity matrix. -------------------------------------------------
LABEL_240:
    j = 3;
    np1 = n + 1;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(wm, j) += 1.0;
        j += np1;
    }
// Do LU decomposition on P. --------------------------------------------
    DGEFA(&ARRAY1D(wm, 3), n, n, &ARRAY1D(iwm, 21), ier);
    if (ier != 0) ierpj = 1;
    return;
// If MITER = 3, construct a diagonal approximation to J and P. ---------
LABEL_300:
    ARRAY1D(wm, 2) = hl0;
    r = el0 * 0.1;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = ARRAY1D(y, i) + r * (h * ARRAY1D(savf, i) - YH(i, 2));
    }
    (*f)(neq, tn, y, &ARRAY1D(wm, 3), user_data);
    nfe++;
    for (i = 1; i <= n; ++i) {
        r0 = h * ARRAY1D(savf, i) - YH(i, 2);
        di = 0.1 * r0 - h * (ARRAY1D(wm, i + 2) - ARRAY1D(savf, i));
        ARRAY1D(wm, i + 2) = 1.0;
        if (std::abs(r0) < uround / ARRAY1D(ewt, i)) continue;
        if (std::abs(di) == 0.0) goto LABEL_330;
        ARRAY1D(wm, i + 2) = 0.1 * r0 / di;
    }
    return;
LABEL_330:
    ierpj = 1;
    return;
// If MITER = 4, call JAC and multiply by scalar. -----------------------
LABEL_400:
    ml = ARRAY1D(iwm, 1);
    mu = ARRAY1D(iwm, 2);
    ml3 = ml + 3;
    mband = ml + mu + 1;
    meband = mband + ml;
    lenp = meband * n;
    for (i = 1; i <= lenp; ++i) {
        ARRAY1D(wm, i + 2) = 0.0;
    }
    (*jac)(neq, tn, y, ml, mu, &ARRAY1D(wm, ml3), meband, user_data);
    con = -hl0;
    for (i = 1; i <= lenp; ++i) {
        ARRAY1D(wm, i + 2) *= con;
    }
    goto LABEL_570;
// If MITER = 5, make MBAND calls to F to approximate J. ----------------
LABEL_500:
    ml = ARRAY1D(iwm, 1);
    mu = ARRAY1D(iwm, 2);
    mband = ml + mu + 1;
    mba = std::min(mband, n);
    meband = mband + ml;
    meb1 = meband - 1;
    srur = ARRAY1D(wm, 1);
    fac = DVNORM(n, savf, ewt);
    r0 = 1000.0 * std::abs(h) * uround * static_cast<odepack_cpp_real>(n) * fac;
    if (r0 == 0.0) r0 = 1.0;
    for (j = 1; j <= mba; ++j) {
        for (i = j; i <= n; i += mband) {
            yi = ARRAY1D(y, i);
            r = std::max(srur * std::abs(yi), r0 / ARRAY1D(ewt, i));
            ARRAY1D(y, i) += r;
        }
        (*f)(neq, tn, y, ftem, user_data);
        for (jj = j; jj <= n; jj += mband) {
            ARRAY1D(y, jj) = YH(jj, 1);
            yjj = ARRAY1D(y, jj);
            r = std::max(srur * std::abs(yjj), r0 / ARRAY1D(ewt, jj));
            fac = -hl0 / r;
            i1 = std::max(jj - mu, 1);
            i2 = std::min(jj + ml, n);
            ii = jj * meb1 - ml + 2;
            for (i = i1; i <= i2; ++i) {
                ARRAY1D(wm, ii + i) = (ARRAY1D(ftem, i) - ARRAY1D(savf, i)) * fac;
            }
        }
    }
    nfe += mba;
// Add identity matrix. -------------------------------------------------
LABEL_570:
    ii = mband + 2;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(wm, ii) += 1.0;
        ii += meband;
    }
// Do LU decomposition of P. --------------------------------------------
    DGBFA(&ARRAY1D(wm, 3), meband, n, ml, mu, &ARRAY1D(iwm, 21), ier);
    if (ier != 0) ierpj = 1;
    return;
//
#ifdef YH
#undef YH
#endif
}


/**
 * PURPOSE  ODEPACK linear system solver.
 * AUTHOR  Hindmarsh, Alan C., (LLNL)
 * DESCRIPTION
 * 
 * This routine manages the solution of the linear system arising from
 * a chord iteration.  It is called if MITER .ne. 0.
 * If MITER is 1 or 2, it calls DGESL to accomplish this.
 * If MITER = 3 it updates the coefficient h*EL0 in the diagonal
 * matrix, and then computes the solution.
 * If MITER is 4 or 5, it calls DGBSL.
 * Communication with DSOLSY uses the following variables:
 * WM    = real work space containing the inverse diagonal matrix if
 *         MITER = 3 and the LU decomposition of the matrix otherwise.
 *         Storage of matrix elements starts at WM(3).
 *         WM also contains the following matrix-related data:
 *         WM(1) = SQRT(UROUND) (not used here),
 *         WM(2) = HL0, the previous value of h*EL0, used if MITER = 3.
 * IWM   = integer work space containing pivot information, starting at
 *         IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
 *         parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
 * X     = the right-hand side vector on input, and the solution vector
 *         on output, of length N.
 * TEM   = vector of work space of length N, not used in this version.
 * IERSL = output flag (in COMMON).  IERSL = 0 if no trouble occurred.
 *         IERSL = 1 if a singular matrix arose with MITER = 3.
 * This routine also uses the COMMON variables EL0, H, MITER, and N.
 */
void Odepack::DSOLSY(odepack_cpp_real *wm, int *iwm, odepack_cpp_real *x, odepack_cpp_real *tem) {
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
//
    int i, meband, ml, mu;
    odepack_cpp_real di, hl0, phl0, r;
//
//***FIRST EXECUTABLE STATEMENT  DSOLSY
    iersl = 0;
    if (miter == 1 || miter == 2) {
        goto LABEL_100;
    } else if (miter == 3) {
        goto LABEL_300;
    } else if (miter == 4 || miter == 5) {
        goto LABEL_400;
    }
LABEL_100:
    DGESL(&ARRAY1D(wm, 3), n, n, &ARRAY1D(iwm, 21), x, 0);
    return;
//
LABEL_300:
    phl0 = ARRAY1D(wm, 2);
    hl0 = h * el0;
    ARRAY1D(wm, 2) = hl0;
    if (hl0 == phl0) goto LABEL_330;
    r = hl0 / phl0;
    for (i = 1; i <= n; ++i) {
        di = 1.0 - r * (1.0 - 1.0 / ARRAY1D(wm, i + 2));
        if (std::abs(di) == 0.0) goto LABEL_390;
        ARRAY1D(wm, i + 2) = 1.0 / di;
    }
LABEL_330:
    for (i = 1; i <= n; ++i) {
        ARRAY1D(x, i) *= ARRAY1D(wm, i + 2);
    }
    return;
LABEL_390:
    iersl = 1;
    return;
//
LABEL_400:
    ml = ARRAY1D(iwm, 1);
    mu = ARRAY1D(iwm, 2);
    meband = 2 * ml + mu + 1;
    DGBSL(&ARRAY1D(wm, 3), meband, n, ml, mu, &ARRAY1D(iwm, 21), x, 0);
    return;
}

/**
 * 
 * PURPOSE  Save/restore ODEPACK COMMON blocks.
 * AUTHOR  Hindmarsh, Alan C., (LLNL)
 * DESCRIPTION
 * 
 * This routine saves or restores (depending on JOB) the contents of
 * the COMMON block DLS001, which is used internally
 * by one or more ODEPACK solvers.
 * 
 * RSAV = real array of length 218 or more.
 * ISAV = integer array of length 37 or more.
 * JOB  = flag indicating to save or restore the COMMON blocks:
 *        JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV)
 *        JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV)
 *        A call with JOB = 2 presumes a prior call with JOB = 1.
 */
void Odepack::DSRCOM(odepack_cpp_real *rsav, int *isav, int job) {
//
    int i;
// DLS001
    odepack_cpp_real *rls = dls1_.rls;
    int *ils = dls1_.ils;
//
    int lenrls = 218;
    int lenils = 37;
//
//***FIRST EXECUTABLE STATEMENT  DSRCOM
    if (job == 2) goto LABEL_100;
//
    for (i = 1; i <= lenrls; ++i) {
        ARRAY1D(rsav, i) = ARRAY1D(rls, i);
    }
    for (i = 1; i <= lenils; ++i) {
        ARRAY1D(isav, i) = ARRAY1D(ils, i);
    }
    return;
//
LABEL_100:
    for (i = 1; i <= lenrls; ++i) {
        ARRAY1D(rls, i) = ARRAY1D(rsav, i);
    }
    for (i = 1; i <= lenils; ++i) {
        ARRAY1D(ils, i) = ARRAY1D(isav, i);
    }
    return;
}

/**
 * @fn DSTODE
 * 
 * PURPOSE  Performs one step of an ODEPACK integration.
 * AUTHOR  Hindmarsh, Alan C., (LLNL)
 * DESCRIPTION
 * 
 * DSTODE performs one step of the integration of an initial value
 * problem for a system of ordinary differential equations.
 * Note:  DSTODE is independent of the value of the iteration method
 * indicator MITER, when this is .ne. 0, and hence is independent
 * of the type of chord method used, or the Jacobian structure.
 * Communication with DSTODE is done with the following variables:
 * 
 * NEQ    = integer array containing problem size in NEQ(1), and
 *          passed as the NEQ argument in all calls to F and JAC.
 * Y      = an array of length .ge. N used as the Y argument in
 *          all calls to F and JAC.
 * YH     = an NYH by LMAX array containing the dependent variables
 *          and their approximate scaled derivatives, where
 *          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
 *          j-th derivative of y(i), scaled by h**j/factorial(j)
 *          (j = 0,1,...,NQ).  on entry for the first step, the first
 *          two columns of YH must be set from the initial values.
 * NYH    = a constant integer .ge. N, the first dimension of YH.
 * YH1    = a one-dimensional array occupying the same space as YH.
 * EWT    = an array of length N containing multiplicative weights
 *          for local error measurements.  Local errors in Y(i) are
 *          compared to 1.0/EWT(i) in various error tests.
 * SAVF   = an array of working storage, of length N.
 *          Also used for input of YH(*,MAXORD+2) when JSTART = -1
 *          and MAXORD .lt. the current order NQ.
 * ACOR   = a work array of length N, used for the accumulated
 *          corrections.  On a successful return, ACOR(i) contains
 *          the estimated one-step local error in Y(i).
 * WM,IWM = real and integer work arrays associated with matrix
 *          operations in chord iteration (MITER .ne. 0).
 * PJAC   = name of routine to evaluate and preprocess Jacobian matrix
 *          and P = I - h*el0*JAC, if a chord method is being used.
 * SLVS   = name of routine to solve linear system in chord iteration.
 * CCMAX  = maximum relative change in h*el0 before PJAC is called.
 * H      = the step size to be attempted on the next step.
 *          H is altered by the error control algorithm during the
 *          problem.  H can be either positive or negative, but its
 *          sign must remain constant throughout the problem.
 * HMIN   = the minimum absolute value of the step size h to be used.
 * HMXI   = inverse of the maximum absolute value of h to be used.
 *          HMXI = 0.0 is allowed and corresponds to an infinite hmax.
 *          HMIN and HMXI may be changed at any time, but will not
 *          take effect until the next change of h is considered.
 * TN     = the independent variable. TN is updated on each step taken.
 * JSTART = an integer used for input only, with the following
 *          values and meanings:
 *               0  perform the first step.
 *           .gt.0  take a new step continuing from the last.
 *              -1  take the next step with a new value of H, MAXORD,
 *                    N, METH, MITER, and/or matrix parameters.
 *              -2  take the next step with a new value of H,
 *                    but with other inputs unchanged.
 *          On return, JSTART is set to 1 to facilitate continuation.
 * KFLAG  = a completion code with the following meanings:
 *               0  the step was succesful.
 *              -1  the requested error could not be achieved.
 *              -2  corrector convergence could not be achieved.
 *              -3  fatal error in PJAC or SLVS.
 *          A return with KFLAG = -1 or -2 means either
 *          abs(H) = HMIN or 10 consecutive failures occurred.
 *          On a return with KFLAG negative, the values of TN and
 *          the YH array are as of the beginning of the last
 *          step, and H is the last step size attempted.
 * MAXORD = the maximum order of integration method to be allowed.
 * MAXCOR = the maximum number of corrector iterations allowed.
 * MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
 * MXNCF  = maximum number of convergence failures allowed.
 * METH/MITER = the method flags.  See description in driver.
 * N      = the number of first-order differential equations.
 * The values of CCMAX, H, HMIN, HMXI, TN, JSTART, KFLAG, MAXORD,
 * MAXCOR, MSBP, MXNCF, METH, MITER, and N are communicated via COMMON.
 */
template <typename ODEPACK_JACOBIAN>
void Odepack::DSTODE(int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *yh1, odepack_cpp_real *ewt, odepack_cpp_real *savf, odepack_cpp_real *acor,
    odepack_cpp_real *wm, void *iwm_in, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN jac, FUNC_PJAC<ODEPACK_JACOBIAN> pjac,
    FUNC_SLVS slvs, void *user_data) {
#ifndef YH
#define YH(i, j) ARRAY2D(yh, nyh, i, j)
#endif
#ifndef ELCO
#define ELCO(i, j) ARRAY2D(elco, 13, i, j)
#endif
#ifndef TESCO
#define TESCO(i, j) ARRAY2D(tesco, 3, i, j)
#endif
//
    int i, i1, iredo, iret, j, jb, m, ncf, newq;
    odepack_cpp_real dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup,
        r, rh, rhdn, rhsm, rhup, told;
// DLS001
    odepack_cpp_real &conit = dls1_.conit, &crate = dls1_.crate, *el = dls1_.el, *elco = dls1_.elco,
        &hold = dls1_.hold, &rmax = dls1_.rmax, *tesco = dls1_.tesco,
        &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &ialth = dls1_.ialth, &ipup = dls1_.ipup, &lmax = dls1_.lmax, &meo = dls1_.meo, &nqnyh = dls1_.nqnyh, &nslp = dls1_.nslp,
        &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
//
    int *iwm = static_cast<int*>(iwm_in);
//
//***FIRST EXECUTABLE STATEMENT  DSTODE
    kflag = 0;
    told = tn;
    ncf = 0;
    ierpj = 0;
    iersl = 0;
    jcur = 0;
    icf = 0;
    delp = 0.0;
    if (jstart > 0) goto LABEL_200;
    if (jstart == -1) goto LABEL_100;
    if (jstart == -2) goto LABEL_160;
//-----------------------------------------------------------------------
// On the first call, the order is set to 1, and other variables are
// initialized.  RMAX is the maximum ratio by which H can be increased
// in a single step.  It is initially 1.E4 to compensate for the small
// initial H, but then is normally equal to 10.  If a failure
// occurs (in corrector convergence or error test), RMAX is set to 2
// for the next increase.
//-----------------------------------------------------------------------
    lmax = maxord + 1;
    nq = 1;
    l = 2;
    ialth = 2;
    rmax = 10000.0;
    rc = 0.0;
    el0 = 1.0;
    crate = 0.7;
    hold = h;
    meo = meth;
    nslp = 0;
    ipup = miter;
    iret = 3;
    goto LABEL_140;
//-----------------------------------------------------------------------
// The following block handles preliminaries needed when JSTART = -1.
// IPUP is set to MITER to force a matrix update.
// If an order increase is about to be considered (IALTH = 1),
// IALTH is reset to 2 to postpone consideration one more step.
// If the caller has changed METH, DCFODE is called to reset
// the coefficients of the method.
// If the caller has changed MAXORD to a value less than the current
// order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
// If H is to be changed, YH must be rescaled.
// If H or METH is being changed, IALTH is reset to L = NQ + 1
// to prevent further changes in H for that many steps.
//-----------------------------------------------------------------------
LABEL_100:
    ipup = miter;
    lmax = maxord + 1;
    if (ialth == 1) ialth = 2;
    if (meth == meo) goto LABEL_110;
    DCFODE(meth, elco, tesco);
    meo = meth;
    if (nq > maxord) goto LABEL_120;
    ialth = l;
    iret = 1;
    goto LABEL_150;
LABEL_110:
    if (nq <= maxord) goto LABEL_160;
LABEL_120:
    nq = maxord;
    l = lmax;
    for (i = 1; i <= l; ++i) {
        ARRAY1D(el, i) = ELCO(i, nq);
    }
    nqnyh = nq * nyh;
    rc = rc * ARRAY1D(el, 1) / el0;
    el0 = ARRAY1D(el, 1);
    conit = 0.5 / static_cast<odepack_cpp_real>(nq + 2);
    ddn = DVNORM(n, savf, ewt) / TESCO(1, l);
    exdn = 1.0 / static_cast<odepack_cpp_real>(l);
    rhdn = 1.0 / (1.3 * std::pow(ddn, exdn) + 0.0000013);
    rh = std::min(rhdn, 1.0);
    iredo = 3;
    if (h == hold) goto LABEL_170;
    rh = std::min(rh, std::abs(h / hold));
    h = hold;
    goto LABEL_175;
//-----------------------------------------------------------------------
// DCFODE is called to get all the integration coefficients for the
// current METH.  Then the EL vector and related constants are reset
// whenever the order NQ is changed, or at the start of the problem.
//-----------------------------------------------------------------------
LABEL_140:
    DCFODE(meth, elco, tesco);
LABEL_150:
    for (i = 1; i <= l; ++i) {
        ARRAY1D(el, i) = ELCO(i, nq);
    }
    nqnyh = nq * nyh;
    rc = rc * ARRAY1D(el, 1) / el0;
    el0 = ARRAY1D(el, 1);
    conit = 0.5 / static_cast<odepack_cpp_real>(nq + 2);
    if (iret == 1) {
        goto LABEL_160;
    } else if (iret == 2) {
        goto LABEL_170;
    } else if (iret == 3) {
        goto LABEL_200;
    }
//-----------------------------------------------------------------------
// If H is being changed, the H ratio RH is checked against
// RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
// L = NQ + 1 to prevent a change of H for that many steps, unless
// forced by a convergence or error test failure.
//-----------------------------------------------------------------------
LABEL_160:
    if (h == hold) goto LABEL_200;
    rh = h / hold;
    h = hold;
    iredo = 3;
    goto LABEL_175;
LABEL_170:
    rh = std::max(rh, hmin / std::abs(h));
LABEL_175:
    rh = std::min(rh, rmax);
    rh = rh / std::max(1.0, std::abs(h) * hmxi * rh);
    r = 1.0;
    for (j = 2; j <= l; ++j) {
        r *= rh;
        for (i = 1; i <= n; ++i) {
            YH(i, j) *= r;
        }
    }
    h *= rh;
    rc *= rh;
    ialth = l;
    if (iredo == 0) goto LABEL_690;
//-----------------------------------------------------------------------
// This section computes the predicted values by effectively
// multiplying the YH array by the Pascal Triangle matrix.
// RC is the ratio of new to old values of the coefficient  H*EL(1).
// When RC differs from 1 by more than CCMAX, IPUP is set to MITER
// to force PJAC to be called, if a Jacobian is involved.
// In any case, PJAC is called at least every MSBP steps.
//-----------------------------------------------------------------------
LABEL_200:
    if (std::abs(rc - 1.0) > ccmax) ipup = miter; 
    if (nst >= nslp + msbp) ipup = miter;
    tn += h;
    i1 = nqnyh + 1;
    for (jb = 1; jb <= nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= nqnyh; ++i) {
            ARRAY1D(yh1, i) += ARRAY1D(yh1, i + nyh);
        }
    }
//-----------------------------------------------------------------------
// Up to MAXCOR corrector iterations are taken.  A convergence test is
// made on the R.M.S. norm of each correction, weighted by the error
// weight vector EWT.  The sum of the corrections is accumulated in the
// vector ACOR(i).  The YH array is not altered in the corrector loop.
//-----------------------------------------------------------------------
LABEL_220:
    m = 0;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = YH(i, 1);
    }
    (*f)(neq, tn, y, savf, user_data);
    nfe++;
    if (ipup <= 0) goto LABEL_250;
//-----------------------------------------------------------------------
// If indicated, the matrix P = I - h*el(1)*J is reevaluated and
// preprocessed before starting the corrector iteration.  IPUP is set
// to 0 as an indicator that this has been done.
//-----------------------------------------------------------------------
    (this->*pjac)(neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac, user_data);
    ipup = 0;
    rc = 1.0;
    nslp = nst;
    crate = 0.7;
    if (ierpj != 0) goto LABEL_430;
LABEL_250:
    for (i = 1; i <= n; ++i) {
        ARRAY1D(acor, i) = 0.0;
    }
LABEL_270:
    if (miter != 0) goto LABEL_350;
//-----------------------------------------------------------------------
// In the case of functional iteration, update Y directly from
// the result of the last function evaluation.
//-----------------------------------------------------------------------
    for (i = 1; i <= n; ++i) {
        ARRAY1D(savf, i) = h * ARRAY1D(savf, i) - YH(i, 2);
        ARRAY1D(y, i) = ARRAY1D(savf, i) - ARRAY1D(acor, i);
    }
    del = DVNORM(n, y, ewt);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = YH(i, 1) + ARRAY1D(el, 1) * ARRAY1D(savf, i);
        ARRAY1D(acor, i) = ARRAY1D(savf, i);
    }
    goto LABEL_400;
//-----------------------------------------------------------------------
// In the case of the chord method, compute the corrector error,
// and solve the linear system with that as right-hand side and
// P as coefficient matrix.
//-----------------------------------------------------------------------
LABEL_350:
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = h * ARRAY1D(savf, i) - (YH(i, 2) + ARRAY1D(acor, i));
    }
    (this->*slvs)(wm, iwm, y, savf);
    if (iersl < 0) goto LABEL_430;
    if (iersl > 0) goto LABEL_410;
    del = DVNORM(n, y, ewt);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(acor, i) += ARRAY1D(y, i);
        ARRAY1D(y, i) = YH(i, 1) + ARRAY1D(el, 1) * ARRAY1D(acor, i);
    }
//-----------------------------------------------------------------------
// Test for convergence.  If M.gt.0, an estimate of the convergence
// rate constant is stored in CRATE, and this is used in the test.
//-----------------------------------------------------------------------
LABEL_400:
    if (m != 0) crate = std::max(0.2 * crate, del / delp);
    dcon = del * std::min(1.0, 1.5 * crate) / (TESCO(2, nq) * conit);
    if (dcon <= 1.0) goto LABEL_450;
    m++;
    if (m == maxcor) goto LABEL_410;
    if (m >= 2 && del > 2.0 * delp) goto LABEL_410;
    delp = del;
    (*f)(neq, tn, y, savf, user_data);
    nfe++;
    goto LABEL_270;
//-----------------------------------------------------------------------
// The corrector iteration failed to converge.
// If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
// the next try.  Otherwise the YH array is retracted to its values
// before prediction, and H is reduced, if possible.  If H cannot be
// reduced or MXNCF failures have occurred, exit with KFLAG = -2.
//-----------------------------------------------------------------------
LABEL_410:
    if (miter == 0 || jcur == 1) goto LABEL_430;
    icf = 1;
    ipup = miter;
    goto LABEL_220;
LABEL_430:
    icf = 2;
    ncf++;
    rmax = 2.0;
    tn = told;
    i1 = nqnyh + 1;
    for (jb = 1; jb <= nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= nqnyh; ++i) {
            ARRAY1D(yh1, i) -= ARRAY1D(yh1, i + nyh);
        }
    }
    if (ierpj < 0 || iersl < 0) goto LABEL_680;
    if (std::abs(h) <= hmin * 1.00001) goto LABEL_670;
    if (ncf == mxncf) goto LABEL_670;
    rh = 0.25;
    ipup = miter;
    iredo = 1;
    goto LABEL_170;
//-----------------------------------------------------------------------
// The corrector has converged.  JCUR is set to 0
// to signal that the Jacobian involved may need updating later.
// The local error test is made and control passes to statement 500
// if it fails.
//-----------------------------------------------------------------------
LABEL_450:
    jcur = 0;
    if (m == 0) dsm = del / TESCO(2, nq);
    if (m > 0) dsm = DVNORM(n, acor, ewt) / TESCO(2, nq);
    if (dsm > 1.0) goto LABEL_500;
//-----------------------------------------------------------------------
// After a successful step, update the YH array.
// Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
// If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
// use in a possible order increase on the next step.
// If a change in H is considered, an increase or decrease in order
// by one is considered also.  A change in H is made only if it is by a
// factor of at least 1.1.  If not, IALTH is set to 3 to prevent
// testing for that many steps.
//-----------------------------------------------------------------------
    kflag = 0;
    iredo = 0;
    nst++;
    hu = h;
    nqu = nq;
    for (j = 1; j <= l; ++j) {
        for (i = 1; i <= n; ++i) {
            YH(i, j) += ARRAY1D(el, j) * ARRAY1D(acor, i);
        }
    }
    ialth--;
    if (ialth == 0) goto LABEL_520;
    if (ialth > 1) goto LABEL_700;
    if (l == lmax) goto LABEL_700;
    for (i = 1; i <= n; ++i) {
        YH(i, lmax) = ARRAY1D(acor, i);
    }
    goto LABEL_700;
//-----------------------------------------------------------------------
// The error test failed.  KFLAG keeps track of multiple failures.
// Restore TN and the YH array to their previous values, and prepare
// to try the step again.  Compute the optimum step size for this or
// one lower order.  After 2 or more failures, H is forced to decrease
// by a factor of 0.2 or less.
//-----------------------------------------------------------------------
LABEL_500:
    kflag--;
    tn = told;
    i1 = nqnyh + 1;
    for (jb = 1; jb <= nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= nqnyh; ++i) {
            ARRAY1D(yh1, i) -= ARRAY1D(yh1, i + nyh);
        }
    }
    rmax = 2.0;
    if (std::abs(h) <= hmin * 1.00001) goto LABEL_660;
    if (kflag <= -3) goto LABEL_640;
    iredo = 2;
    rhup = 0.0;
    goto LABEL_540;
//-----------------------------------------------------------------------
// Regardless of the success or failure of the step, factors
// RHDN, RHSM, and RHUP are computed, by which H could be multiplied
// at order NQ - 1, order NQ, or order NQ + 1, respectively.
// In the case of failure, RHUP = 0.0 to avoid an order increase.
// The largest of these is determined and the new order chosen
// accordingly.  If the order is to be increased, we compute one
// additional scaled derivative.
//-----------------------------------------------------------------------
LABEL_520:
    rhup = 0.0;
    if (l == lmax) goto LABEL_540;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(savf, i) = ARRAY1D(acor, i) - YH(i, lmax);
    }
    dup = DVNORM(n, savf, ewt) / TESCO(3, nq);
    exup = 1.0 / static_cast<odepack_cpp_real>(l + 1);
    rhup = 1.0 / (1.4 * std::pow(dup, exup) + 0.0000014);
LABEL_540:
    exsm = 1.0 / static_cast<odepack_cpp_real>(l);
    rhsm = 1.0 / (1.2 * std::pow(dsm, exsm) + 0.0000012);
    rhdn = 0.0;
    if (nq == 1) goto LABEL_560;
    ddn = DVNORM(n, &YH(1, l), ewt) / TESCO(1, nq);
    exdn = 1.0 / static_cast<odepack_cpp_real>(nq);
    rhdn = 1.0 / (1.3 * std::pow(ddn, exdn) + 0.0000013);
LABEL_560:
    if (rhsm >= rhup) goto LABEL_570;
    if (rhup > rhdn) goto LABEL_590;
    goto LABEL_580;
LABEL_570:
    if (rhsm < rhdn) goto LABEL_580;
    newq = nq;
    rh = rhsm;
    goto LABEL_620;
LABEL_580:
    newq = nq - 1;
    rh = rhdn;
    if (kflag < 0 && rh > 1.0) rh = 1.0;
    goto LABEL_620;
LABEL_590:
    newq = l;
    rh = rhup;
    if (rh < 1.1) goto LABEL_610;
    r = ARRAY1D(el, l) / static_cast<odepack_cpp_real>(l);
    for (i = 1; i <= n; ++i) {
        YH(i, newq + 1) = ARRAY1D(acor, i) * r;
    }
    goto LABEL_630;
LABEL_610:
    ialth = 3;
    goto LABEL_700;
LABEL_620:
    if (kflag == 0 && rh < 1.1) goto LABEL_610;
    if (kflag <= -2) rh = std::min(rh, 0.2);
//-----------------------------------------------------------------------
// If there is a change of order, reset NQ, l, and the coefficients.
// In any case H is reset according to RH and the YH array is rescaled.
// Then exit from 690 if the step was OK, or redo the step otherwise.
//-----------------------------------------------------------------------
    if (newq == nq) goto LABEL_170;
LABEL_630:
    nq = newq;
    l = nq + 1;
    iret = 2;
    goto LABEL_150;
//-----------------------------------------------------------------------
// Control reaches this section if 3 or more failures have occured.
// If 10 failures have occurred, exit with KFLAG = -1.
// It is assumed that the derivatives that have accumulated in the
// YH array have errors of the wrong order.  Hence the first
// derivative is recomputed, and the order is set to 1.  Then
// H is reduced by a factor of 10, and the step is retried,
// until it succeeds or H reaches HMIN.
//-----------------------------------------------------------------------
LABEL_640:
    if (kflag == -10) goto LABEL_660;
    rh = 0.1;
    rh = std::max(hmin / std::abs(h), rh);
    h *= rh;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = YH(i, 1);
    }
    (*f)(neq, tn, y, savf, user_data);
    nfe++;
    for (i = 1; i <= n; ++i) {
        YH(i, 2) = h * ARRAY1D(savf, i);
    }
    ipup = miter;
    ialth = 5;
    if (nq == 1) goto LABEL_200;
    nq = 1;
    l = 2;
    iret = 3;
    goto LABEL_150;
//-----------------------------------------------------------------------
// All returns are made through this section.  H is saved in HOLD
// to allow the caller to change H on the next step.
//-----------------------------------------------------------------------
LABEL_660:
    kflag = -1;
    goto LABEL_720;
LABEL_670:
    kflag = -2;
    goto LABEL_720;
LABEL_680:
    kflag = -3;
    goto LABEL_720;
LABEL_690:
    rmax = 10.0;
LABEL_700:
    r = 1.0 / TESCO(2, nqu);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(acor, i) *= r;
    }
LABEL_720:
    hold = h;
    jstart = 1;
    return;
//
#ifdef YH
#undef YH
#endif
#ifdef ELCO
#undef ELCO
#endif
#ifdef TESCO
#undef TESCO
#endif
}

// Explicit instantiation
template void Odepack::DSTODE<ODEPACK_JACOBIAN1>(
    int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real*yh1, odepack_cpp_real *ewt, odepack_cpp_real *savf, odepack_cpp_real *acor,
    odepack_cpp_real *wm, void *iwm_in, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN1 jac, FUNC_PJAC<ODEPACK_JACOBIAN1> pjac, FUNC_SLVS slvs, void *user_data
);

template void Odepack::DSTODE<ODEPACK_JACOBIAN2>(
    int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real*yh1, odepack_cpp_real *ewt, odepack_cpp_real *savf, odepack_cpp_real *acor,
    odepack_cpp_real *wm, void *iwm_in, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN2 jac, FUNC_PJAC<ODEPACK_JACOBIAN2> pjac, FUNC_SLVS slvs, void *user_data
);


/**
 * @fn DEWSET
 * 
 * PURPOSE  Set error weight vector.
 * AUTHOR  Hindmarsh, Alan C., (LLNL)
 * DESCRIPTION
 * 
 * This subroutine sets the error weight vector EWT according to
 *     EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
 * with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
 * depending on the value of ITOL.
 */
void Odepack::DEWSET(
    int n, int itol, odepack_cpp_real *rtol, odepack_cpp_real *atol, odepack_cpp_real *ycur, odepack_cpp_real *ewt
)
{
    int i;
//**FIRST EXECUTABLE STATEMENT  DEWSET
    if (itol == 1) {
        goto LABEL_10;
    } else if (itol == 2) {
        goto LABEL_20;
    } else if (itol == 3) {
        goto LABEL_30;
    } else if (itol == 4) {
        goto LABEL_40;
    }
LABEL_10:
    for (i = 1; i <= n; ++i) {
        ARRAY1D(ewt, i) = ARRAY1D(rtol, 1) * std::abs(ARRAY1D(ycur, i)) + ARRAY1D(atol, 1);
    }
    return;  
LABEL_20:
    for (i = 1; i <= n; ++i) {   
        ARRAY1D(ewt, i) = ARRAY1D(rtol, 1) * std::abs(ARRAY1D(ycur, i)) + ARRAY1D(atol, i);
    }
    return;
LABEL_30:
    for (i = 1; i <= n; ++i) {   
        ARRAY1D(ewt, i) = ARRAY1D(rtol, i) * std::abs(ARRAY1D(ycur, i)) + ARRAY1D(atol, 1);
    }
    return;
LABEL_40:
    for (i = 1; i <= n; ++i) {   
        ARRAY1D(ewt, i) = ARRAY1D(rtol, i) * std::abs(ARRAY1D(ycur, i)) + ARRAY1D(atol, i);
    }
    return;
}


/**
 * @fn DVNORM
 * 
 * PURPOSE  Weighted root-mean-square vector norm.
 * AUTHOR  Hindmarsh, Alan C., (LLNL)
 * DESCRIPTION
 * 
 * This function routine computes the weighted root-mean-square norm
 * of the vector of length N contained in the array V, with weights
 * contained in the array W of length N:
 *   DVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 )
 */
odepack_cpp_real Odepack::DVNORM(
    int n, odepack_cpp_real *v, odepack_cpp_real *w
)
{
    odepack_cpp_real sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += v[i] * v[i] * w[i] * w[i];
    }
    return std::sqrt(sum / static_cast<odepack_cpp_real>(n));
}


/**
 * @fn DIPREP
 * 
 * This routine serves as an interface between the driver and
 * Subroutine DPREP.  It is called only if MITER is 1 or 2.
 * Tasks performed here are:
 *  * call DPREP,
 *  * reset the required WM segment length LENWK,
 *  * move YH back to its final location (following WM in RWORK),
 *  * reset pointers for YH, SAVF, EWT, and ACOR, and
 *  * move EWT to its new position if ISTATE = 1.
 * IPFLAG is an output error indication flag.  IPFLAG = 0 if there was
 * no trouble, and IPFLAG is the value of the DPREP error flag IPPER
 * if there was trouble in Subroutine DPREP.
 */
void Odepack::DIPREP(
    int neq, odepack_cpp_real *y, odepack_cpp_real *rwork, int *ia, int *ja,
    int &ipflag, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN2 jac, void *user_data
)
{
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
// DLSS01
    int &iplost = dlss_.iplost, &iesp = dlss_.iesp, &istatc = dlss_.istatc, &iys = dlss_.iys, &iba = dlss_.iba, &ibian = dlss_.ibian, &ibjan = dlss_.ibjan, &ibjgp = dlss_.ibjgp,
        &ipian = dlss_.ipian, &ipjan = dlss_.ipjan, &ipjgp = dlss_.ipjgp, &ipigp = dlss_.ipigp, &ipr = dlss_.ipr, &ipc = dlss_.ipc, &ipic = dlss_.ipic, &ipisp = dlss_.ipisp, &iprsp = dlss_.iprsp, &ipa = dlss_.ipa,
        &lenyh = dlss_.lenyh, &lenyhm = dlss_.lenyhm, &lenwk = dlss_.lenwk, &lreq = dlss_.lreq, &lrat = dlss_.lrat, &lrest = dlss_.lrest, &lwmin = dlss_.lwmin, &moss = dlss_.moss, &msbj = dlss_.msbj,
        &nslj = dlss_.nslj, &ngp = dlss_.ngp, &nlu = dlss_.nlu, &nnz = dlss_.nnz, &nsp = dlss_.nsp, &nzl = dlss_.nzl, &nzu = dlss_.nzu;
//
    int i, imax, lewtn, lyhd, lyhn;
//
    ipflag = 0;
// Call DPREP to do matrix preprocessing operations. --------------------
    DPREP(neq, y, &ARRAY1D(rwork, lyh), &ARRAY1D(rwork, lsavf), &ARRAY1D(rwork, lewt),
        &ARRAY1D(rwork, lacor), ia, ja, &ARRAY1D(rwork, lwm), &ARRAY1D(rwork, lwm), ipflag, f, jac, user_data);
    lenwk = std::max(lreq, lwmin);
    if (ipflag < 0) return;
// If DPREP was successful, move YH to end of required space for WM. ----
    lyhn = lwm + lenwk;
    if (lyhn > lyh) return;
    lyhd = lyh - lyhn;
    if (lyhd == 0) goto LABEL_20;
    imax = lyhn - 1 + lenyhm;
    for (i = lyhn; i <= imax; ++i) {
        ARRAY1D(rwork, i) = ARRAY1D(rwork, i + lyhd);
    }
    lyh = lyhn;
// Reset pointers for SAVF, EWT, and ACOR. ------------------------------
LABEL_20:
    lsavf = lyh + lenyh;
    lewtn = lsavf + n;
    lacor = lewtn + n;
    if (istatc == 3) goto LABEL_40;
// If ISTATE = 1, move EWT (left) to its new position. ------------------
    if (lewtn > lewt) return;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(rwork, i + lewtn - 1) = ARRAY1D(rwork, i + lewt - 1);
    }
LABEL_40:
    lewt = lewtn;
    return;
}


/**
 * @fn DPREP
 * 
 * This routine performs preprocessing related to the sparse linear
 * systems that must be solved if MITER = 1 or 2.
 * The operations that are performed here are:
 *  * compute sparseness structure of Jacobian according to MOSS,
 *  * compute grouping of column indices (MITER = 2),
 *  * compute a new ordering of rows and columns of the matrix,
 *  * reorder JA corresponding to the new ordering,
 *  * perform a symbolic LU factorization of the matrix, and
 *  * set pointers for segments of the IWK/WK array.
 * In addition to variables described previously, DPREP uses the
 * following for communication:
 * YH     = the history array.  Only the first column, containing the
 *          current Y vector, is used.  Used only if MOSS .ne. 0.
 * SAVF   = a work array of length NEQ, used only if MOSS .ne. 0.
 * EWT    = array of length NEQ containing (inverted) error weights.
 *          Used only if MOSS = 2 or if ISTATE = MOSS = 1.
 * FTEM   = a work array of length NEQ, identical to ACOR in the driver,
 *          used only if MOSS = 2.
 * WK     = a real work array of length LENWK, identical to WM in
 *          the driver.
 * IWK    = integer work array, assumed to occupy the same space as WK.
 * LENWK  = the length of the work arrays WK and IWK.
 * ISTATC = a copy of the driver input argument ISTATE (= 1 on the
 *          first call, = 3 on a continuation call).
 * IYS    = flag value from ODRV or CDRV.
 * IPPER  = output error flag with the following values and meanings:
 *          0  no error.
 *         -1  insufficient storage for internal structure pointers.
 *         -2  insufficient storage for JGROUP.
 *         -3  insufficient storage for ODRV.
 *         -4  other error flag from ODRV (should never occur).
 *         -5  insufficient storage for CDRV.
 *         -6  other error flag from CDRV.
 */
void Odepack::DPREP(
    int neq, odepack_cpp_real *y, odepack_cpp_real *yh, odepack_cpp_real *savf, odepack_cpp_real *ewt, odepack_cpp_real *ftem, int *ia, int *ja,
    odepack_cpp_real *wk, void *iwk_in, int &ipper, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN2 jac, void *user_data
)
{
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
// DLSS01
    odepack_cpp_real &con0 = dlss_.con0, &conmin = dlss_.conmin, &ccmxj = dlss_.ccmxj, &psmall = dlss_.psmall, &rbig = dlss_.rbig, &seth = dlss_.seth;
    int &iplost = dlss_.iplost, &iesp = dlss_.iesp, &istatc = dlss_.istatc, &iys = dlss_.iys, &iba = dlss_.iba, &ibian = dlss_.ibian, &ibjan = dlss_.ibjan, &ibjgp = dlss_.ibjgp,
        &ipian = dlss_.ipian, &ipjan = dlss_.ipjan, &ipjgp = dlss_.ipjgp, &ipigp = dlss_.ipigp, &ipr = dlss_.ipr, &ipc = dlss_.ipc, &ipic = dlss_.ipic, &ipisp = dlss_.ipisp, &iprsp = dlss_.iprsp, &ipa = dlss_.ipa,
        &lenyh = dlss_.lenyh, &lenyhm = dlss_.lenyhm, &lenwk = dlss_.lenwk, &lreq = dlss_.lreq, &lrat = dlss_.lrat, &lrest = dlss_.lrest, &lwmin = dlss_.lwmin, &moss = dlss_.moss, &msbj = dlss_.msbj,
        &nslj = dlss_.nslj, &ngp = dlss_.ngp, &nlu = dlss_.nlu, &nnz = dlss_.nnz, &nsp = dlss_.nsp, &nzl = dlss_.nzl, &nzu = dlss_.nzu;
//
    int i, ibr, ier, ipil, ipiu, iptt1, iptt2, j, jfound, k, 
        knew, kmax, kmin, ldif, lenigp, liwk, maxg, np1, nzsut;
    odepack_cpp_real dq, dyj, erwt, fac, yj;
//
    int *iwk = static_cast<int*>(iwk_in);
//
    ibian = lrat * 2;
    ipian = ibian + 1;
    np1 = n + 1;
    ipjan = ipian + np1;
    ibjan = ipjan - 1;
    liwk = lenwk * lrat;
    if (ipjan + n - 1 > liwk) goto LABEL_210;
    if (moss == 0) goto LABEL_30;
//
    if (istatc == 3) goto LABEL_20;
// ISTATE = 1 and MOSS .ne. 0.  Perturb Y for structure determination. --
    for (i = 1; i <= n; ++i) {
        erwt = 1.0 / ARRAY1D(ewt, i);
        fac = 1.0 + 1.0 / (static_cast<odepack_cpp_real>(i) + 1.0);
        ARRAY1D(y, i) += fac * std::abs(erwt) * (ARRAY1D(y, i) >= 0.0 ? 1.0 : -1.0);
    }
    if (moss == 1) {
        goto LABEL_70;
    } else if (moss == 2) {
        goto LABEL_100;
    }
//
LABEL_20:
// ISTATE = 3 and MOSS .ne. 0.  Load Y from YH(*,1). --------------------
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = ARRAY1D(yh, i);
    }
    if (moss == 1) {
        goto LABEL_70;
    } else if (moss == 2) {
        goto LABEL_100;
    }
// MOSS = 0.  Process user's IA,JA.  Add diagonal entries if necessary. -
LABEL_30:
    knew = ipjan;
    kmin = ARRAY1D(ia, 1);
    ARRAY1D(iwk, ipian) = 1;
    for (j = 1; j <= n; ++j) {
        jfound = 0;
        kmax = ARRAY1D(ia, j + 1) - 1;
        if (kmin > kmax) goto LABEL_45;
        for (k = kmin; k <= kmax; ++k) {
            i = ARRAY1D(ja, k);
            if (i == j) jfound = 1;
            if (knew > liwk) goto LABEL_210;
            ARRAY1D(iwk, knew) = i;
            knew++;
        }
        if (jfound == 1) goto LABEL_50;
LABEL_45:
        if (knew > liwk) goto LABEL_210;
        ARRAY1D(iwk, knew) = j;
        knew++;
LABEL_50:
        ARRAY1D(iwk, ipian + j) = knew + 1 - ipjan;
        kmin = kmax + 1;
    }
    goto LABEL_140;
//
// MOSS = 1.  Compute structure from user-supplied Jacobian routine JAC.
LABEL_70:
// A dummy call to F allows user to create temporaries for use in JAC. --
    (*f)(neq, tn, y, savf, user_data);
    k = ipjan;
    ARRAY1D(iwk, ipian) = 1;
    for (j = 1; j <= n; ++j) {
        if (k > liwk) goto LABEL_210;
        ARRAY1D(iwk, k) = j;
        k++;
        for (i = 1; i <= n; ++i) {
            ARRAY1D(savf, i) = 0.0;
        }
        (*jac)(neq, tn, y, j, &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan), savf, user_data);
        for (i = 1; i <= n; ++i) {
            if (std::abs(ARRAY1D(savf, i)) <= seth) continue;
            if (i == j) continue;
            if (k > liwk) goto LABEL_210;
            ARRAY1D(iwk, k) = i;
            k++;
        }
        ARRAY1D(iwk, ipian + j) = k + 1 - ipjan;
    }
    goto LABEL_140;
//
// MOSS = 2.  Compute structure from results of N + 1 calls to F. -------
LABEL_100:
    k = ipjan;
    ARRAY1D(iwk, ipian) = 1;
    (*f)(neq, tn, y, savf, user_data);
    for (j = 1; j <= n; ++j) {
        if (k > liwk) goto LABEL_210;
        ARRAY1D(iwk, k) = j;
        k++;
        yj = ARRAY1D(y, j);
        erwt = 1.0 / ARRAY1D(ewt, j);
        dyj = std::abs(erwt) * (yj >= 0.0 ? 1.0 : -1.0);
        ARRAY1D(y, j) = yj + dyj;
        (*f)(neq, tn, y, ftem, user_data);
        ARRAY1D(y, j) = yj;
        for (i = 1; i <= n; ++i) {
            dq = (ARRAY1D(ftem, i) - ARRAY1D(savf, i)) / dyj;
            if (std::abs(dq) <= seth) continue;
            if (i == j) continue;
            if (k > liwk) goto LABEL_210;
            ARRAY1D(iwk, k) = i;
            k++;
        }
        ARRAY1D(iwk, ipian + j) = k + 1 - ipjan;
    }
//
LABEL_140:
    if (moss == 0 || istatc != 1) goto LABEL_150;
// If ISTATE = 1 and MOSS .ne. 0, restore Y from YH. --------------------
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = ARRAY1D(yh, i);
    }
LABEL_150:
    nnz = ARRAY1D(iwk, ipian + n) - 1;
    lenigp = 0;
    ipigp = ipjan + nnz;
    if (miter != 2) goto LABEL_160;
//
// Compute grouping of column indices (MITER = 2). ----------------------
    maxg = np1;
    ipjgp = ipjan + nnz;
    ibjgp = ipjgp - 1;
    ipigp = ipjgp + n;
    iptt1 = ipigp + np1;
    iptt2 = iptt1 + n;
    lreq = iptt2 + n - 1;
    if (lreq > liwk) goto LABEL_220;
    JGROUP(n, &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan), maxg, ngp, &ARRAY1D(iwk, ipigp), 
        &ARRAY1D(iwk, ipjgp), &ARRAY1D(iwk, iptt1), &ARRAY1D(iwk, iptt2), ier);
    if (ier != 0) goto LABEL_220;
    lenigp = ngp + 1;
//
// Compute new ordering of rows/columns of Jacobian. --------------------
LABEL_160:
    ipr = ipigp + lenigp;
    ipc = ipr;
    ipic = ipc + n;
    ipisp = ipic + n;
    iprsp = (ipisp - 2) / lrat + 2;
    iesp = lenwk + 1 - iprsp;
    if (iesp < 0) goto LABEL_230;
    ibr = ipr - 1;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(iwk, ibr + i) = i;
    }
    nsp = liwk + 1 - ipisp;
    ODRV(n, &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan), wk, &ARRAY1D(iwk, ipr), &ARRAY1D(iwk, ipic), 
        nsp, &ARRAY1D(iwk, ipisp), 1, iys);
    if (iys == 11 * n + 1) goto LABEL_240;
    if (iys != 0) goto LABEL_230;
//
// Reorder JAN and do symbolic LU factorization of matrix. --------------
    ipa = lenwk + 1 - nnz;
    nsp = ipa - iprsp;
    lreq = std::max(12 * n / lrat, 6 * n / lrat + 2 * n + nnz) + 3;
    lreq = lreq + iprsp - 1 + nnz;
    if (lreq > lenwk) goto LABEL_250;
    iba = ipa - 1;
    for (i = 1; i <= nnz; ++i) {
        ARRAY1D(wk, iba + i) = 0.0;
    }
    ipisp = lrat * (iprsp - 1) + 1;
    CDRV(n, &ARRAY1D(iwk, ipr), &ARRAY1D(iwk, ipc), &ARRAY1D(iwk, ipic), &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan),
        &ARRAY1D(wk, ipa), &ARRAY1D(wk, ipa), &ARRAY1D(wk, ipa), nsp, &ARRAY1D(iwk, ipisp), &ARRAY1D(wk, iprsp), iesp, 5, iys);
    lreq = lenwk - iesp;
    if (iys == 10 * n + 1) goto LABEL_250;
    if (iys != 0) goto LABEL_260;
    ipil = ipisp;
    ipiu = ipil + 2 * n + 1;
    nzu = ARRAY1D(iwk, ipil + n) - ARRAY1D(iwk, ipil);
    nzl = ARRAY1D(iwk, ipiu + n) - ARRAY1D(iwk, ipiu);
    if (lrat > 1) goto LABEL_190;
    ADJLR(n, &ARRAY1D(iwk, ipisp), ldif);
    lreq += ldif;
LABEL_190:
    if (lrat == 2 && nnz == n) lreq++;
    nsp = nsp + lreq - lenwk;
    ipa = lreq + 1 - nnz;
    iba = ipa  - 1;
    ipper = 0;
    return;
//
LABEL_210:
    ipper = -1;
    lreq = 2 + (2 * n + 1) / lrat;
    lreq = std::max(lenwk + 1, lreq);
    return;
//
LABEL_220:
    ipper = -2;
    lreq = (lreq - 1) / lrat + 1;
    return;
//
LABEL_230:
    ipper = -3;
    CNTNZU(n, &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan), nzsut);
    lreq = lenwk - iesp + (3 * n + 4 * nzsut - 1) / lrat + 1;
    return;
//
LABEL_240:
    ipper = -4;
    return;
//
LABEL_250:
    ipper = -5;
    return;
//
LABEL_260:
    ipper = -6;
    lreq = lenwk;
    return;
}


/**
 * @fn JGROUP
 * 
 * This subroutine constructs groupings of the column indices of
 * the Jacobian matrix, used in the numerical evaluation of the
 * Jacobian by finite differences.
 * 
 * Input:
 * N      = the order of the matrix.
 * IA,JA  = sparse structure descriptors of the matrix by rows.
 * MAXG   = length of available storage in the IGP array.
 * 
 * Output:
 * NGRP   = number of groups.
 * JGP    = array of length N containing the column indices by groups.
 * IGP    = pointer array of length NGRP + 1 to the locations in JGP
 *          of the beginning of each group.
 * IER    = error indicator.  IER = 0 if no error occurred, or 1 if
 *          MAXG was insufficient.
 * 
 * INCL and JDONE are working arrays of length N.
 */
void Odepack::JGROUP(
    int n, int *ia, int *ja, int maxg, int &ngrp, int *igp, int *jgp, int *incl, int *jdone, int &ier
)
{
    int i, j, k, kmin, kmax, ncol, ng;
    bool is_goto_50 = false;
//
    ier = 0;
    for (j = 1; j <= n; ++j) {
        ARRAY1D(jdone, j) = 0;
    }
    ncol = 1;
    for (ng = 1; ng <= maxg; ++ng) {
        ARRAY1D(igp, ng) = ncol;
        for (i = 1; i <= n; ++i) {
            ARRAY1D(incl, i) = 0;
        }
        for (j = 1; j <= n; ++j) {
            is_goto_50 = false;
// Reject column J if it is already in a group.--------------------------
            if (ARRAY1D(jdone, j) == 1) continue;
            kmin = ARRAY1D(ia, j);
            kmax = ARRAY1D(ia, j + 1) - 1;
            for (k = kmin; k <= kmax; ++k) {
// Reject column J if it overlaps any column already in this group.------
                i = ARRAY1D(ja, k);
                if (ARRAY1D(incl, i) == 1) {
                    is_goto_50 = true;
                    break;
                };
            }
            if (is_goto_50) continue;
// Accept column J into group NG.----------------------------------------
            ARRAY1D(jgp, ncol) = j;
            ncol++;
            ARRAY1D(jdone, j) = 1;
            for (k = kmin; k <= kmax; ++k) {
                i = ARRAY1D(ja, k);
                ARRAY1D(incl, i) = 1;
            }
        }
// Stop if this group is empty (grouping is complete).-------------------
        if (ncol == ARRAY1D(igp, ng)) goto LABEL_70;
    }
// Error return if not all columns were chosen (MAXG too small).---------
    if (ncol <= n) goto LABEL_80;
    ng = maxg;
LABEL_70:
    ngrp = ng - 1;
    return;
LABEL_80:
    ier = 1;
    return;
}


/**
 * @fn ADJLR
 * 
 * This routine computes an adjustment, LDIF, to the required
 * integer storage space in IWK (sparse matrix work space).
 * It is called only if the word length ratio is LRAT = 1.
 * This is to account for the possibility that the symbolic LU phase
 * may require more storage than the numerical LU and solution phases.
 */
void Odepack::ADJLR(
    int n, int *isp, int &ldif
)
{
    int ip, jlmax, jumax, lnfc, lsfc, nzlu;
//
    ip = 2 * n + 1;
// Get JLMAX = IJL(N) and JUMAX = IJU(N) (sizes of JL and JU). ----------
    jlmax = ARRAY1D(isp, ip);
    jumax = ARRAY1D(isp, ip + ip);
// NZLU = (size of L) + (size of U) = (IL(N+1)-IL(1)) + (IU(N+1)-IU(1)).
    nzlu = ARRAY1D(isp, n + 1) - ARRAY1D(isp, 1) + ARRAY1D(isp, ip + n + 1) - ARRAY1D(isp, ip + 1);
    lsfc = 12 * n + 3 + 2 * std::max(jlmax, jumax);
    lnfc = 9 * n + 2 + jlmax + jumax + nzlu;
    ldif = std::max(0, lsfc - lnfc);
    return;
}


/**
 * @fn CNTNZU
 * 
 * This routine counts the number of nonzero elements in the strict
 * upper triangle of the matrix M + M(transpose), where the sparsity
 * structure of M is given by pointer arrays IA and JA.
 * This is needed to compute the storage requirements for the
 * sparse matrix reordering operation in ODRV.
 */
void Odepack::CNTNZU(
    int n, int *ia, int *ja, int &nzsut
)
{
    int ii, jj, j, jmin, jmax, k, kmin, kmax, num;
    bool is_goto_40 = false;
//
    num = 0;
    for (ii = 1; ii <= n; ++ii) {
        jmin = ARRAY1D(ia, ii);
        jmax = ARRAY1D(ia, ii + 1) - 1;
        if (jmin > jmax) continue;
        for (j = jmin; j <= jmax; ++j) {
            if (ARRAY1D(ja, j) < ii) {
                goto LABEL_10;
            } else if (ARRAY1D(ja, j) == ii) {
                continue;
            } else {
                goto LABEL_30;
            }
LABEL_10:
            jj = ARRAY1D(ja,  j);
            kmin = ARRAY1D(ia, jj);
            kmax = ARRAY1D(ia, jj + 1) - 1;
            if (kmin > kmax) goto LABEL_30;
            is_goto_40 = false;
            for (k = kmin; k <= kmax; ++k) {
                if (ARRAY1D(ja, k) == ii) {
                    is_goto_40 = true;
                    break;
                }
            }
            if (is_goto_40) continue;
LABEL_30:
            num++;
        }
    }
    nzsut = num;
    return;
}


/**
 * @fn DPRJS
 * 
 * DPRJS is called to compute and process the matrix
 * P = I - H*EL(1)*J , where J is an approximation to the Jacobian.
 * J is computed by columns, either by the user-supplied routine JAC
 * if MITER = 1, or by finite differencing if MITER = 2.
 * if MITER = 3, a diagonal approximation to J is used.
 * if MITER = 1 or 2, and if the existing value of the Jacobian
 * (as contained in P) is considered acceptable, then a new value of
 * P is reconstructed from the old value.  In any case, when MITER
 * is 1 or 2, the P matrix is subjected to LU decomposition in CDRV.
 * P and its LU decomposition are stored (separately) in WK.
 * 
 * In addition to variables described previously, communication
 * with DPRJS uses the following:
 * Y     = array containing predicted values on entry.
 * FTEM  = work array of length N (ACOR in DSTODE).
 * SAVF  = array containing f evaluated at predicted y.
 * WK    = real work space for matrices.  On output it contains the
 *         inverse diagonal matrix if MITER = 3, and P and its sparse
 *         LU decomposition if MITER is 1 or 2.
 *         Storage of matrix elements starts at WK(3).
 *         WK also contains the following matrix-related data:
 *         WK(1) = SQRT(UROUND), used in numerical Jacobian increments.
 *         WK(2) = H*EL0, saved for later use if MITER = 3.
 * IWK   = integer work space for matrix-related data, assumed to
 *         be equivalenced to WK.  In addition, WK(IPRSP) and IWK(IPISP)
 *         are assumed to have identical locations.
 * EL0   = EL(1) (input).
 * IERPJ = output error flag (in Common).
 *       = 0 if no error.
 *       = 1  if zero pivot found in CDRV.
 *       = 2  if a singular matrix arose with MITER = 3.
 *       = -1 if insufficient storage for CDRV (should not occur here).
 *       = -2 if other error found in CDRV (should not occur here).
 * JCUR  = output flag showing status of (approximate) Jacobian matrix:
 *          = 1 to indicate that the Jacobian is now current, or
 *          = 0 to indicate that a saved value was used.
 * This routine also uses other variables in Common.
 */
void Odepack::DPRJS(
    int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *ewt, odepack_cpp_real *ftem, odepack_cpp_real *savf, odepack_cpp_real *wk, int *iwk, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN2 jac, void *user_data
)
{
#ifndef YH
#define YH(i, j) ARRAY2D(yh, nyh, i, j)
#endif
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
// DLSS01
    odepack_cpp_real &con0 = dlss_.con0, &conmin = dlss_.conmin, &ccmxj = dlss_.ccmxj, &psmall = dlss_.psmall, &rbig = dlss_.rbig, &seth = dlss_.seth;
    int &iplost = dlss_.iplost, &iesp = dlss_.iesp, &istatc = dlss_.istatc, &iys = dlss_.iys, &iba = dlss_.iba, &ibian = dlss_.ibian, &ibjan = dlss_.ibjan, &ibjgp = dlss_.ibjgp,
        &ipian = dlss_.ipian, &ipjan = dlss_.ipjan, &ipjgp = dlss_.ipjgp, &ipigp = dlss_.ipigp, &ipr = dlss_.ipr, &ipc = dlss_.ipc, &ipic = dlss_.ipic, &ipisp = dlss_.ipisp, &iprsp = dlss_.iprsp, &ipa = dlss_.ipa,
        &lenyh = dlss_.lenyh, &lenyhm = dlss_.lenyhm, &lenwk = dlss_.lenwk, &lreq = dlss_.lreq, &lrat = dlss_.lrat, &lrest = dlss_.lrest, &lwmin = dlss_.lwmin, &moss = dlss_.moss, &msbj = dlss_.msbj,
        &nslj = dlss_.nslj, &ngp = dlss_.ngp, &nlu = dlss_.nlu, &nnz = dlss_.nnz, &nsp = dlss_.nsp, &nzl = dlss_.nzl, &nzu = dlss_.nzu;
//
    int i, imul, j, jj, jok, jmax, jmin, k, kmax, kmin, ng;
    odepack_cpp_real con, di, fac, hl0, pij, r, r0, rcon, rcont,
        srur;
//
    hl0 = h * el0;
    con = -hl0;
    if (miter == 3) goto LABEL_300;
// See whether J should be reevaluated (JOK = 0) or not (JOK = 1). ------
    jok = 1;
    if (nst == 0 || nst >= (nslj + msbj)) jok = 0;
    if (icf == 1 && std::abs(rc - 1.0) < ccmxj) jok = 0;
    if (icf == 2) jok = 0;
    if (jok == 1) goto LABEL_250;
//
// MITER = 1 or 2, and the Jacobian is to be reevaluated. ---------------
LABEL_20:
    jcur = 1;
    nje++;
    nslj = nst;
    iplost = 0;
    conmin = std::abs(con);
    if (miter == 1) {
        goto LABEL_100;
    } else if (miter == 2) {
        goto LABEL_200;
    }
//
// If MITER = 1, call JAC, multiply by scalar, and add identity. --------
LABEL_100:
    kmin = ARRAY1D(iwk, ipian);
    for (j = 1; j <= n; ++j) {
        kmax = ARRAY1D(iwk, ipian + j) - 1;
        for (i = 1; i <= n; ++i) {
            ARRAY1D(ftem, i) = 0.0;
        }
        (*jac)(neq, tn, y, j, &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan), ftem, user_data);
        for (k = kmin; k <= kmax; ++k) {
            i = ARRAY1D(iwk, ibjan + k);
            ARRAY1D(wk, iba + k) = ARRAY1D(ftem, i) * con;
            if (i == j) ARRAY1D(wk, iba + k) += 1.0;
        }
        kmin = kmax + 1;
    }
    goto LABEL_290;
//
// If MITER = 2, make NGP calls to F to approximate J and P. ------------
LABEL_200:
    fac = DVNORM(n, savf, ewt);
    r0 = 1000.0 * std::abs(h) * uround * static_cast<odepack_cpp_real>(n) * fac;
    if (r0 == 0.0) r0 = 1.0;
    srur = ARRAY1D(wk, 1);
    jmin = ARRAY1D(iwk, ipigp);
    for (ng = 1; ng <= ngp; ++ng) {
        jmax = ARRAY1D(iwk, ipigp + ng) - 1;
        for (j = jmin; j <= jmax; ++j) {
            jj = ARRAY1D(iwk, ibjgp + j);
            r  = std::max(srur * std::abs(ARRAY1D(y, jj)), r0 / ARRAY1D(ewt, jj));
            ARRAY1D(y, jj) += r;
        }
        (*f)(neq, tn, y, ftem, user_data);
        for (j = jmin; j <= jmax; ++j) {
            jj = ARRAY1D(iwk, ibjgp + j);
            ARRAY1D(y, jj) = YH(jj, 1);
            r = std::max(srur * std::abs(ARRAY1D(y, jj)), r0 / ARRAY1D(ewt, jj));
            fac  = - hl0 / r;
            kmin = ARRAY1D(iwk, ibian + jj);
            kmax = ARRAY1D(iwk, ibian + jj + 1) - 1;
            for (k = kmin; k <= kmax; ++k) {
                i = ARRAY1D(iwk, ibjan + k);
                ARRAY1D(wk, iba + k) = (ARRAY1D(ftem, i) - ARRAY1D(savf, i)) * fac;
                if (i == jj) ARRAY1D(wk, iba + k) += 1.0;
            }
        }
        jmin = jmax + 1;
    }
    nfe += ngp;
    goto LABEL_290;
//
// If JOK = 1, reconstruct new P from old P. ----------------------------
LABEL_250:
    jcur = 0;
    rcon = con / con0;
    rcont = std::abs(con) / conmin;
    if (rcont > rbig && iplost == 1) goto LABEL_20;
    kmin = ARRAY1D(iwk, ipian);
    for (j = 1; j <= n; ++j) {
        kmax = ARRAY1D(iwk, ipian + j) - 1;
        for (k = kmin; k <= kmax; ++k) {
            i = ARRAY1D(iwk, ibjan + k);
            pij = ARRAY1D(wk, iba + k);
            if (i != j) goto LABEL_260;
            pij -= 1.0;
            if (std::abs(pij) >= psmall) goto LABEL_260;
            iplost = 1;
            conmin = std::min(std::abs(con0), conmin);
LABEL_260:
            pij *= rcon;
            if (i == j) pij += 1.0;
            ARRAY1D(wk, iba + k) = pij;
        }
        kmin = kmax + 1;
    }
//
// Do numerical factorization of P matrix. ------------------------------
LABEL_290:
    nlu++;
    con0 = con;
    ierpj = 0;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(ftem, i) = 0.0;
    }
    CDRV(n, &ARRAY1D(iwk, ipr), &ARRAY1D(iwk, ipc), &ARRAY1D(iwk, ipic), &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan),
        &ARRAY1D(wk, ipa), ftem, ftem, nsp, &ARRAY1D(iwk, ipisp), &ARRAY1D(wk, iprsp), iesp, 2, iys);
    if (iys == 0) return;
    imul = (iys - 1) / n;
    ierpj = -2;
    if (imul == 8) ierpj = 1;
    if (imul == 10) ierpj = -1;
    return;
//
// If MITER = 3, construct a diagonal approximation to J and P. ---------
LABEL_300:
    jcur = 1;
    nje++;
    ARRAY1D(wk, 2) = hl0;
    ierpj = 0;
    r = el0 * 0.1;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = ARRAY1D(y, i) + r * (h * ARRAY1D(savf, i) - YH(i, 2));
    }
    (*f)(neq, tn, y, &ARRAY1D(wk, 3), user_data);
    nfe++;
    for (i = 1; i <= n; ++i) {
        r0 = h * ARRAY1D(savf, i) - YH(i, 2);
        di = 0.1 * r0 - h * (ARRAY1D(wk, i + 2) - ARRAY1D(savf, i));
        ARRAY1D(wk, i + 2) = 1.0;
        if (std::abs(r0) < uround / ARRAY1D(ewt, i)) continue;
        if (std::abs(di) == 0.0) goto LABEL_330;
        ARRAY1D(wk, i + 2) = 0.1 * r0 / di;
    }
    return;
LABEL_330:
    ierpj = 2;
    return;
//
#ifdef YH
#undef YH
#endif
}


/**
 * @fn DSOLSS
 * 
 * This routine manages the solution of the linear system arising from
 * a chord iteration.  It is called if MITER .ne. 0.
 * If MITER is 1 or 2, it calls CDRV to accomplish this.
 * If MITER = 3 it updates the coefficient H*EL0 in the diagonal
 * matrix, and then computes the solution.
 * communication with DSOLSS uses the following variables:
 * WK    = real work space containing the inverse diagonal matrix if
 *         MITER = 3 and the LU decomposition of the matrix otherwise.
 *         Storage of matrix elements starts at WK(3).
 *         WK also contains the following matrix-related data:
 *         WK(1) = SQRT(UROUND) (not used here),
 *         WK(2) = HL0, the previous value of H*EL0, used if MITER = 3.
 * IWK   = integer work space for matrix-related data, assumed to
 *         be equivalenced to WK.  In addition, WK(IPRSP) and IWK(IPISP)
 *         are assumed to have identical locations.
 * X     = the right-hand side vector on input, and the solution vector
 *         on output, of length N.
 * TEM   = vector of work space of length N, not used in this version.
 * IERSL = output flag (in Common).
 *         IERSL = 0  if no trouble occurred.
 *         IERSL = -1 if CDRV returned an error flag (MITER = 1 or 2).
 *                    This should never occur and is considered fatal.
 *         IERSL = 1  if a singular matrix arose with MITER = 3.
 * This routine also uses other variables in Common.
 */
void Odepack::DSOLSS(
    odepack_cpp_real *wk, int *iwk, odepack_cpp_real *x, odepack_cpp_real *tem
)
{
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
// DLSS01
    int &iplost = dlss_.iplost, &iesp = dlss_.iesp, &istatc = dlss_.istatc, &iys = dlss_.iys, &iba = dlss_.iba, &ibian = dlss_.ibian, &ibjan = dlss_.ibjan, &ibjgp = dlss_.ibjgp,
        &ipian = dlss_.ipian, &ipjan = dlss_.ipjan, &ipjgp = dlss_.ipjgp, &ipigp = dlss_.ipigp, &ipr = dlss_.ipr, &ipc = dlss_.ipc, &ipic = dlss_.ipic, &ipisp = dlss_.ipisp, &iprsp = dlss_.iprsp, &ipa = dlss_.ipa,
        &lenyh = dlss_.lenyh, &lenyhm = dlss_.lenyhm, &lenwk = dlss_.lenwk, &lreq = dlss_.lreq, &lrat = dlss_.lrat, &lrest = dlss_.lrest, &lwmin = dlss_.lwmin, &moss = dlss_.moss, &msbj = dlss_.msbj,
        &nslj = dlss_.nslj, &ngp = dlss_.ngp, &nlu = dlss_.nlu, &nnz = dlss_.nnz, &nsp = dlss_.nsp, &nzl = dlss_.nzl, &nzu = dlss_.nzu;
//
    int i;
    odepack_cpp_real di, hl0, phl0, r;
//
    iersl = 0;
    if (miter == 1 || miter == 2) {
        goto LABEL_100;
    } else if (miter == 3) {
        goto LABEL_300;
    }
LABEL_100:
    CDRV(n, &ARRAY1D(iwk, ipr), &ARRAY1D(iwk, ipc), &ARRAY1D(iwk, ipic), &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan),
        &ARRAY1D(wk, ipa), x, x, nsp, &ARRAY1D(iwk, ipisp), &ARRAY1D(wk, iprsp), iesp, 4, iersl);
    if (iersl != 0) iersl = -1;
    return;
//
LABEL_300:
    phl0 = ARRAY1D(wk, 2);
    hl0 = h * el0;
    ARRAY1D(wk, 2) = hl0;
    if (hl0 == phl0) goto LABEL_330;
    r = hl0 / phl0;
    for (i = 1; i <= n; ++i) {
        di = 1.0 - r * (1.0 - 1.0 / ARRAY1D(wk, i + 2));
        if (std::abs(di) == 0.0) goto LABEL_390;
        ARRAY1D(wk, i + 2) = 1.0 / di;
    }
LABEL_330:
    for (i = 1; i <= n; ++i) {
        ARRAY1D(x, i) *= ARRAY1D(wk, i + 2);
    }
    return;
LABEL_390:
    iersl = 1;
    return;
}


/**
 * @fn DSRCMS
 * 
 * This routine saves or restores (depending on JOB) the contents of
 * the Common blocks DLS001, DLSS01, which are used
 * internally by one or more ODEPACK solvers.
 * 
 * RSAV = real array of length 224 or more.
 * ISAV = integer array of length 71 or more.
 * JOB  = flag indicating to save or restore the Common blocks:
 *        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
 *        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
 *        A call with JOB = 2 presumes a prior call with JOB = 1.
 */
void Odepack::DSRCMS(
    odepack_cpp_real *rsav, int *isav, int job
)
{
    int i;
// DLS001
    odepack_cpp_real *rls = dls1_.rls;
    int *ils = dls1_.ils;
// DLSS01
    odepack_cpp_real *rlss = dlss_.rlss;
    int *ilss = dlss_.ilss;
//
    int lenrls = 218;
    int lenils = 37;
    int lenrss = 6;
    int leniss = 34;
//
    if (job == 2) goto LABEL_100;
    for (i = 1; i <= lenrls; ++i) {
        ARRAY1D(rsav, i) = ARRAY1D(rls, i);
    }
    for (i = 1; i <= lenrss; ++i) {
        ARRAY1D(rsav, lenrls + i) = ARRAY1D(rlss, i);
    }
//
    for (i = 1; i <= lenils; ++i) {
        ARRAY1D(isav, i) = ARRAY1D(ils, i);
    }
    for (i = 1; i <= leniss; ++i) {
        ARRAY1D(isav, lenils + i) = ARRAY1D(ilss, i);
    }
//
    return;
//
LABEL_100:
    for (i = 1; i <= lenrls; ++i) {
        ARRAY1D(rls, i) = ARRAY1D(rsav, i);
    }
    for (i = 1; i <= lenrss; ++i) {
        ARRAY1D(rlss, i) = ARRAY1D(rsav, lenrls + i);
    }
//
    for (i = 1; i <= lenils; ++i) {
        ARRAY1D(ils, i) = ARRAY1D(isav, i);
    }
    for (i = 1; i <= leniss; ++i) {
        ARRAY1D(ilss, i) = ARRAY1D(isav, lenils + i);
    }
//
    return;
}


/**
 * @fn ODRV -- driver for sparse matrix reordering routines
 * 
 * description
 * 
 * odrv finds a minimum degree ordering of the rows and columns
 * of a matrix m stored in (ia,ja,a) format (see below).  for the
 * reordered matrix, the work and storage required to perform
 * gaussian elimination is (usually) significantly less.
 * 
 * note.. odrv and its subordinate routines have been modified to
 * compute orderings for general matrices, not necessarily having any
 * symmetry.  the miminum degree ordering is computed for the
 * structure of the symmetric matrix  m + m-transpose.
 * modifications to the original odrv module have been made in
 * the coding in subroutine mdi, and in the initial comments in
 * subroutines odrv and md.
 * 
 * if only the nonzero entries in the upper triangle of m are being
 * stored, then odrv symmetrically reorders (ia,ja,a), (optionally)
 * with the diagonal entries placed first in each row.  this is to
 * ensure that if m(i,j) will be in the upper triangle of m with
 * respect to the new ordering, then m(i,j) is stored in row i (and
 * thus m(j,i) is not stored),  whereas if m(i,j) will be in the
 * strict lower triangle of m, then m(j,i) is stored in row j (and
 * thus m(i,j) is not stored).
 * 
 * storage of sparse matrices
 * 
 * the nonzero entries of the matrix m are stored row-by-row in the
 * array a.  to identify the individual nonzero entries in each row,
 * we need to know in which column each entry lies.  these column
 * indices are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
 * ja(k) = j.  to identify the individual rows, we need to know where
 * each row starts.  these row pointers are stored in the array ia.
 * i.e., if m(i,j) is the first nonzero entry (stored) in the i-th row
 * and  a(k) = m(i,j),  then  ia(i) = k.  moreover, ia(n+1) points to
 * the first location following the last element in the last row.
 * thus, the number of entries in the i-th row is  ia(i+1) - ia(i),
 * the nonzero entries in the i-th row are stored consecutively in
 *
 *         a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
 * 
 * and the corresponding column indices are stored consecutively in
 * 
 *         ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
 * 
 * when the coefficient matrix is symmetric, only the nonzero entries
 * in the upper triangle need be stored.  for example, the matrix
 * 
 *          ( 1  0  2  3  0 )
 *          ( 0  4  0  0  0 )
 *      m = ( 2  0  5  6  0 )
 *          ( 3  0  6  7  8 )
 *          ( 0  0  0  8  9 )
 * 
 * could be stored as
 * 
 *           - 1  2  3  4  5  6  7  8  9 10 11 12 13
 *        ---+--------------------------------------
 *        ia - 1  4  5  8 12 14
 *        ja - 1  3  4  2  1  3  4  1  3  4  5  4  5
 *         a - 1  2  3  4  2  5  6  3  6  7  8  8  9
 * 
 *   or (symmetrically) as
 * 
 *            - 1  2  3  4  5  6  7  8  9
 *         ---+--------------------------
 *         ia - 1  4  5  7  9 10
 *         ja - 1  3  4  2  3  4  4  5  5
 *          a - 1  2  3  4  5  6  7  8  9          .
 * 
 * 
 * parameters
 * 
 *   n    - order of the matrix
 * 
 *   ia   - integer one-dimensional array containing pointers to delimit
 *          rows in ja and a.  dimension = n+1
 * 
 *   ja   - integer one-dimensional array containing the column indices
 *          corresponding to the elements of a.  dimension = number of
 *          nonzero entries in (the upper triangle of) m
 * 
 *   a    - real one-dimensional array containing the nonzero entries in
 *          (the upper triangle of) m, stored by rows.  dimension =
 *          number of nonzero entries in (the upper triangle of) m
 * 
 *   p    - integer one-dimensional array used to return the permutation
 *          of the rows and columns of m corresponding to the minimum
 *          degree ordering.  dimension = n
 * 
 *   ip   - integer one-dimensional array used to return the inverse of
 *          the permutation returned in p.  dimension = n
 * 
 *   nsp  - declared dimension of the one-dimensional array isp.  nsp
 *          must be at least  3n+4k,  where k is the number of nonzeroes
 *          in the strict upper triangle of m
 * 
 *   isp  - integer one-dimensional array used for working storage.
 *          dimension = nsp
 * 
 *   path - integer path specification.  values and their meanings are -
 *            1  find minimum degree ordering only
 *            2  find minimum degree ordering and reorder symmetrically
 *                 stored matrix (used when only the nonzero entries in
 *                 the upper triangle of m are being stored)
 *            3  reorder symmetrically stored matrix as specified by
 *                 input permutation (used when an ordering has already
 *                 been determined and only the nonzero entries in the
 *                 upper triangle of m are being stored)
 *            4  same as 2 but put diagonal entries at start of each row
 *            5  same as 3 but put diagonal entries at start of each row
 * 
 *   flag - integer error flag.  values and their meanings are -
 *              0    no errors detected
 *             9n+k  insufficient storage in md
 *            10n+1  insufficient storage in odrv
 *            11n+1  illegal path specification
 * 
 * 
 * conversion from real to odepack_cpp_real precision
 *
 *   change the real declarations in odrv and sro to odepack_cpp_real precision
 *   declarations.
 */
void Odepack::ODRV(
    int n, int *ia, int *ja, odepack_cpp_real *a, int *p, int *ip, int nsp, int *isp, int path, int &flag
)
{
    int max, v, l, head, next, tmp, q;
    bool dflag;
//
//----initialize error flag and validate path specification
    flag = 0;
    if (path < 1 || 5 < path) goto LABEL_111;
//
//----allocate storage and find minimum degree ordering
    if ((path - 1) * (path - 2) * (path - 4) != 0) goto LABEL_1;
    max = (nsp - n) / 2;
    v = 1;
    l = v + max;
    head = l + max;
    next = head + n;
    if (max < n) goto LABEL_110;
//
    MD(n, ia, ja, max, &ARRAY1D(isp, v), &ARRAY1D(isp, l), &ARRAY1D(isp, head), p, ip, &ARRAY1D(isp, v), flag);
    if (flag != 0) goto LABEL_100;
//
//----allocate storage and symmetrically reorder matrix
LABEL_1:
    if ((path - 2) * (path - 3) * (path - 4) * (path - 5) != 0) goto LABEL_2;
    tmp = (nsp + 1) - n;
    q = tmp - (ARRAY1D(ia, n + 1) - 1);
    if (q < 1) goto LABEL_110;
//
    dflag = (path == 4 || path == 5);
    SRO(n, ip, ia, ja, a, &ARRAY1D(isp, tmp), &ARRAY1D(isp, q), dflag);
//
LABEL_2:
    return;
//
// ** error -- error detected in md
LABEL_100:
    return;
// ** error -- insufficient storage
LABEL_110:
    flag = 10 * n + 1;
    return;
// ** error -- illegal path specified
LABEL_111:
    flag = 11 * n + 1;
    return;
}


/**
 * @fn MD -- minimum degree algorithm (based on element model)
 * 
 * description
 * 
 *   md finds a minimum degree ordering of the rows and columns of a
 *   general sparse matrix m stored in (ia,ja,a) format.
 *   when the structure of m is nonsymmetric, the ordering is that
 *   obtained for the symmetric matrix  m + m-transpose.
 * 
 * additional parameters
 * 
 *   max  - declared dimension of the one-dimensional arrays v and l.
 *          max must be at least  n+2k,  where k is the number of
 *          nonzeroes in the strict upper triangle of m + m-transpose
 * 
 *   v    - integer one-dimensional work array.  dimension = max
 * 
 *   l    - integer one-dimensional work array.  dimension = max
 * 
 *   head - integer one-dimensional work array.  dimension = n
 * 
 *   last - integer one-dimensional array used to return the permutation
 *          of the rows and columns of m corresponding to the minimum
 *          degree ordering.  dimension = n
 * 
 *   next - integer one-dimensional array used to return the inverse of
 *          the permutation returned in last.  dimension = n
 * 
 *   mark - integer one-dimensional work array (may be the same as v).
 *          dimension = n
 * 
 *   flag - integer error flag.  values and their meanings are -
 *            0     no errors detected
 *            9n+k  insufficient storage in md
 * 
 * 
 * definitions of internal parameters
 * 
 *   ---------+---------------------------------------------------------
 *   v(s)     - value field of list entry
 *   ---------+---------------------------------------------------------
 *   l(s)     - link field of list entry  (0 =) end of list)
 *   ---------+---------------------------------------------------------
 *   l(vi)    - pointer to element list of uneliminated vertex vi
 *   ---------+---------------------------------------------------------
 *   l(ej)    - pointer to boundary list of active element ej
 *   ---------+---------------------------------------------------------
 *   head(d)  - vj =) vj head of d-list d
 *            -  0 =) no vertex in d-list d
 * 
 * 
 *            -                  vi uneliminated vertex
 *            -          vi in ek           -       vi not in ek
 *   ---------+-----------------------------+---------------------------
 *   next(vi) - undefined but nonnegative   - vj =) vj next in d-list
 *            -                             -  0 =) vi tail of d-list
 *   ---------+-----------------------------+---------------------------
 *   last(vi) - (not set until mdp)         - -d =) vi head of d-list d
 *            --vk =) compute degree        - vj =) vj last in d-list
 *            - ej =) vi prototype of ej    -  0 =) vi not in any d-list
 *            -  0 =) do not compute degree -
 *   ---------+-----------------------------+---------------------------
 *   mark(vi) - mark(vk)                    - nonneg. tag .lt. mark(vk)
 * 
 * 
 *            -                   vi eliminated vertex
 *            -      ei active element      -           otherwise
 *   ---------+-----------------------------+---------------------------
 *   next(vi) - -j =) vi was j-th vertex    - -j =) vi was j-th vertex
 *            -       to be eliminated      -       to be eliminated
 *   ---------+-----------------------------+---------------------------
 *   last(vi) -  m =) size of ei = m        - undefined
 *   ---------+-----------------------------+---------------------------
 *   mark(vi) - -m =) overlap count of ei   - undefined
 *            -       with ek = m           -
 *            - otherwise nonnegative tag   -
 *            -       .lt. mark(vk)         -
 * 
 * -----------------------------------------------------------------------
 */
void Odepack::MD(
    int n, int *ia, int *ja, int max, int *v, int *l, int *head, int *last, int *next, int *mark, int &flag
)
{
    int tag, k, dmin, vk, tail;
    int *ek = &vk; // equivalence  (vk,ek)
//
//----initialization
    tag = 0;
    MDI(n, ia, ja, max, v, l, head, last, next, mark, tag, flag);
    if (flag != 0) return;
//
    k = 0;
    dmin = 1;
//
//----while  k .lt. n  do
LABEL_1:
    if (k >= n) goto LABEL_4;
//
//------search for vertex of minimum degree
LABEL_2:
    if (ARRAY1D(head, dmin) > 0) goto LABEL_3;
    dmin++;
    goto LABEL_2;
//
//------remove vertex vk of minimum degree from degree list
LABEL_3:
    vk = ARRAY1D(head, dmin);
    ARRAY1D(head, dmin) = ARRAY1D(next, vk);
    if (ARRAY1D(head, dmin) > 0) ARRAY1D(last, ARRAY1D(head, dmin)) = -dmin;
//
//------number vertex vk, adjust tag, and tag vk
    k++;
    ARRAY1D(next,  vk) = -k;
    ARRAY1D(last, *ek) = dmin - 1;
    tag += ARRAY1D(last, *ek);
    ARRAY1D(mark,  vk) = tag;
//
//------form element ek from uneliminated neighbors of vk
    MDM(vk, tail, v, l, last, next, mark);
//
//------purge inactive elements and do mass elimination
    MDP(k, *ek, tail, v, l, head, last, next, mark);
//
//------update degrees of uneliminated vertices in ek
    MDU(*ek, dmin, v, l, head, last, next, mark);
//
    goto LABEL_1;
//
//----generate inverse permutation from permutation
LABEL_4:
    for (k = 1; k <= n; ++k) {
        ARRAY1D(next, k) = - ARRAY1D(next, k);
        ARRAY1D(last, ARRAY1D(next, k)) = k;
    }
//
    return;
}


/**
 * @fn MDI -- initialization
 */
void Odepack::MDI(
    int n, int *ia, int *ja, int max, int *v, int *l, int *head, int *last, int *next, int *mark, int tag, int &flag
)
{
    int sfs, j, jmin, jmax, vi, dvi, vj, lvk, k, kmax, nextvi;
    bool is_goto_5 = false;
//
//----initialize degrees, element lists, and degree lists
    for (vi = 1; vi <= n; ++vi) {
        ARRAY1D(mark, vi) = 1;
        ARRAY1D(l, vi) = 0;
        ARRAY1D(head, vi) = 0;
    }
    sfs = n + 1;
//
//----create nonzero structure
//----for each nonzero entry a(vi,vj)
    for (vi = 1; vi <= n; ++vi) {
        jmin = ARRAY1D(ia, vi);
        jmax = ARRAY1D(ia, vi + 1) - 1;
        if (jmin > jmax) continue;
        for (j = jmin; j <= jmax; ++j) {
            vj = ARRAY1D(ja, j);
            if (vj < vi) {
                goto LABEL_2;
            } if (vj == vi) {
                continue;
            } else {
                goto LABEL_4;
            }
//
//------if a(vi,vj) is in strict lower triangle
//------check for previous occurrence of a(vj,vi)
LABEL_2:
            lvk  = vi;
            kmax = ARRAY1D(mark, vi) - 1;
            if (kmax == 0) goto LABEL_4;
            is_goto_5 = false;
            for (k = 1; k <= kmax; ++k) {
                lvk = ARRAY1D(l, lvk);
                if (ARRAY1D(v, lvk) == vj) {
                    is_goto_5 = true;
                    break;
                }
            }
            if (is_goto_5) continue;
//----for unentered entries a(vi,vj)
LABEL_4:
            if (sfs >= max) goto LABEL_101;
//
//------enter vj in element list for vi
            ARRAY1D(mark, vi) += 1;
            ARRAY1D(v, sfs) = vj;
            ARRAY1D(l, sfs) = ARRAY1D(l, vi);
            ARRAY1D(l, vi) = sfs;
            sfs++;
//
//------enter vi in element list for vj
            ARRAY1D(mark, vj) += 1;
            ARRAY1D(v, sfs) = vi;
            ARRAY1D(l, sfs) = ARRAY1D(l, vj);
            ARRAY1D(l, vj) = sfs;
            sfs++;
        }
    }
//
//----create degree lists and initialize mark vector
    for (vi = 1; vi <= n; ++vi) {
        dvi = ARRAY1D(mark, vi);
        ARRAY1D(next, vi) = ARRAY1D(head, dvi);
        ARRAY1D(head, dvi) = vi;
        ARRAY1D(last, vi) = -dvi;
        nextvi = ARRAY1D(next,  vi);
        if (nextvi > 0) ARRAY1D(last, nextvi) = vi;
        ARRAY1D(mark, vi) = tag;
    }
//
    return;
// ** error-  insufficient storage
LABEL_101:
    flag = 9 * n + vi;
    return;
}


/**
 * @fn MDM -- purge inactive elements and do mass elimination
 */
void Odepack::MDM(
    int vk, int &tail, int *v, int *l, int *last, int *next, int *mark
)
{
    int tag, s, ls, vs, b, lb, vb, blp, blpmax;
    int *es = &vs;
//
//----initialize tag and list of uneliminated neighbors
    tag = ARRAY1D(mark, vk);
    tail = vk;
//
//----for each vertex/element vs/es in element list of vk
    ls = ARRAY1D(l, vk);
LABEL_1:
    s = ls;
    if (s == 0) goto LABEL_5;
    ls = ARRAY1D(l, s);
    vs = ARRAY1D(v, s);
    if (ARRAY1D(next, vs) < 0) goto LABEL_2;
//
//------if vs is uneliminated vertex, then tag and append to list of
//------uneliminated neighbors
    ARRAY1D(mark, vs) = tag;
    ARRAY1D(l, tail) = s;
    tail = s;
    goto LABEL_4;
//
//------if es is active element, then ...
//--------for each vertex vb in boundary list of element es
LABEL_2:
    lb = ARRAY1D(l, *es);
    blpmax = ARRAY1D(last, *es);
    for (blp = 1; blp <= blpmax; ++blp) {
        b = lb;
        lb = ARRAY1D(l, b);
        vb = ARRAY1D(v, b);
//
//----------if vb is untagged vertex, then tag and append to list of
//----------uneliminated neighbors
        if (ARRAY1D(mark, vb) >= tag) continue;
        ARRAY1D(mark, vb) = tag;
        ARRAY1D(l, tail) = b;
        tail = b;
    }
//
//--------mark es inactive
    ARRAY1D(mark, *es) = tag;
//
LABEL_4:
    goto LABEL_1;
//
//----terminate list of uneliminated neighbors
LABEL_5:
    ARRAY1D(l, tail) = 0;
//
    return;
}


/**
 * @fn MDP -- purge inactive elements and do mass elimination
 */
void Odepack::MDP(
    int &k, int ek, int &tail, int *v, int *l, int *head, int *last, int *next, int *mark
)
{
    int tag, free, li, vi, lvi, evi, s, ls, es, ilp, ilpmax, i;
//
//----initialize tag
    tag = ARRAY1D(mark, ek);
//
//----for each vertex vi in ek
    li = ek;
    ilpmax = ARRAY1D(last, ek);
    if (ilpmax <= 0) goto LABEL_12;
    for (ilp = 1; ilp <= ilpmax; ++ilp) {
        i = li;
        li = ARRAY1D(l, i);
        vi = ARRAY1D(v, li);
//
//------remove vi from degree list
        if (ARRAY1D(last, vi) == 0) goto LABEL_3;
        if (ARRAY1D(last, vi) > 0) goto LABEL_1;
        ARRAY1D(head, -ARRAY1D(last, vi)) = ARRAY1D(next, vi);
        goto LABEL_2;
LABEL_1:
        ARRAY1D(next, ARRAY1D(last, vi)) = ARRAY1D(next, vi);
LABEL_2:
        if (ARRAY1D(next, vi) > 0) ARRAY1D(last, ARRAY1D(next, vi)) = ARRAY1D(last, vi);
//
//------remove inactive items from element list of vi
LABEL_3:
        ls = vi;
LABEL_4:
        s = ls;
        ls = ARRAY1D(l, s);
        if (ls == 0) goto LABEL_6;
        es = ARRAY1D(v, ls);
        if (ARRAY1D(mark, es) < tag) goto LABEL_5;
        free = ls;
        ARRAY1D(l, s) = ARRAY1D(l, ls);
        ls = s;
LABEL_5:
        goto LABEL_4;
//
//------if vi is interior vertex, then remove from list and eliminate
LABEL_6:
        lvi = ARRAY1D(l, vi);
        if (lvi != 0) goto LABEL_7;
        ARRAY1D(l, i) = ARRAY1D(l, li);
        li = i;
//
        k++;
        ARRAY1D(next, vi) = -k;
        ARRAY1D(last, ek) -= 1;
        continue;
//
//------else ...
//--------classify vertex vi
LABEL_7:
        if (ARRAY1D(l, lvi) != 0) goto LABEL_9;
        evi = ARRAY1D(v, lvi);
        if (ARRAY1D(next, evi) >= 0) goto LABEL_9;
        if (ARRAY1D(mark, evi) < 0) goto LABEL_8;
//
//----------if vi is prototype vertex, then mark as such, initialize
//----------overlap count for corresponding element, and move vi to end
//----------of boundary list
        ARRAY1D(last, vi) = evi;
        ARRAY1D(mark, evi) = -1;
        ARRAY1D(l, tail) = li;
        tail = li;
        ARRAY1D(l, i) = ARRAY1D(l, li);
        li = i;
        goto LABEL_10;
//
//----------else if vi is duplicate vertex, then mark as such and adjust
//----------overlap count for corresponding element
LABEL_8:
        ARRAY1D(last, vi) = 0;
        ARRAY1D(mark, evi) -= 1;
        goto LABEL_10;
//
//----------else mark vi to compute degree
LABEL_9:
        ARRAY1D(last, vi) = -ek;
//
//--------insert ek in element list of vi
LABEL_10:
        ARRAY1D(v, free) = ek;
        ARRAY1D(l, free) = ARRAY1D(l, vi);
        ARRAY1D(l, vi) = free;
    }
//
//----terminate boundary list
LABEL_12:
    ARRAY1D(l, tail) = 0;
//
    return;
}


/**
 * @fn MDU -- update degrees of uneliminated vertices in ek
 */
void Odepack::MDU(
    int ek, int &dmin, int *v, int *l, int *head, int *last, int *next, int *mark
)
{
    int tag, vi, evi, dvi, s, vs, b, vb, ilp, ilpmax, blp, blpmax, i;
    int *es = &vs;
//
//----initialize tag
    tag = ARRAY1D(mark, ek) - ARRAY1D(last, ek);
//
//----for each vertex vi in ek
    i = ek;
    ilpmax = ARRAY1D(last, ek);
    if (ilpmax <= 0) goto LABEL_11;
    for (ilp = 1; ilp <= ilpmax; ++ilp) {
        i = ARRAY1D(l, i);
        vi = ARRAY1D(v, i);
        if (ARRAY1D(last, vi) < 0) {
            goto LABEL_1;
        } else if (ARRAY1D(last, vi) == 0) {
            continue;
        } else {
            goto LABEL_8;
        }
//
//------if vi neither prototype nor duplicate vertex, then merge elements
//------to compute degree
LABEL_1:
        tag++;
        dvi = ARRAY1D(last, ek);
//
//--------for each vertex/element vs/es in element list of vi
        s = ARRAY1D(l, vi);
LABEL_2:
        s = ARRAY1D(l, s);
        if (s == 0) goto LABEL_9;
        vs = ARRAY1D(v, s);
        if (ARRAY1D(next, vs) < 0) goto LABEL_3;
//
//----------if vs is uneliminated vertex, then tag and adjust degree
        ARRAY1D(mark, vs) = tag;
        dvi++;
        goto LABEL_5;
//
//----------if es is active element, then expand
//------------check for outmatched vertex
LABEL_3:
        if (ARRAY1D(mark, *es) < 0) goto LABEL_6;
//
//------------for each vertex vb in es
        b = *es;
        blpmax = ARRAY1D(last, *es);
        for (blp = 1; blp <= blpmax; ++blp) {
            b = ARRAY1D(l, b);
            vb = ARRAY1D(v, b);
//
//--------------if vb is untagged, then tag and adjust degree
            if (ARRAY1D(mark, vb) >= tag) continue;
            ARRAY1D(mark, vb) = tag;
            dvi++;
        }
//
LABEL_5:
        goto LABEL_2;
//
//------else if vi is outmatched vertex, then adjust overlaps but do not
//------compute degree
LABEL_6:
        ARRAY1D(last, vi) = 0;
        ARRAY1D(mark, *es) -= 1;
LABEL_7:
        s = ARRAY1D(l, s);
        if (s == 0) continue;
        *es = ARRAY1D(v, s);
        if (ARRAY1D(mark, *es) < 0) ARRAY1D(mark, *es) -= 1;
        goto LABEL_7;
//
//------else if vi is prototype vertex, then calculate degree by
//------inclusion/exclusion and reset overlap count
LABEL_8:
        evi = ARRAY1D(last, vi);
        dvi = ARRAY1D(last, ek) + ARRAY1D(last, evi) + ARRAY1D(mark, evi);
        ARRAY1D(mark, evi) = 0;
//
//------insert vi in appropriate degree list
LABEL_9:
        ARRAY1D(next, vi) = ARRAY1D(head, dvi);
        ARRAY1D(head, dvi) = vi;
        ARRAY1D(last, vi) = -dvi;
        if (ARRAY1D(next, vi) > 0) ARRAY1D(last, ARRAY1D(next, vi)) = vi;
        if (dvi < dmin) dmin = dvi;
    }
//
LABEL_11:
    return;
}


/**
 * @fn SRO -- -- symmetric reordering of sparse symmetric matrix
 * 
 * description
 * 
 *   the nonzero entries of the matrix m are assumed to be stored
 *   symmetrically in (ia,ja,a) format (i.e., not both m(i,j) and m(j,i)
 *   are stored if i ne j).
 * 
 *   sro does not rearrange the order of the rows, but does move
 *   nonzeroes from one row to another to ensure that if m(i,j) will be
 *   in the upper triangle of m with respect to the new ordering, then
 *   m(i,j) is stored in row i (and thus m(j,i) is not stored),  whereas
 *   if m(i,j) will be in the strict lower triangle of m, then m(j,i) is
 *   stored in row j (and thus m(i,j) is not stored).
 * 
 * 
 * additional parameters
 * 
 *   q     - integer one-dimensional work array.  dimension = n
 * 
 *   r     - integer one-dimensional work array.  dimension = number of
 *           nonzero entries in the upper triangle of m
 * 
 *   dflag - logical variable.  if dflag = .true., then store nonzero
 *           diagonal elements at the beginning of the row
 * 
 * -----------------------------------------------------------------------
 */
void Odepack::SRO(
    int n, int *ip, int *ia, int *ja, odepack_cpp_real *a, int *q, int *r, bool dflag
)
{
    int i, j, jmin, jmax, jdummy, k, ilast, jak;
    odepack_cpp_real ak;
//
//--phase 1 -- find row in which to store each nonzero
//----initialize count of nonzeroes to be stored in each row
    for (i = 1; i <= n; ++i) {
        ARRAY1D(q, i) = 0;
    }
//
//----for each nonzero element a(j)
    for (i = 1; i <= n; ++i) {
        jmin = ARRAY1D(ia, i);
        jmax = ARRAY1D(ia, i + 1) - 1;
        if (jmin > jmax) continue;
        for (j = jmin; j <= jmax; ++j) {
//
//--------find row (=r(j)) and column (=ja(j)) in which to store a(j)
            k = ARRAY1D(ja, j);
            if (ARRAY1D(ip, k) < ARRAY1D(ip, i)) ARRAY1D(ja, j) = i;
            if (ARRAY1D(ip, k) >= ARRAY1D(ip, i)) k = i;
            ARRAY1D(r, j) = k;
//
//--------.... and increment count of nonzeroes (=q(r(j)) in that row
            ARRAY1D(q, k) += 1;
        }
    }
//
//
//--phase 2 -- find new ia and permutation to apply to (ja,a)
//----determine pointers to delimit rows in permuted (ja,a)
    for (i = 1; i <= n; ++i) {
        ARRAY1D(ia, i + 1) = ARRAY1D(ia, i) + ARRAY1D(q, i);
        ARRAY1D(q, i) = ARRAY1D(ia, i + 1);
    }
//
//----determine where each (ja(j),a(j)) is stored in permuted (ja,a)
//----for each nonzero element (in reverse order)
    ilast = 0;
    jmin = ARRAY1D(ia, 1);
    jmax = ARRAY1D(ia, n + 1) - 1;
    j = jmax;
    for (jdummy = jmin; jdummy <= jmax; ++jdummy) {
        i = ARRAY1D(r, j);
        if (!dflag || ARRAY1D(ja, j) != i || i == ilast) goto LABEL_5;
//
//------if dflag, then put diagonal nonzero at beginning of row
        ARRAY1D(r, j) = ARRAY1D(ia, i);
        ilast = i;
        continue;
//
//------put (off-diagonal) nonzero in last unused location in row
LABEL_5:
        ARRAY1D(q, i) -= 1;
        ARRAY1D(r, j) = ARRAY1D(q, i);
//
        j--;
    }
//
//
//--phase 3 -- permute (ja,a) to upper triangular form (wrt new ordering)
    for (j = jmin; j <= jmax; ++j) {
LABEL_7:
        if (ARRAY1D(r, j) == j) continue;
        k = ARRAY1D(r, j);
        ARRAY1D(r, j) = ARRAY1D(r, k);
        ARRAY1D(r, k) = k;
        jak = ARRAY1D(ja, k);
        ARRAY1D(ja, k) = ARRAY1D(ja, j);
        ARRAY1D(ja, j) = jak;
        ak = ARRAY1D(a, k);
        ARRAY1D(a, k) = ARRAY1D(a, j);
        ARRAY1D(a, j) = ak;
        goto LABEL_7;
    }
//
    return;
}


/**
 * @fn CDRV
 * driver for subroutines for solving sparse nonsymmetric systems of
 * linear equations (compressed pointer storage)
 * 
 * parameters
 * class abbreviations are--
 *    n - integer variable
 *    f - real variable
 *    v - supplies a value to the driver
 *    r - returns a result from the driver
 *    i - used internally by the driver
 *    a - array
 * 
 * class - parameter
 * ------+----------
 *       -
 *        the nonzero entries of the coefficient matrix m are stored
 *   row-by-row in the array a.  to identify the individual nonzero
 *   entries in each row, we need to know in which column each entry
 *   lies.  the column indices which correspond to the nonzero entries
 *   of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
 *   ja(k) = j.  in addition, we need to know where each row starts and
 *   how long it is.  the index positions in ja and a where the rows of
 *   m begin are stored in the array ia.  i.e., if m(i,j) is the first
 *   nonzero entry (stored) in the i-th row and a(k) = m(i,j),  then
 *   ia(i) = k.  moreover, the index in ja and a of the first location
 *   following the last element in the last row is stored in ia(n+1).
 *   thus, the number of entries in the i-th row is given by
 *   ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
 *   consecutively in
 *           a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
 *   and the corresponding column indices are stored consecutively in
 *           ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
 *   for example, the 5 by 5 matrix
 *               ( 1. 0. 2. 0. 0.)
 *               ( 0. 3. 0. 0. 0.)
 *           m = ( 0. 4. 5. 6. 0.)
 *               ( 0. 0. 0. 7. 0.)
 *               ( 0. 0. 0. 8. 9.)
 *   would be stored as
 *              - 1  2  3  4  5  6  7  8  9
 *           ---+--------------------------
 *           ia - 1  3  4  7  8 10
 *           ja - 1  3  2  2  3  4  4  4  5
 *            a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
 * 
 * nv    - n     - number of variables/equations.
 * fva   - a     - nonzero entries of the coefficient matrix m, stored
 *       -           by rows.
 *       -           size = number of nonzero entries in m.
 * nva   - ia    - pointers to delimit the rows in a.
 *       -           size = n+1.
 * nva   - ja    - column numbers corresponding to the elements of a.
 *       -           size = size of a.
 * fva   - b     - right-hand side b.  b and z can the same array.
 *       -           size = n.
 * fra   - z     - solution x.  b and z can be the same array.
 *       -           size = n.
 * 
 *        the rows and columns of the original matrix m can be
 *   reordered (e.g., to reduce fillin or ensure numerical stability)
 *   before calling the driver.  if no reordering is done, then set
 *   r(i) = c(i) = ic(i) = i  for i=1,...,n.  the solution z is returned
 *   in the original order.
 *        if the columns have been reordered (i.e.,  c(i).ne.i  for some
 *   i), then the driver will call a subroutine (nroc) which rearranges
 *   each row of ja and a, leaving the rows in the original order, but
 *   placing the elements of each row in increasing order with respect
 *   to the new ordering.  if  path.ne.1,  then nroc is assumed to have
 *   been called already.
 * 
 * nva   - r     - ordering of the rows of m.
 *       -           size = n.
 * nva   - c     - ordering of the columns of m.
 *       -           size = n.
 * nva   - ic    - inverse of the ordering of the columns of m.  i.e.,
 *       -           ic(c(i)) = i  for i=1,...,n.
 *       -           size = n.
 * 
 *         the solution of the system of linear equations is divided into
 *    three stages --
 *      nsfc -- the matrix m is processed symbolically to determine where
 *              fillin will occur during the numeric factorization.
 *      nnfc -- the matrix m is factored numerically into the product ldu
 *              of a unit lower triangular matrix l, a diagonal matrix
 *              d, and a unit upper triangular matrix u, and the system
 *              mx = b  is solved.
 *     nnsc -- the linear system  mx = b  is solved using the ldu
 * or           factorization from nnfc.
 *     nntc -- the transposed linear system  mt x = b  is solved using
 *              the ldu factorization from nnf.
 *   for several systems whose coefficient matrices have the same
 *   nonzero structure, nsfc need be done only once (for the first
 *   system).  then nnfc is done once for each additional system.  for
 *   several systems with the same coefficient matrix, nsfc and nnfc
 *   need be done only once (for the first system).  then nnsc or nntc
 *   is done once for each additional right-hand side.
 * 
 * nv    - path  - path specification.  values and their meanings are --
 *       -           1  perform nroc, nsfc, and nnfc.
 *       -           2  perform nnfc only  (nsfc is assumed to have been
 *       -               done in a manner compatible with the storage
 *       -               allocation used in the driver).
 *       -           3  perform nnsc only  (nsfc and nnfc are assumed to
 *       -               have been done in a manner compatible with the
 *       -               storage allocation used in the driver).
 *       -           4  perform nntc only  (nsfc and nnfc are assumed to
 *       -               have been done in a manner compatible with the
 *       -               storage allocation used in the driver).
 *       -           5  perform nroc and nsfc.
 * 
 *         various errors are detected by the driver and the individual
 *    subroutines.
 * 
 * nr    - flag  - error flag.  values and their meanings are --
 *       -             0     no errors detected
 *       -             n+k   null row in a  --  row = k
 *       -            2n+k   duplicate entry in a  --  row = k
 *       -            3n+k   insufficient storage in nsfc  --  row = k
 *       -            4n+1   insufficient storage in nnfc
 *       -            5n+k   null pivot  --  row = k
 *       -            6n+k   insufficient storage in nsfc  --  row = k
 *       -            7n+1   insufficient storage in nnfc
 *       -            8n+k   zero pivot  --  row = k
 *       -           10n+1   insufficient storage in cdrv
 *       -           11n+1   illegal path specification
 * 
 *         working storage is needed for the factored form of the matrix
 *    m plus various temporary vectors.  the arrays isp and rsp should be
 *    equivalenced.  integer storage is allocated from the beginning of
 *    isp and real storage from the end of rsp.
 * 
 * nv    - nsp   - declared dimension of rsp.  nsp generally must
 *       -           be larger than  8n+2 + 2k  (where  k = (number of
 *       -           nonzero entries in m)).
 * nvira - isp   - integer working storage divided up into various arrays
 *       -           needed by the subroutines.  isp and rsp should be
 *       -           equivalenced.
 *       -           size = lratio*nsp.
 * fvira - rsp   - real working storage divided up into various arrays
 *       -           needed by the subroutines.  isp and rsp should be
 *       -           equivalenced.
 *       -           size = nsp.
 * nr    - esp   - if sufficient storage was available to perform the
 *       -           symbolic factorization (nsfc), then esp is set to
 *       -           the amount of excess storage provided (negative if
 *       -           insufficient storage was available to perform the
 *       -           numeric factorization (nnfc)).
 * 
 * 
 *  conversion to odepack_cpp_real precision
 * 
 *    to convert these routines for odepack_cpp_real precision arrays..
 *    (1) use the odepack_cpp_real precision declarations in place of the real
 *    declarations in each subprogram, as given in comment cards.
 *    (2) change the data-loaded value of the integer  lratio
 *    in subroutine cdrv, as indicated below.
 *    (3) change e0 to d0 in the constants in statement number 10
 *    in subroutine nnfc and the line following that.
 */
void Odepack::CDRV(
    int n, int *r, int *c, int *ic, int *ia, int *ja, odepack_cpp_real *a, odepack_cpp_real *b, odepack_cpp_real *z, int nsp, int *isp, odepack_cpp_real *rsp, int &esp, int path, int &flag
)
{
    int i, il, ijl, iu, iju, irl, jrl, jl, max, jlmax, ira, jra, irac, iru, jru, ju, jutmp, jumax, j, l, lmax;
    int d, u, q, row, tmp, ar, umax;
//
    int lratio = 2;
//
    if (path < 1 || 5 < path) goto LABEL_111;
//******initialize and divide up temporary storage  *******************
    il  = 1;
    ijl = il  + n + 1;
    iu  = ijl + n;
    iju = iu  + n + 1;
    irl = iju + n;
    jrl = irl + n;
    jl  = jrl + n;
//
// ******  reorder a if necessary, call nsfc if flag is set  ***********
    if ((path - 1) * (path - 5) != 0) goto LABEL_5;
    max   = lratio * nsp + 1 - jl - (n + 1) - 5 * n;
    jlmax = max / 2;
    q     = jl   + jlmax;
    ira   = q    + n + 1;
    jra   = ira  + n;
    irac  = jra  + n;
    iru   = irac + n;
    jru   = iru  + n;
    jutmp = jru  + n;
    jumax = lratio * nsp + 1 - jutmp;
    esp   = max / lratio;
    if (jlmax <= 0 || jumax <= 0) goto LABEL_110;
//
    for (i = 1; i <= n; ++i) {
        if (ARRAY1D(c, i) != i) goto LABEL_2;
    }
    goto LABEL_3;
LABEL_2:
    ar = nsp + 1 - n;
    NROC(n, ic, ia, ja, a, &ARRAY1D(isp, il), &ARRAY1D(rsp, ar), &ARRAY1D(isp, iu), flag);
    if (flag != 0) goto LABEL_100;
//
LABEL_3:
    NSFC(n, r, ic, ia, ja, 
        jlmax, &ARRAY1D(isp, il), &ARRAY1D(isp, jl), &ARRAY1D(isp, ijl),
        jumax, &ARRAY1D(isp, iu), &ARRAY1D(isp, jutmp), &ARRAY1D(isp, iju),
        &ARRAY1D(isp, q), &ARRAY1D(isp, ira), &ARRAY1D(isp, jra), &ARRAY1D(isp, irac),
        &ARRAY1D(isp, irl), &ARRAY1D(isp, jrl), &ARRAY1D(isp, iru), &ARRAY1D(isp, jru), flag);
    if (flag != 0) goto LABEL_100;
//  ******  move ju next to jl  *****************************************
    jlmax = ARRAY1D(isp, ijl + n - 1);
    ju    = jl + jlmax;
    jumax = ARRAY1D(isp, iju + n - 1);
    if (jumax <= 0) goto LABEL_5;
    for (j = 1; j <= jumax; ++j) {
        ARRAY1D(isp, ju + j - 1) = ARRAY1D(isp, jutmp + j - 1);
    }
//
//  ******  call remaining subroutines  *********************************
LABEL_5:
    jlmax = ARRAY1D(isp, ijl + n - 1);
    ju    = jl + jlmax;
    jumax = ARRAY1D(isp, iju + n - 1);
    l     = (ju + jumax - 2 + lratio) / lratio + 1;
    lmax  = ARRAY1D(isp, il + n) - 1;
    d     = l + lmax;
    u     = d + n;
    row   = nsp + 1 - n;
    tmp   = row - n;
    umax  = tmp - u;
    esp   = umax - (ARRAY1D(isp, iu + n) - 1);
//
    if ((path - 1) * (path - 2) != 0) goto LABEL_6;
    if (umax < 0) goto LABEL_110;
    NNFC(n, r, c, ic, ia, ja, a, z, b, 
        lmax, &ARRAY1D(isp, il), &ARRAY1D(isp, jl), &ARRAY1D(isp, ijl), &ARRAY1D(rsp, l), &ARRAY1D(rsp, d),
        umax, &ARRAY1D(isp, iu), &ARRAY1D(isp, ju), &ARRAY1D(isp, iju), &ARRAY1D(rsp, u),
        &ARRAY1D(rsp, row), &ARRAY1D(rsp, tmp), &ARRAY1D(isp, irl), &ARRAY1D(isp, jrl), flag);
    if (flag != 0) goto LABEL_100;
//
LABEL_6:
    if ((path - 3) != 0) goto LABEL_7;
    NNSC(n, r, c, &ARRAY1D(isp, il), &ARRAY1D(isp, jl), &ARRAY1D(isp, ijl), &ARRAY1D(rsp, l),
        &ARRAY1D(rsp, d), &ARRAY1D(isp, iu), &ARRAY1D(isp, ju), &ARRAY1D(isp, iju), &ARRAY1D(rsp, u),
        z, b, &ARRAY1D(rsp, tmp));
//
LABEL_7:
    if ((path - 4) != 0) goto LABEL_8;
    NNTC(n, r, c, &ARRAY1D(isp, il), &ARRAY1D(isp, jl), &ARRAY1D(isp, ijl), &ARRAY1D(rsp, l),
        &ARRAY1D(rsp, d), &ARRAY1D(isp, iu), &ARRAY1D(isp, ju), &ARRAY1D(isp, iju), &ARRAY1D(rsp, u),
        z, b, &ARRAY1D(rsp, tmp));
LABEL_8:
    return;
//
// ** error.. error detected in nroc, nsfc, nnfc, or nnsc
LABEL_100:
    return;
// ** error.. insufficient storage
LABEL_110:
    flag = 10*n + 1;
    return;
// ** error.. illegal path specification
LABEL_111:
    flag = 11*n + 1;
    return;
}


/**
 * @fn NROC
c       ----------------------------------------------------------------
c
c               yale sparse matrix package - nonsymmetric codes
c                    solving the system of equations mx = b
c
c    i.   calling sequences
c         the coefficient matrix can be processed by an ordering routine
c    (e.g., to reduce fillin or ensure numerical stability) before using
c    the remaining subroutines.  if no reordering is done, then set
c    r(i) = c(i) = ic(i) = i  for i=1,...,n.  if an ordering subroutine
c    is used, then nroc should be used to reorder the coefficient matrix
c    the calling sequence is --
c        (       (matrix ordering))
c        (nroc   (matrix reordering))
c         nsfc   (symbolic factorization to determine where fillin will
c                  occur during numeric factorization)
c         nnfc   (numeric factorization into product ldu of unit lower
c                  triangular matrix l, diagonal matrix d, and unit
c                  upper triangular matrix u, and solution of linear
c                  system)
c         nnsc   (solution of linear system for additional right-hand
c                  side using ldu factorization from nnfc)
c    (if only one system of equations is to be solved, then the
c    subroutine trk should be used.)
c
c    ii.  storage of sparse matrices
c         the nonzero entries of the coefficient matrix m are stored
c    row-by-row in the array a.  to identify the individual nonzero
c    entries in each row, we need to know in which column each entry
c    lies.  the column indices which correspond to the nonzero entries
c    of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
c    ja(k) = j.  in addition, we need to know where each row starts and
c    how long it is.  the index positions in ja and a where the rows of
c    m begin are stored in the array ia.  i.e., if m(i,j) is the first
c    (leftmost) entry in the i-th row and  a(k) = m(i,j),  then
c    ia(i) = k.  moreover, the index in ja and a of the first location
c    following the last element in the last row is stored in ia(n+1).
c    thus, the number of entries in the i-th row is given by
c    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
c    consecutively in
c            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
c    and the corresponding column indices are stored consecutively in
c            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
c    for example, the 5 by 5 matrix
c                ( 1. 0. 2. 0. 0.)
c                ( 0. 3. 0. 0. 0.)
c            m = ( 0. 4. 5. 6. 0.)
c                ( 0. 0. 0. 7. 0.)
c                ( 0. 0. 0. 8. 9.)
c    would be stored as
c               - 1  2  3  4  5  6  7  8  9
c            ---+--------------------------
c            ia - 1  3  4  7  8 10
c            ja - 1  3  2  2  3  4  4  4  5
c             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
c
c         the strict upper (lower) triangular portion of the matrix
c    u (l) is stored in a similar fashion using the arrays  iu, ju, u
c    (il, jl, l)  except that an additional array iju (ijl) is used to
c    compress storage of ju (jl) by allowing some sequences of column
c    (row) indices to used for more than one row (column)  (n.b., l is
c    stored by columns).  iju(k) (ijl(k)) points to the starting
c    location in ju (jl) of entries for the kth row (column).
c    compression in ju (jl) occurs in two ways.  first, if a row
c    (column) i was merged into the current row (column) k, and the
c    number of elements merged in from (the tail portion of) row
c    (column) i is the same as the final length of row (column) k, then
c    the kth row (column) and the tail of row (column) i are identical
c    and iju(k) (ijl(k)) points to the start of the tail.  second, if
c    some tail portion of the (k-1)st row (column) is identical to the
c    head of the kth row (column), then iju(k) (ijl(k)) points to the
c    start of that tail portion.  for example, the nonzero structure of
c    the strict upper triangular part of the matrix
c            d 0 x x x
c            0 d 0 x x
c            0 0 d x 0
c            0 0 0 d x
c            0 0 0 0 d
c    would be represented as
c                - 1 2 3 4 5 6
c            ----+------------
c             iu - 1 4 6 7 8 8
c             ju - 3 4 5 4
c            iju - 1 2 4 3           .
c    the diagonal entries of l and u are assumed to be equal to one and
c    are not stored.  the array d contains the reciprocals of the
c    diagonal entries of the matrix d.
c
c    iii. additional storage savings
c         in nsfc, r and ic can be the same array in the calling
c    sequence if no reordering of the coefficient matrix has been done.
c         in nnfc, r, c, and ic can all be the same array if no
c    reordering has been done.  if only the rows have been reordered,
c    then c and ic can be the same array.  if the row and column
c    orderings are the same, then r and c can be the same array.  z and
c    row can be the same array.
c         in nnsc or nntc, r and c can be the same array if no
c    reordering has been done or if the row and column orderings are the
c    same.  z and b can be the same array.  however, then b will be
c    destroyed.
c
c    iv.  parameters
c         following is a list of parameters to the programs.  names are
c    uniform among the various subroutines.  class abbreviations are --
c       n - integer variable
c       f - real variable
c       v - supplies a value to a subroutine
c       r - returns a result from a subroutine
c       i - used internally by a subroutine
c       a - array
c
c class - parameter
c ------+----------
c fva   - a     - nonzero entries of the coefficient matrix m, stored
c       -           by rows.
c       -           size = number of nonzero entries in m.
c fva   - b     - right-hand side b.
c       -           size = n.
c nva   - c     - ordering of the columns of m.
c       -           size = n.
c fvra  - d     - reciprocals of the diagonal entries of the matrix d.
c       -           size = n.
c nr    - flag  - error flag.  values and their meanings are --
c       -            0     no errors detected
c       -            n+k   null row in a  --  row = k
c       -           2n+k   duplicate entry in a  --  row = k
c       -           3n+k   insufficient storage for jl  --  row = k
c       -           4n+1   insufficient storage for l
c       -           5n+k   null pivot  --  row = k
c       -           6n+k   insufficient storage for ju  --  row = k
c       -           7n+1   insufficient storage for u
c       -           8n+k   zero pivot  --  row = k
c nva   - ia    - pointers to delimit the rows of a.
c       -           size = n+1.
c nvra  - ijl   - pointers to the first element in each column in jl,
c       -           used to compress storage in jl.
c       -           size = n.
c nvra  - iju   - pointers to the first element in each row in ju, used
c       -           to compress storage in ju.
c       -           size = n.
c nvra  - il    - pointers to delimit the columns of l.
c       -           size = n+1.
c nvra  - iu    - pointers to delimit the rows of u.
c       -           size = n+1.
c nva   - ja    - column numbers corresponding to the elements of a.
c       -           size = size of a.
c nvra  - jl    - row numbers corresponding to the elements of l.
c       -           size = jlmax.
c nv    - jlmax - declared dimension of jl.  jlmax must be larger than
c       -           the number of nonzeros in the strict lower triangle
c       -           of m plus fillin minus compression.
c nvra  - ju    - column numbers corresponding to the elements of u.
c       -           size = jumax.
c nv    - jumax - declared dimension of ju.  jumax must be larger than
c       -           the number of nonzeros in the strict upper triangle
c       -           of m plus fillin minus compression.
c fvra  - l     - nonzero entries in the strict lower triangular portion
c       -           of the matrix l, stored by columns.
c       -           size = lmax.
c nv    - lmax  - declared dimension of l.  lmax must be larger than
c       -           the number of nonzeros in the strict lower triangle
c       -           of m plus fillin  (il(n+1)-1 after nsfc).
c nv    - n     - number of variables/equations.
c nva   - r     - ordering of the rows of m.
c       -           size = n.
c fvra  - u     - nonzero entries in the strict upper triangular portion
c       -           of the matrix u, stored by rows.
c       -           size = umax.
c nv    - umax  - declared dimension of u.  umax must be larger than
c       -           the number of nonzeros in the strict upper triangle
c       -           of m plus fillin  (iu(n+1)-1 after nsfc).
c fra   - z     - solution x.
c       -           size = n.
c
c       ----------------------------------------------------------------
c
c*** subroutine nroc
c*** reorders rows of a, leaving row order unchanged
c
c
c       input parameters.. n, ic, ia, ja, a
c       output parameters.. ja, a, flag
c
c       parameters used internally..
c nia   - p     - at the kth step, p is a linked list of the reordered
c       -           column indices of the kth row of a.  p(n+1) points
c       -           to the first entry in the list.
c       -           size = n+1.
c nia   - jar   - at the kth step,jar contains the elements of the
c       -           reordered column indices of a.
c       -           size = n.
c fia   - ar    - at the kth step, ar contains the elements of the
c       -           reordered row of a.
c       -           size = n.
 */
void Odepack::NROC(int n, int *ic, int *ia, int *ja, odepack_cpp_real *a, int *jar, 
    odepack_cpp_real *ar, int *p, int &flag)
{
    int i, j, jmin, jmax, newj, k;
//
//  ******  for each nonempty row  *******************************
    for (k = 1; k <= n; ++k) {
        jmin = ARRAY1D(ia, k);
        jmax = ARRAY1D(ia, k + 1) - 1;
        if (jmin > jmax) continue;
        ARRAY1D(p, n + 1) = n + 1;
//  ******  insert each element in the list  *********************
        for (j = jmin; j <= jmax; ++j) {
            newj = ARRAY1D(ic, ARRAY1D(ja, j));
            i = n + 1;
LABEL_1:
            if (ARRAY1D(p, i) >= newj) goto LABEL_2;
            i = ARRAY1D(p, i);
            goto LABEL_1;
LABEL_2:
            if (ARRAY1D(p, i) == newj) goto LABEL_102;
            ARRAY1D(p, newj) = ARRAY1D( p, i);
            ARRAY1D(p, i) = newj;
            ARRAY1D(jar, newj) = ARRAY1D(ja, j);
            ARRAY1D(ar, newj) = ARRAY1D(a, j);
        }
//  ******  replace old row in ja and a  *************************
        i = n + 1;
        for (j = jmin; j <= jmax; ++j) {
            i = ARRAY1D(p, i);
            ARRAY1D(ja, j) = ARRAY1D(jar, i);
            ARRAY1D(a, j) = ARRAY1D(ar, i);
        }
    }
    flag = 0;
    return;
//
// ** error.. duplicate entry in a
LABEL_102:
    flag = n + k;
    return;
}


/**
 * @fn NSFC
 * symbolic ldu-factorization of nonsymmetric sparse matrix
 * (compressed pointer storage)
 * 
 * 
c       input variables.. n, r, ic, ia, ja, jlmax, jumax.
c       output variables.. il, jl, ijl, iu, ju, iju, flag.
c
c       parameters used internally..
c nia   - q     - suppose  m*  is the result of reordering  m.  if
c       -           processing of the ith row of  m*  (hence the ith
c       -           row of  u) is being done,  q(j)  is initially
c       -           nonzero if  m*(i,j) is nonzero (j.ge.i).  since
c       -           values need not be stored, each entry points to the
c       -           next nonzero and  q(n+1)  points to the first.  n+1
c       -           indicates the end of the list.  for example, if n=9
c       -           and the 5th row of  m*  is
c       -              0 x x 0 x 0 0 x 0
c       -           then  q  will initially be
c       -              a a a a 8 a a 10 5           (a - arbitrary).
c       -           as the algorithm proceeds, other elements of  q
c       -           are inserted in the list because of fillin.
c       -           q  is used in an analogous manner to compute the
c       -           ith column of  l.
c       -           size = n+1.
c nia   - ira,  - vectors used to find the columns of  m.  at the kth
c nia   - jra,      step of the factorization,  irac(k)  points to the
c nia   - irac      head of a linked list in  jra  of row indices i
c       -           such that i .ge. k and  m(i,k)  is nonzero.  zero
c       -           indicates the end of the list.  ira(i)  (i.ge.k)
c       -           points to the smallest j such that j .ge. k and
c       -           m(i,j)  is nonzero.
c       -           size of each = n.
c nia   - irl,  - vectors used to find the rows of  l.  at the kth step
c nia   - jrl       of the factorization,  jrl(k)  points to the head
c       -           of a linked list in  jrl  of column indices j
c       -           such j .lt. k and  l(k,j)  is nonzero.  zero
c       -           indicates the end of the list.  irl(j)  (j.lt.k)
c       -           points to the smallest i such that i .ge. k and
c       -           l(i,j)  is nonzero.
c       -           size of each = n.
c nia   - iru,  - vectors used in a manner analogous to  irl and jrl
c nia   - jru       to find the columns of  u.
c       -           size of each = n.
c
c  internal variables..
c    jlptr - points to the last position used in  jl.
c    juptr - points to the last position used in  ju.
c    jmin,jmax - are the indices in  a or u  of the first and last
c                elements to be examined in a given row.
c                for example,  jmin=ia(k), jmax=ia(k+1)-1.
 */
void Odepack::NSFC(
    int n, int *r, int *ic, int *ia, int *ja, int jlmax, int *il, int *jl, int *ijl, int jumax, int *iu, int *ju, int *iju, int *q,
    int *ira, int *jra, int *irac, int *irl, int *jrl, int *iru, int *jru, int &flag
)
{
    int i, i1, np1, j, jlmin, jlptr, jumin, juptr, k, rk, iak, jaiak, luk, vj, qm, m, lastid, lasti,
        jmin, jmax, llong, jtmp, irll, cend, irul, rend, irai, jairai;
//
//  ******  initialize pointers  ****************************************
    np1 = n + 1;
    jlmin = 1;
    jlptr = 0;
    ARRAY1D(il, 1) = 1;
    jumin = 1;
    juptr = 0;
    ARRAY1D(iu, 1) = 1;
    for (k = 1; k <= n; ++k) {
        ARRAY1D(irac, k) = 0;
        ARRAY1D(jra, k) = 0;
        ARRAY1D(jrl, k) = 0;
        ARRAY1D(jru, k) = 0;
    }
//  ******  initialize column pointers for a  ***************************
    for (k = 1; k <= n; ++k) {
        rk = ARRAY1D(r, k);
        iak = ARRAY1D(ia, rk);
        if (iak >= ARRAY1D(ia, rk + 1)) goto LABEL_101;
        jaiak = ARRAY1D(ic, ARRAY1D(ja, iak));
        if (jaiak > k) goto LABEL_105;
        ARRAY1D(jra, k) = ARRAY1D(irac, jaiak);
        ARRAY1D(irac, jaiak) = k;
        ARRAY1D(ira, k) = iak;
    }
//
//  ******  for each column of l and row of u  **************************
    for (k = 1; k <= n; ++k) {
//
//  ******  initialize q for computing kth column of l  *****************
        ARRAY1D(q, np1) = np1;
        luk = -1;
//  ******  by filling in kth column of a  ******************************
        vj = ARRAY1D(irac, k);
        if (vj == 0) goto LABEL_5;
LABEL_3:
        qm = np1;
LABEL_4:
        m = qm;
        qm = ARRAY1D(q, m);
        if (qm  < vj) goto LABEL_4;
        if (qm == vj) goto LABEL_102;
        luk++;
        ARRAY1D(q, m) = vj;
        ARRAY1D(q, vj) = qm;
        vj = ARRAY1D(jra, vj);
        if (vj != 0) goto LABEL_3;
//  ******  link through jru  *******************************************
LABEL_5:
        lastid = 0;
        lasti = 0;
        ARRAY1D(ijl, k) = jlptr;
        i = k;
LABEL_6:
        i = ARRAY1D(jru, i);
        if (i == 0) goto LABEL_10;
        qm = np1;
        jmin = ARRAY1D(irl, i);
        jmax = ARRAY1D(ijl, i) + ARRAY1D(il, i + 1) - ARRAY1D(il, i) - 1;
        llong = jmax - jmin; // llong = long
        if (llong < 0) goto LABEL_6;
        jtmp  = ARRAY1D(jl, jmin);
        if (jtmp != k) llong++;
        if (jtmp == k) ARRAY1D(r, i) = - ARRAY1D(r, i);
        if (lastid >= llong) goto LABEL_7;
        lasti = i;
        lastid = llong;
//  ******  and merge the corresponding columns into the kth column  ****
LABEL_7:
        for (j = jmin; j <= jmax; ++j) {
            vj = ARRAY1D(jl, j);
LABEL_8:
            m = qm;
            qm = ARRAY1D(q, m);
            if (qm  < vj) goto LABEL_8;
            if (qm == vj) continue;
            luk++;
            ARRAY1D(q, m) = vj;
            ARRAY1D(q, vj) = qm;
            qm = vj;
        }
        goto LABEL_6;
//  ******  lasti is the longest column merged into the kth  ************
//  ******  see if it equals the entire kth column  *********************
LABEL_10:
        qm = ARRAY1D(q, np1);
        if (qm  != k) goto LABEL_105;
        if (luk == 0) goto LABEL_17;
        if (lastid != luk) goto LABEL_11;
//  ******  if so, jl can be compressed  ********************************
        irll = ARRAY1D(irl, lasti);
        ARRAY1D(ijl, k) = irll + 1;
        if (ARRAY1D(jl, irll) != k) ARRAY1D(ijl, k) -= 1;
        goto LABEL_17;
//  ******  if not, see if kth column can overlap the previous one  *****
LABEL_11:
        if (jlmin > jlptr) goto LABEL_15;
        qm = ARRAY1D(q, qm);
        for (j = jlmin; j <= jlptr; ++j) {
            if (ARRAY1D(jl, j) < qm) {
                continue;
            } else if (ARRAY1D(jl, j) == qm) {
                goto LABEL_13;
            } else {
                goto LABEL_15;
            }
        }
        goto LABEL_15;
LABEL_13:
        ARRAY1D(ijl, k) = j;
        for (i = j; i <= jlptr; ++i) {
            if (ARRAY1D(jl, i) != qm) goto LABEL_15;
            qm = ARRAY1D(q, qm);
            if (qm > n) goto LABEL_17;
        }
        jlptr = j - 1;
//  ******  move column indices from q to jl, update vectors  ***********
LABEL_15:
        jlmin = jlptr + 1;
        ARRAY1D(ijl, k) = jlmin;
        if (luk == 0) goto LABEL_17;
        jlptr += luk;
        if (jlptr > jlmax) goto LABEL_103;
        qm = ARRAY1D(q, np1);
        for (j = jlmin; j <= jlptr; ++j) {
            qm = ARRAY1D(q, qm);
            ARRAY1D(jl, j) = qm;
        }
LABEL_17:
        ARRAY1D(irl, k) = ARRAY1D(ijl, k);
        ARRAY1D(il, k + 1) = ARRAY1D(il, k) + luk;
//
//  ******  initialize q for computing kth row of u  ********************
        ARRAY1D(q, np1) = np1;
        luk = -1;
//  ******  by filling in kth row of reordered a  ***********************
        rk = ARRAY1D(r, k);
        jmin = ARRAY1D(ira, k);
        jmax = ARRAY1D(ia, rk + 1) - 1;
        if (jmin > jmax) goto LABEL_20;
        for (j = jmin; j <= jmax; ++j) {
            vj = ARRAY1D(ic, ARRAY1D(ja, j));
            qm = np1;
LABEL_18:
            m = qm;
            qm = ARRAY1D(q, m);
            if (qm  < vj) goto LABEL_18;
            if (qm == vj) goto LABEL_102;
            luk++;
            ARRAY1D(q, m) = vj;
            ARRAY1D(q, vj) = qm;
        }
//  ******  link through jrl,  ******************************************
LABEL_20:
        lastid = 0;
        lasti = 0;
        ARRAY1D(iju, k) = juptr;
        i = k;
        i1 = ARRAY1D(jrl, k);
LABEL_21:
        i = i1;
        if (i == 0) goto LABEL_26;
        i1 = ARRAY1D(jrl, i);
        qm = np1;
        jmin = ARRAY1D(iru, i);
        jmax = ARRAY1D(iju, i) + ARRAY1D(iu, i + 1) - ARRAY1D(iu, i) - 1;
        llong = jmax - jmin;
        if (llong < 0) goto LABEL_21;
        jtmp = ARRAY1D(ju, jmin);
        if (jtmp == k) goto LABEL_22;
//  ******  update irl and jrl, *****************************************
        llong++;
        cend = ARRAY1D(ijl, i) + ARRAY1D(il, i + 1) - ARRAY1D(il, i);
        ARRAY1D(irl, i) += 1;
        if (ARRAY1D(irl, i) >= cend) goto LABEL_22;
        j = ARRAY1D(jl, ARRAY1D(irl, i));
        ARRAY1D(jrl, i) = ARRAY1D(jrl, j);
        ARRAY1D(jrl, j) = i;
LABEL_22:
        if (lastid >= llong) goto LABEL_23;
        lasti = i;
        lastid = llong;
//  ******  and merge the corresponding rows into the kth row  **********
LABEL_23:
        for (j = jmin; j <= jmax; ++j) {
            vj = ARRAY1D(ju, j);
LABEL_24:
            m = qm;
            qm = ARRAY1D(q, m);
            if (qm  < vj) goto LABEL_24;
            if (qm == vj) continue;;
            luk++;
            ARRAY1D(q, m) = vj;
            ARRAY1D(q, vj) = qm;
            qm = vj;
        }
        goto LABEL_21;
//  ******  update jrl(k) and irl(k)  ***********************************
LABEL_26:
        if (ARRAY1D(il, k + 1) <= ARRAY1D(il, k)) goto LABEL_27;
        j = ARRAY1D(jl, ARRAY1D(irl, k));
        ARRAY1D(jrl, k) = ARRAY1D(jrl, j);
        ARRAY1D(jrl, j) = k;
//  ******  lasti is the longest row merged into the kth  ***************
//  ******  see if it equals the entire kth row  ************************
LABEL_27:
        qm = ARRAY1D(q, np1);
        if (qm != k) goto LABEL_105;
        if (luk == 0) goto LABEL_34;
        if (lastid != luk) goto LABEL_28;
//  ******  if so, ju can be compressed  ********************************
        irul = ARRAY1D(iru, lasti);
        ARRAY1D(iju, k) = irul + 1;
        if (ARRAY1D(ju, irul) != k) ARRAY1D(iju, k) -= 1;
        goto LABEL_34;
//  ******  if not, see if kth row can overlap the previous one  ********
LABEL_28:
        if (jumin > juptr) goto LABEL_32;
        qm = ARRAY1D(q, qm);
        for (j = jumin; j <= juptr; ++j) {
            if (ARRAY1D(ju, j) < qm) {
                continue;
            } else if (ARRAY1D(ju, j) == qm) {
                goto LABEL_30;
            } else {
                goto LABEL_32;
            }
        }
        goto LABEL_32;
LABEL_30:
        ARRAY1D(iju, k) = j;
        for (i = j; i <= juptr; ++i) {
            if (ARRAY1D(ju, i) != qm) goto LABEL_32;
            qm = ARRAY1D(q, qm);
            if (qm > n) goto LABEL_34;
        }
        juptr = j - 1;
//  ******  move row indices from q to ju, update vectors  **************
LABEL_32:
        jumin = juptr + 1;
        ARRAY1D(iju, k) = jumin;
        if (luk == 0) goto LABEL_34;
        juptr += luk;
        if (juptr > jumax) goto LABEL_106;
        qm = ARRAY1D(q, np1);
        for (j = jumin; j <= juptr; ++j) {
            qm = ARRAY1D(q, qm);
            ARRAY1D(ju, j) = qm;
        }
LABEL_34:
        ARRAY1D(iru, k) = ARRAY1D(iju, k);
        ARRAY1D(iu, k + 1) = ARRAY1D(iu, k) + luk;
//
//  ******  update iru, jru  ********************************************
        i = k;
LABEL_35:
        i1 = ARRAY1D(jru, i);
        if (ARRAY1D(r, i) < 0) goto LABEL_36;
        rend = ARRAY1D(iju, i) + ARRAY1D(iu, i + 1) - ARRAY1D(iu, i);
        if (ARRAY1D(iru, i) >= rend) goto LABEL_37;
        j = ARRAY1D(ju, ARRAY1D(iru, i));
        ARRAY1D(jru, i) = ARRAY1D(jru, j);
        ARRAY1D(jru, j) = i;
        goto LABEL_37;
LABEL_36:
        ARRAY1D(r, i) = - ARRAY1D(r, i);
LABEL_37:
        i = i1;
        if (i == 0) goto LABEL_38;
        ARRAY1D(iru, i) += 1;
        goto LABEL_35;
//
//  ******  update ira, jra, irac  **************************************
LABEL_38:
        i = ARRAY1D(irac, k);
        if (i == 0) continue;
LABEL_39:
        i1 = ARRAY1D(jra, i);
        ARRAY1D(ira, i) += 1;
        if (ARRAY1D(ira, i) >= ARRAY1D(ia, ARRAY1D(r, i) + 1)) goto LABEL_40;
        irai = ARRAY1D(ira, i);
        jairai = ARRAY1D( ic, ARRAY1D(ja, irai));
        if (jairai > i) goto LABEL_40;
        ARRAY1D(jra, i) = ARRAY1D(irac, jairai);
        ARRAY1D(irac, jairai) = i;
LABEL_40:
        i = i1;
        if (i != 0) goto LABEL_39;
    }
//
    ARRAY1D(ijl, n) = jlptr;
    ARRAY1D(iju, n) = juptr;
    flag = 0;
    return;
//
// ** error.. null row in a
LABEL_101:
    flag = n + rk;
    return;
// ** error.. duplicate entry in a
LABEL_102:
    flag = 2*n + rk;
    return;
// ** error.. insufficient storage for jl
LABEL_103:
    flag = 3*n + k;
    return;
// ** error.. null pivot
LABEL_105:
    flag = 5*n + k;
    return;
// ** error.. insufficient storage for ju
LABEL_106:
    flag = 6*n + k;
    return;
}


/**
 * @fn NNCF
 * 
 * numerical ldu-factorization of sparse nonsymmetric matrix and
 * solution of system of linear equations (compressed pointer
 * storage)
 * 
 * 
c       input variables..  n, r, c, ic, ia, ja, a, b,
c                          il, jl, ijl, lmax, iu, ju, iju, umax
c       output variables.. z, l, d, u, flag
c
c       parameters used internally..
c nia   - irl,  - vectors used to find the rows of  l.  at the kth step
c nia   - jrl       of the factorization,  jrl(k)  points to the head
c       -           of a linked list in  jrl  of column indices j
c       -           such j .lt. k and  l(k,j)  is nonzero.  zero
c       -           indicates the end of the list.  irl(j)  (j.lt.k)
c       -           points to the smallest i such that i .ge. k and
c       -           l(i,j)  is nonzero.
c       -           size of each = n.
c fia   - row   - holds intermediate values in calculation of  u and l.
c       -           size = n.
c fia   - tmp   - holds new right-hand side  b*  for solution of the
c       -           equation ux = b*.
c       -           size = n.
c
c  internal variables..
c    jmin, jmax - indices of the first and last positions in a row to
c      be examined.
c    sum - used in calculating  tmp.
 */
void Odepack::NNFC(
    int n, int *r, int *c, int *ic, int *ia, int *ja, odepack_cpp_real *a, odepack_cpp_real *z, odepack_cpp_real *b,
    int lmax, int *il, int *jl, int *ijl, odepack_cpp_real *l, odepack_cpp_real *d, int umax, int *iu, int *ju, int *iju, odepack_cpp_real *u,
    odepack_cpp_real *row, odepack_cpp_real *tmp, int *irl, int *jrl, int &flag
)
{
    int i, i1, i2, k, j, jmin, jmax, rk, mu, ijlb;
    odepack_cpp_real lki, sum, dk;
//
//  ******  initialize pointers and test storage  ***********************
    if (ARRAY1D(il, n + 1) - 1 > lmax) goto LABEL_104;
    if (ARRAY1D(iu, n + 1) - 1 > umax) goto LABEL_107;
    for (k = 1; k <= n; ++k) {
        ARRAY1D(irl, k) = ARRAY1D(il, k);
        ARRAY1D(jrl, k) = 0;
    }
//
//  ******  for each row  ***********************************************
    for (k = 1; k <= n; ++k) {
//  ******  reverse jrl and zero row where kth row of l will fill in  ***
        ARRAY1D(row, k) = 0;
        i1 = 0;
        if (ARRAY1D(jrl, k) == 0) goto LABEL_3;
        i = ARRAY1D(jrl, k);
LABEL_2:
        i2 = ARRAY1D(jrl, i);
        ARRAY1D(jrl, i) = i1;
        i1 = i;
        ARRAY1D(row, i) = 0;
        i = i2;
        if (i != 0) goto LABEL_2;
//  ******  set row to zero where u will fill in  ***********************
LABEL_3:
        jmin = ARRAY1D(iju, k);
        jmax = jmin + ARRAY1D(iu, k + 1) - ARRAY1D(iu, k) - 1;
        if (jmin > jmax) goto LABEL_5;
        for (j = jmin; j <= jmax; ++j) {
            ARRAY1D(row, ARRAY1D(ju, j)) = 0;
        }
//  ******  place kth row of a in row  **********************************
LABEL_5:
        rk = ARRAY1D(r, k);
        jmin = ARRAY1D(ia, rk);
        jmax = ARRAY1D(ia, rk + 1) - 1;
        for (j = jmin; j <= jmax; ++j) {
            ARRAY1D(row, ARRAY1D(ic, ARRAY1D(ja, j))) = ARRAY1D(a, j);
        }
//  ******  initialize sum, and link through jrl  ***********************
        sum = ARRAY1D(b, rk);
        i = i1;
        if (i == 0) goto LABEL_10;
//  ******  assign the kth row of l and adjust row, sum  ****************
LABEL_7:
        lki = - ARRAY1D(row, i);
//  ******  if l is not required, then comment out the following line  **
        ARRAY1D(l, ARRAY1D(irl, i)) = -lki;
        sum = sum + lki * ARRAY1D(tmp, i);
        jmin = ARRAY1D(iu, i);
        jmax = ARRAY1D(iu, i + 1) - 1;
        if (jmin > jmax) goto LABEL_9;
        mu = ARRAY1D(iju, i) - jmin;
        for (j = jmin; j <= jmax; ++j) {
            ARRAY1D(row, ARRAY1D(ju, mu + j)) += lki * ARRAY1D(u, j);
        }
LABEL_9:
        i = ARRAY1D(jrl, i);
        if (i != 0) goto LABEL_7;
//
//  ******  assign kth row of u and diagonal d, set tmp(k)  *************
LABEL_10:
        if (ARRAY1D(row, k) == 0.0) goto LABEL_108;
        dk = 1.0 / ARRAY1D(row, k);
        ARRAY1D(  d, k) = dk;
        ARRAY1D(tmp, k) = sum * dk;
        if (k == n) continue;
        jmin = ARRAY1D(iu, k);
        jmax = ARRAY1D(iu, k + 1) - 1;
        if (jmin > jmax) goto LABEL_12;
        mu = ARRAY1D(iju, k) - jmin;
        for (j = jmin; j <= jmax; ++j) {
            ARRAY1D(u, j) = ARRAY1D(row, ARRAY1D(ju, mu + j)) * dk;
        }
LABEL_12:
//
//  ******  update irl and jrl, keeping jrl in decreasing order  ********
        i = i1;
        if (i == 0) goto LABEL_18;
LABEL_14:
        ARRAY1D(irl, i) += 1;
        i1 = ARRAY1D(jrl, i);
        if (ARRAY1D(irl, i) >= ARRAY1D(il, i + 1)) goto LABEL_17;
        ijlb = ARRAY1D(irl, i) - ARRAY1D(il, i) + ARRAY1D(ijl, i);
        j = ARRAY1D(jl, ijlb);
LABEL_15:
        if (i > ARRAY1D(jrl, j)) goto LABEL_16;
        j = ARRAY1D(jrl, j);
        goto LABEL_15;
LABEL_16:
        ARRAY1D(jrl, i) = ARRAY1D(jrl, j);
        ARRAY1D(jrl, j) = i;
LABEL_17:
        i = i1;
        if (i != 0) goto LABEL_14;
LABEL_18:
        if (ARRAY1D(irl, k) >= ARRAY1D(il, k + 1)) continue;
        j = ARRAY1D(jl, ARRAY1D(ijl, k));
        ARRAY1D(jrl, k) = ARRAY1D(jrl, j);
        ARRAY1D(jrl, j) = k;
    }
//
//  ******  solve  ux = tmp  by back substitution  **********************
    k = n;
    for (i = 1; i <= n; ++i) {
        sum = ARRAY1D(tmp, k);
        jmin = ARRAY1D( iu, k);
        jmax = ARRAY1D( iu, k + 1) - 1;
        if (jmin > jmax) goto LABEL_21;
        mu = ARRAY1D(iju, k) - jmin;
        for (j = jmin; j <= jmax; ++j) {
            sum = sum - ARRAY1D(u, j) * ARRAY1D(tmp, ARRAY1D(ju, mu + j));
        }
LABEL_21:
        ARRAY1D(tmp, k) = sum;
        ARRAY1D(z, ARRAY1D(c, k)) = sum;
        k--;
    }
    flag = 0;
    return;
//
// ** error.. insufficient storage for l
LABEL_104:
    flag = 4* + 1;
    return;
// ** error.. insufficient storage for u
LABEL_107:
    flag = 7*n + 1;
    return;
// ** error.. zero pivot
LABEL_108:
    flag = 8*n + k;
    return;
}


/**
 * @fn NNSC
 * 
 * numerical solution of sparse nonsymmetric system of linear
 * equations given ldu-factorization (compressed pointer storage)
 * 
 *       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
 *       output variables.. z
 * 
 *       parameters used internally..
 * fia   - tmp   - temporary vector which gets result of solving  ly = b.
 *       -           size = n.
 * 
 *  internal variables..
 *    jmin, jmax - indices of the first and last positions in a row of
 *      u or l  to be used.
 */
void Odepack::NNSC(
    int n, int *r, int *c, int *il, int *jl, int *ijl, 
    odepack_cpp_real *l, odepack_cpp_real *d, int *iu, int *ju, int *iju, 
    odepack_cpp_real *u, odepack_cpp_real *z, odepack_cpp_real *b, odepack_cpp_real *tmp
)
{
    int i, k, j, jmin, jmax, ml, mu;
    odepack_cpp_real tmpk, sum;
//
//  ******  set tmp to reordered b  *************************************
    for (k = 1; k <= n; ++k) {
        ARRAY1D(tmp, k) = ARRAY1D(b, ARRAY1D(r, k));
    }
//  ******  solve  ly = b  by forward substitution  *********************
    for (k = 1; k <= n; ++k) {
        jmin = ARRAY1D(il, k);
        jmax = ARRAY1D(il, k + 1) - 1;
        tmpk = - ARRAY1D(d, k) * ARRAY1D(tmp, k);
        ARRAY1D(tmp, k) = - tmpk;
        if (jmin > jmax) continue;
        ml = ARRAY1D(ijl, k) - jmin;
        for (j = jmin; j <= jmax; ++j) {
            ARRAY1D(tmp, ARRAY1D(jl, ml + j)) += tmpk * ARRAY1D(l, j);
        }
    }
//  ******  solve  ux = y  by back substitution  ************************
    k = n;
    for (i = 1; i <= n; ++i) {
        sum = - ARRAY1D(tmp, k);
        jmin = ARRAY1D(iu, k);
        jmax = ARRAY1D(iu, k + 1) - 1;
        if (jmin > jmax) goto LABEL_5;
        mu = ARRAY1D(iju, k) - jmin;
        for (j = jmin; j <= jmax; ++j) {
            sum += ARRAY1D(u, j) * ARRAY1D(tmp, ARRAY1D(ju, mu + j));
        }
LABEL_5:
        ARRAY1D(tmp, k) = - sum;
        ARRAY1D(z, ARRAY1D(c, k)) = - sum;
        k--;
    }
    return;
}


/**
 * @fn NNTC
 * 
 * numeric solution of the transpose of a sparse nonsymmetric system
 * of linear equations given lu-factorization (compressed pointer
 * storage)
 * 
 * 
 *       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
 *       output variables.. z
 * 
 *       parameters used internally..
 * fia   - tmp   - temporary vector which gets result of solving ut y = b
 *       -           size = n.
 * 
 * internal variables..
 *   jmin, jmax - indices of the first and last positions in a row of
 *     u or l  to be used.
 */
void Odepack::NNTC(
    int n, int *r, int *c, int *il, int *jl, 
    int *ijl, odepack_cpp_real *l, odepack_cpp_real *d, int *iu, int *ju, 
    int *iju, odepack_cpp_real *u, odepack_cpp_real *z, odepack_cpp_real *b, odepack_cpp_real *tmp
)
{
    int i, j, jmin, jmax, k, mu, ml;
    odepack_cpp_real tmpk, sum;
//
//  ******  set tmp to reordered b  *************************************
    for (k = 1; k <= n; ++k) {
        ARRAY1D(tmp, k) = ARRAY1D(b, ARRAY1D(c, k));
    }
//  ******  solve  ut y = b  by forward substitution  *******************
    for (k = 1; k <= n; ++k) {
        jmin = ARRAY1D(iu, k);
        jmax = ARRAY1D(iu, k + 1) - 1;
        tmpk = - ARRAY1D(tmp, k);
        if (jmin > jmax) continue;
        mu   = ARRAY1D(iju, k) - jmin;
        for (j = jmin; j <= jmax; ++j) {
            ARRAY1D(tmp, ARRAY1D(ju, mu + j)) += tmpk * ARRAY1D(u, j);
        }
    }
//  ******  solve  lt x = y  by back substitution  **********************
    k = n;
    for (i = 1; i <= n; ++i) {
        sum = - ARRAY1D(tmp, k);
        jmin = ARRAY1D(il, k);
        jmax = ARRAY1D(il, k + 1) - 1;
        if (jmin > jmax) goto LABEL_5;
        ml   = ARRAY1D(ijl, k) - jmin;
        for (j = jmin; j <= jmax; ++j) {
            sum += ARRAY1D(l, j) * ARRAY1D(tmp, ARRAY1D(jl, ml + j));
        }
LABEL_5:
        ARRAY1D(tmp, k)          = - sum * ARRAY1D(d, k);
        ARRAY1D(z, ARRAY1D(r, k)) = ARRAY1D(tmp, k);
        k--;
    }
    return;
}


/**
 * @fn DSTODA
 * 
C DSTODA performs one step of the integration of an initial value
C problem for a system of ordinary differential equations.
C Note: DSTODA is independent of the value of the iteration method
C indicator MITER, when this is .ne. 0, and hence is independent
C of the type of chord method used, or the Jacobian structure.
C Communication with DSTODA is done with the following variables:
C
C Y      = an array of length .ge. N used as the Y argument in
C          all calls to F and JAC.
C NEQ    = integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls to F and JAC.
C YH     = an NYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C          j-th derivative of y(i), scaled by H**j/factorial(j)
C          (j = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C NYH    = a constant integer .ge. N, the first dimension of YH.
C YH1    = a one-dimensional array occupying the same space as YH.
C EWT    = an array of length N containing multiplicative weights
C          for local error measurements.  Local errors in y(i) are
C          compared to 1.0/EWT(i) in various error tests.
C SAVF   = an array of working storage, of length N.
C ACOR   = a work array of length N, used for the accumulated
C          corrections.  On a successful return, ACOR(i) contains
C          the estimated one-step local error in y(i).
C WM,IWM = real and integer work arrays associated with matrix
C          operations in chord iteration (MITER .ne. 0).
C PJAC   = name of routine to evaluate and preprocess Jacobian matrix
C          and P = I - H*EL0*Jac, if a chord method is being used.
C          It also returns an estimate of norm(Jac) in PDNORM.
C SLVS   = name of routine to solve linear system in chord iteration.
C CCMAX  = maximum relative change in H*EL0 before PJAC is called.
C H      = the step size to be attempted on the next step.
C          H is altered by the error control algorithm during the
C          problem.  H can be either positive or negative, but its
C          sign must remain constant throughout the problem.
C HMIN   = the minimum absolute value of the step size H to be used.
C HMXI   = inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C          HMIN and HMXI may be changed at any time, but will not
C          take effect until the next change of H is considered.
C TN     = the independent variable. TN is updated on each step taken.
C JSTART = an integer used for input only, with the following
C          values and meanings:
C               0  perform the first step.
C           .gt.0  take a new step continuing from the last.
C              -1  take the next step with a new value of H,
C                    N, METH, MITER, and/or matrix parameters.
C              -2  take the next step with a new value of H,
C                    but with other inputs unchanged.
C          On return, JSTART is set to 1 to facilitate continuation.
C KFLAG  = a completion code with the following meanings:
C               0  the step was succesful.
C              -1  the requested error could not be achieved.
C              -2  corrector convergence could not be achieved.
C              -3  fatal error in PJAC or SLVS.
C          A return with KFLAG = -1 or -2 means either
C          ABS(H) = HMIN or 10 consecutive failures occurred.
C          On a return with KFLAG negative, the values of TN and
C          the YH array are as of the beginning of the last
C          step, and H is the last step size attempted.
C MAXORD = the maximum order of integration method to be allowed.
C MAXCOR = the maximum number of corrector iterations allowed.
C MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
C MXNCF  = maximum number of convergence failures allowed.
C METH   = current method.
C          METH = 1 means Adams method (nonstiff)
C          METH = 2 means BDF method (stiff)
C          METH may be reset by DSTODA.
C MITER  = corrector iteration method.
C          MITER = 0 means functional iteration.
C          MITER = JT .gt. 0 means a chord iteration corresponding
C          to Jacobian type JT.  (The DLSODA/DLSODAR argument JT is
C          communicated here as JTYP, but is not used in DSTODA
C          except to load MITER following a method switch.)
C          MITER may be reset by DSTODA.
C N      = the number of first-order differential equations.
 */
void Odepack::DSTODA(
    int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *yh1, 
    odepack_cpp_real *ewt, odepack_cpp_real *savf, odepack_cpp_real *acor, odepack_cpp_real *wm, int *iwm, 
    ODEPACK_FUNCTION f, ODEPACK_JACOBIAN1 jac, FUNC_PJAC<ODEPACK_JACOBIAN1> pjac, FUNC_SLVS slvs, 
    void *user_data
)
{
#ifndef YH
#define YH(i, j) ARRAY2D(yh, nyh, i, j)
#endif
#ifndef ELCO
#define ELCO(i, j) ARRAY2D(elco, 13, i, j)
#endif
#ifndef TESCO
#define TESCO(i, j) ARRAY2D(tesco, 3, i, j)
#endif
// DLS001
    odepack_cpp_real &conit = dls1_.conit, &crate = dls1_.crate, *el = dls1_.el, *elco = dls1_.elco,
        &hold = dls1_.hold, &rmax = dls1_.rmax, *tesco = dls1_.tesco,
        &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &ialth = dls1_.ialth, &ipup = dls1_.ipup, &lmax = dls1_.lmax, &meo = dls1_.meo, &nqnyh = dls1_.nqnyh, &nslp = dls1_.nslp,
        &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
// DLSA01
    odepack_cpp_real *cm1 = dlsa_.cm1, *cm2 = dlsa_.cm2, &pdest = dlsa_.pdest, &pdlast = dlsa_.pdlast, &ratio = dlsa_.ratio,
    &pdnorm = dlsa_.pdnorm;
    int &icount = dlsa_.icount, &irflag = dlsa_.irflag, &jtyp = dlsa_.jtyp, &mused = dlsa_.mused, &mxordn = dlsa_.mxordn, &mxords = dlsa_.mxords;
//
    int i, j, i1, jb, iredo, iret, m, ncf, newq;
    int lm1, lm1p1, lm2, lm2p1, nqm1, nqm2;
    odepack_cpp_real dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup,
        r, rh, rhdn, rhsm, rhup, told, dmnorm;
    odepack_cpp_real alpha, dm1, dm2, exm1, exm2,
        pdh, pnorm, rate, rh1, rh1it, rh2, rm;
    odepack_cpp_real sm1[12] = {
        0.5, 0.575, 0.55, 0.45, 0.35, 0.25, 0.20, 0.15, 0.10, 0.075, 0.050, 0.025
    };
//
    kflag = 0;
    told = tn;
    ncf = 0;
    ierpj = 0;
    iersl = 0;
    jcur = 0;
    icf = 0;
    delp = 0.0;
    if (jstart > 0) goto LABEL_200;
    if (jstart == -1) goto LABEL_100;
    if (jstart == -2) goto LABEL_160;
//-----------------------------------------------------------------------
// On the first call, the order is set to 1, and other variables are
// initialized.  RMAX is the maximum ratio by which H can be increased
// in a single step.  It is initially 1.E4 to compensate for the small
// initial H, but then is normally equal to 10.  If a failure
// occurs (in corrector convergence or error test), RMAX is set at 2
// for the next increase.
// DCFODE is called to get the needed coefficients for both methods.
//-----------------------------------------------------------------------
    lmax = maxord + 1;
    nq = 1;
    l = 2;
    ialth = 2;
    rmax = 10000.0;
    rc = 0.0;
    el0 = 1.0;
    crate = 0.7;
    hold = h;
    nslp = 0;
    ipup = miter;
    iret = 3;
// Initialize switching parameters.  METH = 1 is assumed initially. -----
    icount = 20;
    irflag = 0;
    pdest = 0.0;
    pdlast = 0.0;
    ratio = 5.0;
    DCFODE(2, elco, tesco);
    for (i = 1; i <= 5; ++i) {
        ARRAY1D(cm2, i) = TESCO(2, i) * ELCO(i + 1, i);
    }
    DCFODE(1, elco, tesco);
    for (i = 1; i <= 12; ++i) {
        ARRAY1D(cm1, i) = TESCO(2, i) * ELCO(i + 1, i);
    }
    goto LABEL_150;
//-----------------------------------------------------------------------
// The following block handles preliminaries needed when JSTART = -1.
// IPUP is set to MITER to force a matrix update.
// If an order increase is about to be considered (IALTH = 1),
// IALTH is reset to 2 to postpone consideration one more step.
// If the caller has changed METH, DCFODE is called to reset
// the coefficients of the method.
// If H is to be changed, YH must be rescaled.
// If H or METH is being changed, IALTH is reset to L = NQ + 1
// to prevent further changes in H for that many steps.
//-----------------------------------------------------------------------
LABEL_100:
    ipup = miter;
    lmax = maxord + 1;
    if (ialth == 1) ialth = 2;
    if (meth == mused) goto LABEL_160;
    DCFODE(meth, elco, tesco);
    ialth = l;
    iret = 1;
//-----------------------------------------------------------------------
// The el vector and related constants are reset
// whenever the order NQ is changed, or at the start of the problem.
//-----------------------------------------------------------------------
LABEL_150:
    for (i = 1; i <= l; ++i) {
        ARRAY1D(el, i) = ELCO(i, nq);
    }
    nqnyh = nq * nyh;
    rc *= ARRAY1D(el, 1) / el0;
    el0 = ARRAY1D(el, 1);
    conit = 0.5 / static_cast<odepack_cpp_real>(nq + 2);
    if (iret == 1) {
        goto LABEL_160;
    } else if (iret == 2) {
        goto LABEL_170;
    } else if (iret == 3) {
        goto LABEL_200;
    }
//-----------------------------------------------------------------------
// If H is being changed, the H ratio RH is checked against
// RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
// L = NQ + 1 to prevent a change of H for that many steps, unless
// forced by a convergence or error test failure.
//-----------------------------------------------------------------------
LABEL_160:
    if (h == hold) goto LABEL_200;
    rh = h / hold;
    h = hold;
    iredo = 3;
    goto LABEL_175;
LABEL_170:
    rh = std::max(rh, hmin / std::abs(h));
LABEL_175:
    rh = std::min(rh, rmax);
    rh /= std::max(1.0, std::abs(h) * hmxi * rh);
//-----------------------------------------------------------------------
// If METH = 1, also restrict the new step size by the stability region.
// If this reduces H, set IRFLAG to 1 so that if there are roundoff
// problems later, we can assume that is the cause of the trouble.
//-----------------------------------------------------------------------
    if (meth == 2) goto LABEL_178;
    irflag = 0;
    pdh = std::max(std::abs(h) * pdlast, 0.000001);
    if (rh * pdh * 1.00001 < ARRAY1D(sm1, nq)) goto LABEL_178;
    rh = ARRAY1D(sm1, nq) / pdh;
    irflag = 1;
LABEL_178:
    r = 1.0;
    for (j = 2; j <= l; ++j) {
        r *= rh;
        for (i = 1; i <= n; ++i) {
            YH(i, j) *= r;
        }
    }
    h *= rh;
    rc *= rh;
    ialth = l;
    if (iredo == 0) goto LABEL_690;
//-----------------------------------------------------------------------
// This section computes the predicted values by effectively
// multiplying the YH array by the Pascal triangle matrix.
// RC is the ratio of new to old values of the coefficient  H*EL(1).
// When RC differs from 1 by more than CCMAX, IPUP is set to MITER
// to force PJAC to be called, if a Jacobian is involved.
// In any case, PJAC is called at least every MSBP steps.
//-----------------------------------------------------------------------
LABEL_200:
    if (std::abs(rc - 1.0) > ccmax) ipup = miter; 
    if (nst >= nslp + msbp) ipup = miter;
    tn += h;
    i1 = nqnyh + 1;
    for (jb = 1; jb <= nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= nqnyh; ++i) {
            ARRAY1D(yh1, i) += ARRAY1D(yh1, i + nyh);
        }
    }
    pnorm = DMNORM(n, yh1, ewt);
//-----------------------------------------------------------------------
// Up to MAXCOR corrector iterations are taken.  A convergence test is
// made on the RMS-norm of each correction, weighted by the error
// weight vector EWT.  The sum of the corrections is accumulated in the
// vector ACOR(i).  The YH array is not altered in the corrector loop.
//-----------------------------------------------------------------------
LABEL_220:
    m = 0;
    rate = 0.0;
    del = 0.0;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = YH(i, 1);
    }
    (*f)(neq, tn, y, savf, user_data);
    nfe++;
    if (ipup <= 0) goto LABEL_250;
//-----------------------------------------------------------------------
// If indicated, the matrix P = I - H*EL(1)*J is reevaluated and
// preprocessed before starting the corrector iteration.  IPUP is set
// to 0 as an indicator that this has been done.
//-----------------------------------------------------------------------
    (this->*pjac)(neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac, user_data);
    ipup = 0;
    rc = 1.0;
    nslp = nst;
    crate = 0.7;
    if (ierpj != 0) goto LABEL_430;
LABEL_250:
    for (i = 1; i <= n; ++i) {
        ARRAY1D(acor, i) = 0.0;
    }
LABEL_270:
    if (miter != 0) goto LABEL_350;
//-----------------------------------------------------------------------
// In the case of functional iteration, update Y directly from
// the result of the last function evaluation.
//-----------------------------------------------------------------------
    for (i = 1; i <= n; ++i) {
        ARRAY1D(savf, i) = h * ARRAY1D(savf, i) - YH(i, 2);
        ARRAY1D(y, i) = ARRAY1D(savf, i) - ARRAY1D(acor, i);
    }
    del = DMNORM(n, y, ewt);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = YH(i, 1) + ARRAY1D(el, 1) * ARRAY1D(savf, i);
        ARRAY1D(acor, i) = ARRAY1D(savf, i);
    }
    goto LABEL_400;
//-----------------------------------------------------------------------
// In the case of the chord method, compute the corrector error,
// and solve the linear system with that as right-hand side and
// P as coefficient matrix.
//-----------------------------------------------------------------------
LABEL_350:
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = h * ARRAY1D(savf, i) - (YH(i, 2) + ARRAY1D(acor, i));
    }
    (this->*slvs)(wm, iwm, y, savf);
    if (iersl < 0) goto LABEL_430;
    if (iersl > 0) goto LABEL_410;
    del = DMNORM(n, y, ewt);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(acor, i) += ARRAY1D(y, i);
        ARRAY1D(y, i) = YH(i, 1) + ARRAY1D(el, 1) * ARRAY1D(acor, i);
    }
//-----------------------------------------------------------------------
// Test for convergence.  If M .gt. 0, an estimate of the convergence
// rate constant is stored in CRATE, and this is used in the test.
//
// We first check for a change of iterates that is the size of
// roundoff error.  If this occurs, the iteration has converged, and a
// new rate estimate is not formed.
// In all other cases, force at least two iterations to estimate a
// local Lipschitz constant estimate for Adams methods.
// On convergence, form PDEST = local maximum Lipschitz constant
// estimate.  PDLAST is the most recent nonzero estimate.
//-----------------------------------------------------------------------
LABEL_400:
    if (del <= 100.0 * pnorm * uround) goto LABEL_450;
    if (m == 0 && meth == 1) goto LABEL_405;
    if (m == 0) goto LABEL_402;
    rm = 1024.0;
    if (del <= 1024.0 * delp) rm = del / delp;
    rate = std::max(rate, rm);
    crate = std::max(0.2 * crate, rm);
LABEL_402:
    dcon = del * std::min(1.0, 1.5 * crate) / (TESCO(2, nq) * conit);
    if (dcon > 1.0) goto LABEL_405;
    pdest = std::max(pdest, rate / std::abs(h * ARRAY1D(el, 1)));
    if (pdest != 0.0) pdlast = pdest;
    goto LABEL_450;
LABEL_405:
    m++;
    if (m == maxcor) goto LABEL_410;
    if (m >= 2 && del > 2.0 * delp) goto LABEL_410;
    delp = del;
    (*f)(neq, tn, y, savf, user_data);
    nfe++;
    goto LABEL_270;
//-----------------------------------------------------------------------
// The corrector iteration failed to converge.
// If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
// the next try.  Otherwise the YH array is retracted to its values
// before prediction, and H is reduced, if possible.  If H cannot be
// reduced or MXNCF failures have occurred, exit with KFLAG = -2.
//-----------------------------------------------------------------------
LABEL_410:
    if (miter == 0 || jcur == 1) goto LABEL_430;
    icf = 1;
    ipup = miter;
    goto LABEL_220;
LABEL_430:
    icf = 2;
    ncf++;
    rmax = 2.0;
    tn = told;
    i1 = nqnyh + 1;
    for (jb = 1; jb <= nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= nqnyh; ++i) {
            ARRAY1D(yh1, i) -= ARRAY1D(yh1, i + nyh);
        }
    }
    if (ierpj < 0 || iersl < 0) goto LABEL_680;
    if (std::abs(h) <= hmin * 1.00001) goto LABEL_670;
    if (ncf == mxncf) goto LABEL_670;
    rh = 0.25;
    ipup = miter;
    iredo = 1;
    goto LABEL_170;
//-----------------------------------------------------------------------
// The corrector has converged.  JCUR is set to 0
// to signal that the Jacobian involved may need updating later.
// The local error test is made and control passes to statement 500
// if it fails.
//-----------------------------------------------------------------------
LABEL_450:
    jcur = 0;
    if (m == 0) dsm = del / TESCO(2, nq);
    if (m > 0)  dsm = DMNORM(n, acor, ewt) / TESCO(2, nq);
    if (dsm > 1.0) goto LABEL_500;
//-----------------------------------------------------------------------
// After a successful step, update the YH array.
// Decrease ICOUNT by 1, and if it is -1, consider switching methods.
// If a method switch is made, reset various parameters,
// rescale the YH array, and exit.  If there is no switch,
// consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
// If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
// use in a possible order increase on the next step.
// If a change in H is considered, an increase or decrease in order
// by one is considered also.  A change in H is made only if it is by a
// factor of at least 1.1.  If not, IALTH is set to 3 to prevent
// testing for that many steps.
//-----------------------------------------------------------------------
    kflag = 0;
    iredo = 0;
    nst++;
    hu = h;
    nqu = nq;
    mused = meth;
    for (j = 1; j <= l; ++j) {
        for (i = 1; i <= n; ++i) {
            YH(i, j) += ARRAY1D(el, j) * ARRAY1D(acor, i);
        }
    }
    icount--;
    if (icount >= 0) goto LABEL_488;
    if (meth == 2) goto LABEL_480;
//-----------------------------------------------------------------------
// We are currently using an Adams method.  Consider switching to BDF.
// If the current order is greater than 5, assume the problem is
// not stiff, and skip this section.
// If the Lipschitz constant and error estimate are not polluted
// by roundoff, go to 470 and perform the usual test.
// Otherwise, switch to the BDF methods if the last step was
// restricted to insure stability (irflag = 1), and stay with Adams
// method if not.  When switching to BDF with polluted error estimates,
// in the absence of other information, odepack_cpp_real the step size.
//
// When the estimates are OK, we make the usual test by computing
// the step size we could have (ideally) used on this step,
// with the current (Adams) method, and also that for the BDF.
// If NQ .gt. MXORDS, we consider changing to order MXORDS on switching.
// Compare the two step sizes to decide whether to switch.
// The step size advantage must be at least RATIO = 5 to switch.
//-----------------------------------------------------------------------
    if (nq > 5) goto LABEL_488;
    if (dsm > 100.0 * pnorm * uround && pdest != 0.0) goto LABEL_470;
    if (irflag == 0) goto LABEL_488;
    rh2 = 2.0;
    nqm2 = std::min(nq, mxords);
    goto LABEL_478;
LABEL_470:
    exsm = 1.0 / static_cast<odepack_cpp_real>(l);
    rh1 = 1.0 / (1.2 * std::pow(dsm, exsm) + 0.0000012);
    rh1it = 2.0 * rh1;
    pdh = pdlast * std::abs(h);
    if (pdh * rh1 > 0.00001) rh1it = ARRAY1D(sm1, nq) / pdh;
    rh1 = std::min(rh1, rh1it);
    if (nq <= mxords) goto LABEL_474;
    nqm2 = mxords;
    lm2 = mxords + 1;
    exm2 = 1.0 / static_cast<odepack_cpp_real>(lm2);
    lm2p1 = lm2 + 1;
    dm2 = DMNORM(n, &YH(1, lm2p1), ewt) / ARRAY1D(cm2, mxords);
    rh2 = 1.0 / (1.2 * std::pow(dm2, exm2) + 0.0000012);
    goto LABEL_476;
LABEL_474:
    dm2 = dsm * (ARRAY1D(cm1, nq) / ARRAY1D(cm2, nq));
    rh2 = 1.0 / (1.2 * std::pow(dm2, exsm) + 0.0000012);
    nqm2 = nq;
LABEL_476:
    if (rh2 < ratio * rh1) goto LABEL_488;
// THE SWITCH TEST PASSED.  RESET RELEVANT QUANTITIES FOR BDF. ----------
LABEL_478:
    rh = rh2;
    icount = 20;
    meth = 2;
    miter = jtyp;
    pdlast = 0.0;
    nq = nqm2;
    l = nq + 1;
    goto LABEL_170;
//-----------------------------------------------------------------------
// We are currently using a BDF method.  Consider switching to Adams.
// Compute the step size we could have (ideally) used on this step,
// with the current (BDF) method, and also that for the Adams.
// If NQ .gt. MXORDN, we consider changing to order MXORDN on switching.
// Compare the two step sizes to decide whether to switch.
// The step size advantage must be at least 5/RATIO = 1 to switch.
// If the step size for Adams would be so small as to cause
// roundoff pollution, we stay with BDF.
//-----------------------------------------------------------------------
LABEL_480:
    exsm = 1.0 / static_cast<odepack_cpp_real>(l);
    if (mxordn >= nq) goto LABEL_484;
    nqm1 = mxordn;
    lm1 = mxordn + 1;
    exm1 = 1.0 / static_cast<odepack_cpp_real>(lm1);
    lm1p1 = lm1 + 1;
    dm1 = DMNORM(n, &YH(1, lm1p1), ewt) / ARRAY1D(cm1, mxordn);
    rh1 = 1.0 / (1.2 * std::pow(dm1, exm1) + 0.0000012);
    goto LABEL_486;
LABEL_484:
    dm1 = dsm * (ARRAY1D(cm2, nq) / ARRAY1D(cm1, nq));
    rh1 = 1.0 / (1.2 * std::pow(dm1, exsm) + 0.0000012);
    nqm1 = nq;
    exm1 = exsm;
LABEL_486:
    rh1it = 2.0 * rh1;
    pdh = pdnorm * std::abs(h);
    if (pdh * rh1 > 0.00001) rh1it = ARRAY1D(sm1, nqm1) / pdh;
    rh1 = std::min(rh1, rh1it);
    rh2 = 1.0 / (1.2 * std::pow(dsm, exsm) + 0.0000012);
    if (rh1 * ratio < 5.0 * rh2) goto LABEL_488;
    alpha = std::max(0.001, rh1);
    dm1 = std::pow(alpha, exm1) * dm1;
    if (dm1 <= 1000.0 * uround * pnorm) goto LABEL_488;
// The switch test passed.  Reset relevant quantities for Adams. --------
    rh = rh1;
    icount = 20;
    meth = 1;
    miter = 0;
    pdlast = 0.0;
    nq = nqm1;
    l = nq + 1;
    goto LABEL_170;
//
// No method switch is being made.  Do the usual step/order selection. --
LABEL_488:
    ialth--;
    if (ialth == 0) goto LABEL_520;
    if (ialth > 1) goto LABEL_700;
    if (l == lmax) goto LABEL_700;
    for (i = 1; i <= n; ++i) {
        YH(i, lmax) = ARRAY1D(acor, i);
    }
    goto LABEL_700;
//-----------------------------------------------------------------------
// The error test failed.  KFLAG keeps track of multiple failures.
// Restore TN and the YH array to their previous values, and prepare
// to try the step again.  Compute the optimum step size for this or
// one lower order.  After 2 or more failures, H is forced to decrease
// by a factor of 0.2 or less.
//-----------------------------------------------------------------------
LABEL_500:
    kflag--;
    tn = told;
    i1 = nqnyh + 1;
    for (jb = 1; jb <= nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= nqnyh; ++i) {
            ARRAY1D(yh1, i) -= ARRAY1D(yh1, i + nyh);
        }
    }
    rmax = 2.0;
    if (std::abs(h) <= hmin * 1.00001) goto LABEL_660;
    if (kflag <= -3) goto LABEL_640;
    iredo = 2;
    rhup = 0.0;
    goto LABEL_540;
//-----------------------------------------------------------------------
// Regardless of the success or failure of the step, factors
// RHDN, RHSM, and RHUP are computed, by which H could be multiplied
// at order NQ - 1, order NQ, or order NQ + 1, respectively.
// In the case of failure, RHUP = 0.0 to avoid an order increase.
// The largest of these is determined and the new order chosen
// accordingly.  If the order is to be increased, we compute one
// additional scaled derivative.
//-----------------------------------------------------------------------
LABEL_520:
    rhup = 0.0;
    if (l == lmax) goto LABEL_540;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(savf, i) = ARRAY1D(acor, i) - YH(i, lmax);
    }
    dup = DMNORM(n, savf, ewt) / TESCO(3, nq);
    exup = 1.0 / static_cast<odepack_cpp_real>(l + 1);
    rhup = 1.0 / (1.4 * std::pow(dup, exup) + 0.0000014);
LABEL_540:
    exsm = 1.0 / static_cast<odepack_cpp_real>(l);
    rhsm = 1.0 / (1.2 * std::pow(dsm, exsm) + 0.0000012);
    rhdn = 0.0;
    if (nq == 1) goto LABEL_550;
    ddn = DMNORM(n, &YH(1, l), ewt) / TESCO(1, nq);
    exdn = 1.0 / static_cast<odepack_cpp_real>(nq);
    rhdn = 1.0 / (1.3 * std::pow(ddn, exdn) + 0.0000013);
// If METH = 1, limit RH according to the stability region also. --------
LABEL_550:
    if (meth == 2) goto LABEL_560;
    pdh = std::max(std::abs(h) * pdlast, 0.000001);
    if (l < lmax) rhup = std::min(rhup, ARRAY1D(sm1, l) / pdh);
    rhsm = std::min(rhsm, ARRAY1D(sm1, nq) / pdh);
    if (nq > 1) rhdn = std::min(rhdn, ARRAY1D(sm1, nq - 1) / pdh);
    pdest = 0.0;
LABEL_560:
    if (rhsm >= rhup) goto LABEL_570;
    if (rhup > rhdn) goto LABEL_590;
    goto LABEL_580;
LABEL_570:
    if (rhsm < rhdn) goto LABEL_580;
    newq = nq;
    rh = rhsm;
    goto LABEL_620;
LABEL_580:
    newq = nq - 1;
    rh = rhdn;
    if (kflag < 0 && rh > 1.0) rh = 1.0;
    goto LABEL_620;
LABEL_590:
    newq = l;
    rh = rhup;
    if (rh < 1.1) goto LABEL_610;
    r = ARRAY1D(el, l) / static_cast<odepack_cpp_real>(l);
    for (i = 1; i <= n; ++i) {
        YH(i, newq + 1) = ARRAY1D(acor, i) * r;
    }
    goto LABEL_630;
LABEL_610:
    ialth = 3;
    goto LABEL_700;
// If METH = 1 and H is restricted by stability, bypass 10 percent test.
LABEL_620:
    if (meth == 2) goto LABEL_622;
    if (rh * pdh * 1.00001 >= ARRAY1D(sm1, newq)) goto LABEL_625;
LABEL_622:
    if (kflag == 0 && rh < 1.1) goto LABEL_610;
LABEL_625:
    if (kflag <= -2) rh = std::min(rh, 0.2);
//-----------------------------------------------------------------------
// If there is a change of order, reset NQ, L, and the coefficients.
// In any case H is reset according to RH and the YH array is rescaled.
// Then exit from 690 if the step was OK, or redo the step otherwise.
//-----------------------------------------------------------------------
    if (newq == nq) goto LABEL_170;
LABEL_630:
    nq = newq;
    l = nq + 1;
    iret = 2;
    goto LABEL_150;
//-----------------------------------------------------------------------
// Control reaches this section if 3 or more failures have occured.
// If 10 failures have occurred, exit with KFLAG = -1.
// It is assumed that the derivatives that have accumulated in the
// YH array have errors of the wrong order.  Hence the first
// derivative is recomputed, and the order is set to 1.  Then
// H is reduced by a factor of 10, and the step is retried,
// until it succeeds or H reaches HMIN.
//-----------------------------------------------------------------------
LABEL_640:
    if (kflag == -10) goto LABEL_660;
    rh = 0.1;
    rh = std::max(hmin / std::abs(h), rh);
    h *= rh;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = YH(i, 1);
    }
    (*f)(neq, tn, y, savf, user_data);
    nfe++;
    for (i = 1; i <= n; ++i) {
        YH(i, 2) = h * ARRAY1D(savf, i);
    }
    ipup = miter;
    ialth = 5;
    if (nq == 1) goto LABEL_200;
    nq = 1;
    l = 2;
    iret = 3;
    goto LABEL_150;
//-----------------------------------------------------------------------
// All returns are made through this section.  H is saved in HOLD
// to allow the caller to change H on the next step.
//-----------------------------------------------------------------------
LABEL_660:
    kflag = -1;
    goto LABEL_720;
LABEL_670:
    kflag = -2;
    goto LABEL_720;
LABEL_680:
    kflag = -3;
    goto LABEL_720;
LABEL_690:
    rmax = 10.0;
LABEL_700:
    r = 1.0 / TESCO(2, nqu);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(acor, i) *= r;
    }
LABEL_720:
    hold = h;
    jstart = 1;
    return;
//
#ifdef YH
#undef YH
#endif
#ifdef ELCO
#undef ELCO
#endif
#ifdef TESCO
#undef TESCO
#endif
}


/**
 * @fn DPRJA
 * 
C DPRJA is called by DSTODA to compute and process the matrix
C P = I - H*EL(1)*J , where J is an approximation to the Jacobian.
C Here J is computed by the user-supplied routine JAC if
C MITER = 1 or 4 or by finite differencing if MITER = 2 or 5.
C J, scaled by -H*EL(1), is stored in WM.  Then the norm of J (the
C matrix norm consistent with the weighted max-norm on vectors given
C by DMNORM) is computed, and J is overwritten by P.  P is then
C subjected to LU decomposition in preparation for later solution
C of linear systems with P as coefficient matrix.  This is done
C by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
C
C In addition to variables described previously, communication
C with DPRJA uses the following:
C Y     = array containing predicted values on entry.
C FTEM  = work array of length N (ACOR in DSTODA).
C SAVF  = array containing f evaluated at predicted y.
C WM    = real work space for matrices.  On output it contains the
C         LU decomposition of P.
C         Storage of matrix elements starts at WM(3).
C         WM also contains the following matrix-related data:
C         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
C IWM   = integer work space containing pivot information, starting at
C         IWM(21).   IWM also contains the band parameters
C         ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C EL0   = EL(1) (input).
C PDNORM= norm of Jacobian matrix. (Output).
C IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
C         P matrix found to be singular.
C JCUR  = output flag = 1 to indicate that the Jacobian matrix
C         (or approximation) is now current.
C This routine also uses the Common variables EL0, H, TN, UROUND,
C MITER, N, NFE, and NJE.
 */
void Odepack::DPRJA(int neq,  odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *ewt,  
        odepack_cpp_real *ftem, odepack_cpp_real *savf, odepack_cpp_real *wm, int *iwm, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN1 jac, 
        void *user_data)
{
#ifndef YH
#define YH(i, j) ARRAY2D(yh, nyh, i, j)
#endif
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
// DLSA01
    odepack_cpp_real &pdnorm = dlsa_.pdnorm;
    int &jtyp = dlsa_.jtyp, &mused = dlsa_.mused, &mxordn = dlsa_.mxordn, &mxords = dlsa_.mxords;
//
    int i, i1, i2, ier, ii, j, j1, jj, lenp, 
        mba, mband, meb1, meband, ml, ml3, mu, np1;
    odepack_cpp_real con, fac, hl0, r, r0, srur, yi, yj, yjj;
//
    nje++;
    ierpj = 0;
    jcur = 1;
    hl0 = h * el0;
    if (miter == 1) {
        goto LABEL_100;
    } else if (miter == 2) {
        goto LABEL_200;
    } else if (miter == 3) {
        goto LABEL_300;
    } else if (miter == 4) {
        goto LABEL_400;
    } else if (miter == 5) {
        goto LABEL_500;
    }
// If MITER = 1, call JAC and multiply by scalar. -----------------------
LABEL_100:
    lenp = n * n;
    for (i = 1; i <= lenp; ++i) {
        ARRAY1D(wm, i + 2) = 0.0;
    }
    (*jac)(neq, tn, y, 0, 0, &ARRAY1D(wm, 3), n, user_data);
    con = -hl0;
    for (i = 1; i <= lenp; ++i) {
        ARRAY1D(wm, i + 2) *= con;
    }
    goto LABEL_240;
// If MITER = 2, make N calls to F to approximate J. --------------------
LABEL_200:
    fac = DMNORM(n, savf, ewt);
    r0 = 1000.0 * std::abs(h) * uround * static_cast<odepack_cpp_real>(n) * fac;
    if (r0 == 0.0) r0 = 1.0;
    srur = ARRAY1D(wm, 1);
    j1 = 2;
    for (j = 1; j <= n; ++j) {
        yj = ARRAY1D(y, j);
        r = std::max(srur * std::abs(yj), r0 / ARRAY1D(ewt, j));
        ARRAY1D(y, j) += r;
        fac = -hl0 / r;
        (*f)(neq, tn, y, ftem, user_data);
        for (i = 1; i <= n; ++i) {
            ARRAY1D(wm, i + j1) = (ARRAY1D(ftem, i) - ARRAY1D(savf, i)) * fac;
        }
        ARRAY1D(y, j) = yj;
        j1 += n;
    }
    nfe += n;
LABEL_240:
// Compute norm of Jacobian. --------------------------------------------
    pdnorm = DFNORM(n, &ARRAY1D(wm, 3), ewt) / std::abs(hl0);
// Add identity matrix. -------------------------------------------------
    j = 3;
    np1 = n + 1;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(wm, j) += 1.0;
        j += np1;
    }
// Do LU decomposition on P. --------------------------------------------
    DGEFA(&ARRAY1D(wm, 3), n, n, &ARRAY1D(iwm, 21), ier);
    if (ier != 0) ierpj = 1;
    return;
// Dummy block only, since MITER is never 3 in this routine. ------------
LABEL_300:
    return;
// If MITER = 4, call JAC and multiply by scalar. -----------------------
LABEL_400:
    ml = ARRAY1D(iwm, 1);
    mu = ARRAY1D(iwm, 2);
    ml3 = ml + 3;
    mband = ml + mu + 1;
    meband = mband + ml;
    lenp = meband * n;
    for (i = 1; i <= lenp; ++i) {
        ARRAY1D(wm, i + 2) = 0.0;
    }
    (*jac)(neq, tn, y, ml, mu, &ARRAY1D(wm, ml3), meband, user_data);
    con = -hl0;
    for (i = 1; i <= lenp; ++i) {
        ARRAY1D(wm, i + 2) *= con;
    }
    goto LABEL_570;
// If MITER = 5, make MBAND calls to F to approximate J. ----------------
LABEL_500:
    ml = ARRAY1D(iwm, 1);
    mu = ARRAY1D(iwm, 2);
    mband = ml + mu + 1;
    mba = std::min(mband, n);
    meband = mband + ml;
    meb1 = meband - 1;
    srur = ARRAY1D(wm, 1);
    fac = DMNORM(n, savf, ewt);
    r0 = 1000.0 * std::abs(h) * uround * static_cast<odepack_cpp_real>(n) * fac;
    if (r0 == 0.0) r0 = 1.0;
    for (j = 1; j <= mba; ++j) {
        for (i = j; j <= n; j += mband) {
            yi = ARRAY1D(y, i);
            r = std::max(srur * std::abs(yi), r0 / ARRAY1D(ewt, i));
            ARRAY1D(y, i) += r;
        }
        (*f)(neq, tn, y, ftem, user_data);
        for (jj = j; jj <= n; jj += mband) {
            ARRAY1D(y, jj) = YH(jj, 1);
            yjj = ARRAY1D(y, jj);
            r = std::max(srur * std::abs(yjj), r0 / ARRAY1D(ewt, jj));
            fac = -hl0 / r;
            i1 = std::max(jj - mu, 1);
            i2 = std::min(jj + ml, n);
            ii = jj * meb1 - ml + 2;
            for (i = i1; i <= i2; ++i) {
                ARRAY1D(wm, ii + i) = (ARRAY1D(ftem, i) - ARRAY1D(savf, i)) * fac;
            }
        }
    }
    nfe += mba;
LABEL_570:
// Compute norm of Jacobian. --------------------------------------------
    pdnorm = DBNORM(n, &ARRAY1D(wm, ml + 3), meband, ml, mu, ewt) / std::abs(hl0);
// Add identity matrix. -------------------------------------------------
    ii = mband + 2;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(wm, ii) += 1.0;
        ii += mband;
    }
// Do LU decomposition of P. --------------------------------------------
    DGBFA(&ARRAY1D(wm, 3), meband, n, ml, mu, &ARRAY1D(iwm, 3), ier);
    if (ier != 0) ierpj = 1;
    return;
//
#ifdef YH
#undef YH
#endif
}


/**
 * @fn DMNORM
 * 
C This function routine computes the weighted max-norm
C of the vector of length N contained in the array V, with weights
C contained in the array w of length N:
C   DMNORM = MAX(i=1,...,N) ABS(V(i))*W(i)
 */
odepack_cpp_real Odepack::DMNORM(int n, odepack_cpp_real *v, odepack_cpp_real *w)
{
    odepack_cpp_real vm = 0.0;
    for (int i = 0; i < n; ++i) {
        vm = std::max(vm, std::abs(v[i]) * w[i]);
    }
    return vm;
}


/**
 * @fn DFNORM
 * 
C This function computes the norm of a full N by N matrix,
C stored in the array A, that is consistent with the weighted max-norm
C on vectors, with weights stored in the array W:
C   DFNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
 */
odepack_cpp_real Odepack::DFNORM(int n, odepack_cpp_real *a, odepack_cpp_real *w)
{
#ifndef MATA
#define MATA(i, j) ARRAY2D(a, n, i, j)
#endif
//
    int i, j;
    odepack_cpp_real an, sum;
    an = 0.0;
    for (i = 1; i <= n; ++i) {
        sum = 0.0;
        for (j = 1; j <= n; ++j) {
            sum += std::abs(MATA(i, j)) / ARRAY1D(w, j);
        }
        an = std::max(an, sum * ARRAY1D(w, i));
    }
    return an;
//
#ifdef MATA
#undef MATA
#endif
}


/**
 * @fn DMNORM
 * 
C This function computes the norm of a banded N by N matrix,
C stored in the array A, that is consistent with the weighted max-norm
C on vectors, with weights stored in the array W.
C ML and MU are the lower and upper half-bandwidths of the matrix.
C NRA is the first dimension of the A array, NRA .ge. ML+MU+1.
C In terms of the matrix elements a(i,j), the norm is given by:
C   DBNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
 */
odepack_cpp_real Odepack::DBNORM(int n, odepack_cpp_real *a, int nra, int ml, int mu, odepack_cpp_real *w)
{
#ifndef MATA
#define MATA(i, j) ARRAY2D(a, nra, i, j)
#endif
//
    int i, i1, jlo, jhi, j;
    odepack_cpp_real an, sum;
    an = 0.0;
    for (i = 1; i <= n; ++i) {
        sum = 0.0;
        i1 = i + mu + 1;
        jlo = std::max(i - ml, 1);
        jhi = std::min(i + mu, n);
        for (j = jlo; j <= jhi; ++j) {
            sum += std::abs(MATA(i1 - j, j)) / ARRAY1D(w, j);
        }
        an = std::max(an, sum * ARRAY1D(w, j));
    }
    return an;
//
#ifdef MATA
#undef MATA
#endif
}


/**
 * @fn DSRCMA
 * 
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLSA01, which are used
C internally by one or more ODEPACK solvers.
C
C RSAV = real array of length 240 or more.
C ISAV = integer array of length 46 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
 */
void Odepack::DSRCMA(
    odepack_cpp_real *rsav, int *isav, int job
)
{
    int i;
// DLS001
    odepack_cpp_real *rls = dls1_.rls;
    int *ils = dls1_.ils;
// DLSS01
    odepack_cpp_real *rlsa = dlsa_.rlsa;
    int *ilsa = dlsa_.ilsa;
//
    int lenrls = 218;
    int lenils = 37;
    int lenrla = 22;
    int lenila = 9;
//
    if (job == 2) goto LABEL_100;
    for (i = 1; i <= lenrls; ++i) {
        ARRAY1D(rsav, i) = ARRAY1D(rls, i);
    }
    for (i = 1; i <= lenrla; ++i) {
        ARRAY1D(rsav, lenrls + i) = ARRAY1D(rlsa, i);
    }
//
    for (i = 1; i <= lenils; ++i) {
        ARRAY1D(isav, i) = ARRAY1D(ils, i);
    }
    for (i = 1; i <= lenila; ++i) {
        ARRAY1D(isav, lenils + i) = ARRAY1D(ilsa, i);
    }
//
    return;
//
LABEL_100:
    for (i = 1; i <= lenrls; ++i) {
        ARRAY1D(rls, i) = ARRAY1D(rsav, i);
    }
    for (i = 1; i <= lenrla; ++i) {
        ARRAY1D(rlsa, i) = ARRAY1D(rsav, lenrls + i);
    }
//
    for (i = 1; i <= lenils; ++i) {
        ARRAY1D(ils, i) = ARRAY1D(isav, i);
    }
    for (i = 1; i <= lenila; ++i) {
        ARRAY1D(ilsa, i) = ARRAY1D(isav, lenils + i);
    }
//
    return;
}


/**
 * @fn DRCHEK
 * 
C This routine checks for the presence of a root in the vicinity of
C the current T, in a manner depending on the input flag JOB.  It calls
C Subroutine DROOTS to locate the root as precisely as possible.
C
C In addition to variables described previously, DRCHEK
C uses the following for communication:
C JOB    = integer flag indicating type of call:
C          JOB = 1 means the problem is being initialized, and DRCHEK
C                  is to look for a root at or very near the initial T.
C          JOB = 2 means a continuation call to the solver was just
C                  made, and DRCHEK is to check for a root in the
C                  relevant part of the step last taken.
C          JOB = 3 means a successful step was just taken, and DRCHEK
C                  is to look for a root in the interval of the step.
C G0     = array of length NG, containing the value of g at T = T0.
C          G0 is input for JOB .ge. 2, and output in all cases.
C G1,GX  = arrays of length NG for work space.
C IRT    = completion flag:
C          IRT = 0  means no root was found.
C          IRT = -1 means JOB = 1 and a root was found too near to T.
C          IRT = 1  means a legitimate root was found (JOB = 2 or 3).
C                   On return, T0 is the root location, and Y is the
C                   corresponding solution vector.
C T0     = value of T at one endpoint of interval of interest.  Only
C          roots beyond T0 in the direction of integration are sought.
C          T0 is input if JOB .ge. 2, and output in all cases.
C          T0 is updated by DRCHEK, whether a root is found or not.
C TLAST  = last value of T returned by the solver (input only).
C TOUTC  = copy of TOUT (input only).
C IRFND  = input flag showing whether the last step taken had a root.
C          IRFND = 1 if it did, = 0 if not.
C ITASKC = copy of ITASK (input only).
C NGC    = copy of NG (input only).
 */
void Odepack::DRCHEK(int job, ODEPACK_CONSTRAINT g, int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *g0, odepack_cpp_real *g1, odepack_cpp_real *gx, int *jroot, int &irt, void *user_data)
{
#ifndef YH
#define YH(i, j) ARRAY2D(yh, nyh, i, j)
#endif
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
// DLSR01
    odepack_cpp_real &t0 = dlsr_.t0, &tlast = dlsr_.tlast, &toutc = dlsr_.toutc;
    int &irfnd = dlsr_.irfnd, &itaskc = dlsr_.itaskc, &ngc = dlsr_.ngc, &nge = dlsr_.nge;
//
    int i, iflag, jflag;
    odepack_cpp_real hming, t1, temp1, temp2, x;
    bool zroot;
//
    irt = 0;
    for (i = 1; i <= ngc; ++i) {
        ARRAY1D(jroot, i) = 0;
    }
    hming = (std::abs(tn) + std::abs(h)) * uround * 100.0;
//
    if (job == 1) {
        goto LABEL_100;
    } else if (job == 2) {
        goto LABEL_200;
    } else if (job == 3) {
        goto LABEL_300;
    }
//
// Evaluate g at initial T, and check for zero values. ------------------
LABEL_100:
    t0  = tn;
    (*g)(neq, t0, y, ngc, g0, user_data);
    nge = 1;
    zroot = false;
    for (i = 1; i <= ngc; ++i) {
        if (std::abs(ARRAY1D(g0, i)) <= 0.0) zroot = true;
    }
    if (!zroot) goto LABEL_190;
// g has a zero at T.  Look at g at T + (small increment). --------------
    temp2 = std::max(hming / std::abs(h), 0.1);
    temp1 = temp2 * h;
    t0 += temp1;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) += temp2 * YH(i, 2);
    }
    (*g)(neq, t0, y, ngc, g0, user_data);
    nge++;
    zroot = false;
    for (i = 1; i <= ngc; ++i) {
        if (std::abs(ARRAY1D(g0, i)) <= 0.0) zroot = true;
    }
    if (!zroot) goto LABEL_190;
// g has a zero at T and also close to T.  Take error return. -----------
    irt = -1;
    return;
//
LABEL_190:
    return;
//
//
LABEL_200:
    if (irfnd == 0) goto LABEL_260;
// If a root was found on the previous step, evaluate G0 = g(T0). -------
    DINTDY(t0, 0, yh, nyh, y, iflag);
    (*g)(neq, t0, y, ngc, g0, user_data);
    nge++;
    zroot = false;
    for (i = 1; i <= ngc; ++i) {
        if (std::abs(ARRAY1D(g0, i)) <= 0.0) zroot = true;
    }
    if (!zroot) goto LABEL_260;
// g has a zero at T0.  Look at g at T + (small increment). -------------
    temp1 = std::abs(hming) * (h >= 0 ? 1.0: -1.0);
    t0 += temp1;
    if ((t0 - tn) * h < 0.0) goto LABEL_230;
    temp2 = temp1 / h;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) += temp2 * YH(i, 2);
    }
    goto LABEL_240;
LABEL_230:
    DINTDY(t0, 0, yh, nyh, y, iflag);
LABEL_240:
    (*g)(neq, t0, y, ngc, g0, user_data);
    nge++;
    zroot = false;
    for (i = 1; i <= ngc; ++i) {
        if (std::abs(ARRAY1D(g0, i)) > 0.0) continue;
        ARRAY1D(jroot, i) = 1;
        zroot = true;
    }
    if (!zroot) goto LABEL_260;
// g has a zero at T0 and also close to T0.  Return root. ---------------
    irt = 1;
    return;
// G0 has no zero components.  Proceed to check relevant interval. ------
LABEL_260:
    if (tn == tlast) goto LABEL_390;
//
LABEL_300:
// Set T1 to TN or TOUTC, whichever comes first, and get g at T1. -------
    if (itaskc == 2 || itaskc == 3 || itaskc || 5) goto LABEL_310;
    if ((toutc - tn) * h >= 0.0) goto LABEL_310;
    t1 = toutc;
    if ((t1 - t0) * h <= 0.0) goto LABEL_390;
    DINTDY(t1, 0, yh, nyh, y, iflag);
    goto LABEL_330;
LABEL_310:
    t1 = tn;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = YH(i, 1);
    }
LABEL_330:
    (*g)(neq, t1, y, ngc, g1, user_data);
    nge++;
// Call DROOTS to search for root in interval from T0 to T1. ------------
    jflag = 0;
LABEL_350:
    DROOTS(ngc, hming, jflag, t0, t1, g0, g1, gx, x, jroot);
    if (jflag > 1) goto LABEL_360;
    DINTDY(x, 0, yh, nyh, y, iflag);
    (*g)(neq, x, y, ngc, gx, user_data);
    nge++;
    goto LABEL_350;
LABEL_360:
    t0 = x;
    DCOPY(ngc, gx, 1, g0, 1);
    if (jflag == 4) goto LABEL_390;
// Found a root.  interpolate to X and return. --------------------------
    DINTDY(x, 0, yh, nyh, y, iflag);
    irt = 1;
    return;
//
LABEL_390:
    return;
//
#ifdef YH
#undef YH
#endif
}


/**
 * @fn DROOTS
 * 
C This subroutine finds the leftmost root of a set of arbitrary
C functions gi(x) (i = 1,...,NG) in an interval (X0,X1).  Only roots
C of odd multiplicity (i.e. changes of sign of the gi) are found.
C Here the sign of X1 - X0 is arbitrary, but is constant for a given
C problem, and -leftmost- means nearest to X0.
C The values of the vector-valued function g(x) = (gi, i=1...NG)
C are communicated through the call sequence of DROOTS.
C The method used is the Illinois algorithm.
C
C Reference:
C Kathie L. Hiebert and Lawrence F. Shampine, Implicitly Defined
C Output Points for Solutions of ODEs, Sandia Report SAND80-0180,
C February 1980.
C
C Description of parameters.
C
C NG     = number of functions gi, or the number of components of
C          the vector valued function g(x).  Input only.
C
C HMIN   = resolution parameter in X.  Input only.  When a root is
C          found, it is located only to within an error of HMIN in X.
C          Typically, HMIN should be set to something on the order of
C               100 * UROUND * MAX(ABS(X0),ABS(X1)),
C          where UROUND is the unit roundoff of the machine.
C
C JFLAG  = integer flag for input and output communication.
C
C          On input, set JFLAG = 0 on the first call for the problem,
C          and leave it unchanged until the problem is completed.
C          (The problem is completed when JFLAG .ge. 2 on return.)
C
C          On output, JFLAG has the following values and meanings:
C          JFLAG = 1 means DROOTS needs a value of g(x).  Set GX = g(X)
C                    and call DROOTS again.
C          JFLAG = 2 means a root has been found.  The root is
C                    at X, and GX contains g(X).  (Actually, X is the
C                    rightmost approximation to the root on an interval
C                    (X0,X1) of size HMIN or less.)
C          JFLAG = 3 means X = X1 is a root, with one or more of the gi
C                    being zero at X1 and no sign changes in (X0,X1).
C                    GX contains g(X) on output.
C          JFLAG = 4 means no roots (of odd multiplicity) were
C                    found in (X0,X1) (no sign changes).
C
C X0,X1  = endpoints of the interval where roots are sought.
C          X1 and X0 are input when JFLAG = 0 (first call), and
C          must be left unchanged between calls until the problem is
C          completed.  X0 and X1 must be distinct, but X1 - X0 may be
C          of either sign.  However, the notion of -left- and -right-
C          will be used to mean nearer to X0 or X1, respectively.
C          When JFLAG .ge. 2 on return, X0 and X1 are output, and
C          are the endpoints of the relevant interval.
C
C G0,G1  = arrays of length NG containing the vectors g(X0) and g(X1),
C          respectively.  When JFLAG = 0, G0 and G1 are input and
C          none of the G0(i) should be zero.
C          When JFLAG .ge. 2 on return, G0 and G1 are output.
C
C GX     = array of length NG containing g(X).  GX is input
C          when JFLAG = 1, and output when JFLAG .ge. 2.
C
C X      = independent variable value.  Output only.
C          When JFLAG = 1 on output, X is the point at which g(x)
C          is to be evaluated and loaded into GX.
C          When JFLAG = 2 or 3, X is the root.
C          When JFLAG = 4, X is the right endpoint of the interval, X1.
C
C JROOT  = integer array of length NG.  Output only.
C          When JFLAG = 2 or 3, JROOT indicates which components
C          of g(x) have a root at X.  JROOT(i) is 1 if the i-th
C          component has a root, and JROOT(i) = 0 otherwise.
 */
void Odepack::DROOTS(
    int ng, odepack_cpp_real hmin, int &jflag, odepack_cpp_real &x0, odepack_cpp_real &x1, odepack_cpp_real *g0, odepack_cpp_real *g1, odepack_cpp_real *gx, odepack_cpp_real &x, int *jroot
)
{
// DLSR01
    odepack_cpp_real &alpha = dlsr_.alpha, &x2 = dlsr_.x2;
    int &imax = dlsr_.imax, &last = dlsr_.last;
//
    int i, imxold, nxlast;
    odepack_cpp_real t2, tmax, fracint, fracsub;
    bool zroot, sgnchg, xroot;
    const odepack_cpp_real zero  = 0.0;
    const odepack_cpp_real half  = 0.0;
    const odepack_cpp_real tenth = 0.1;
    const odepack_cpp_real five  = 5.0;
//
    if (jflag == 1) goto LABEL_200;
// JFLAG .ne. 1.  Check for change in sign of g or zero at X1. ----------
    imax = 0;
    tmax = zero;
    zroot = false;
    for (i = 1; i <= ng; ++i) {
        if (std::abs(ARRAY1D(g1, i)) > zero) goto LABEL_110;
        zroot = true;
        continue;
// At this point, G0(i) has been checked and cannot be zero. ------------
LABEL_110:
        if (SIGN(ARRAY1D(g0, i)) == SIGN(ARRAY1D(g1, i))) continue;
        t2 = std::abs(ARRAY1D(g1, i) / (ARRAY1D(g1, i) - ARRAY1D(g0, i)));
        if (t2 <= tmax) continue;
        tmax = t2;
        imax = i;
    }
    if (imax > 0) goto LABEL_130;
    sgnchg = false;
    goto LABEL_140;
LABEL_130:
    sgnchg = true;
LABEL_140:
    if (!sgnchg) goto LABEL_400;
// There is a sign change.  Find the first root in the interval. --------
    xroot = false;
    nxlast = 0;
    last = 1;
//
// Repeat until the first root in the interval is found.  Loop point. ---
LABEL_150:
    if (xroot) goto LABEL_300;
    if (nxlast == last) goto LABEL_160;
    alpha = 1.0;
    goto LABEL_180;
LABEL_160:
    if (last == 0) goto LABEL_170;
    alpha *= 0.5;
    goto LABEL_180;
LABEL_170:
    alpha *= 2.0;
LABEL_180:
    x2 = x1 - (x1 - x0) * ARRAY1D(g1, imax) / (ARRAY1D(g1, imax) - alpha * ARRAY1D(g0, imax));
// If X2 is too close to X0 or X1, adjust it inward, by a fractional ----
// distance that is between 0.1 and 0.5. --------------------------------
    if (std::abs(x2 - x0) < half * hmin) {
        fracint = std::abs(x1 - x0) / hmin;
        fracsub = tenth;
        if (fracint <= five) fracsub = half / fracint;
        x2 = x0 + fracsub * (x1 - x0);
    }
    if (std::abs(x1 - x2) < half * hmin) {
        fracint = std::abs(x1 - x0) / hmin;
        fracsub = tenth;
        if (fracint <= five) fracsub = half / fracint;
        x2 = x1 - fracsub * (x1 - x0);
    }
    jflag = 1;
    x = x2;
// Return to the calling routine to get a value of GX = g(X). -----------
    return;
// Check to see in which interval g changes sign. -----------------------
LABEL_200:
    imxold = imax;
    imax = 0;
    tmax = zero;
    zroot = false;
    for (i = 1; i <= ng; ++i) {
        if (std::abs(ARRAY1D(gx, i)) > zero) goto LABEL_210;
        zroot = true;
        continue;
LABEL_210:
        if (SIGN(ARRAY1D(g0, i)) == SIGN(ARRAY1D(gx, i))) continue;
        t2 = std::abs(ARRAY1D(gx, i) / (ARRAY1D(gx, i) - ARRAY1D(g0, i)));
        if (t2 <= tmax) continue;
        tmax = t2;
        imax = i;
    }
    if (imax > 0) goto LABEL_230;
    sgnchg = false;
    imax = imxold;
    goto LABEL_240;
LABEL_230:
    sgnchg = true;
LABEL_240:
    nxlast = last;
    if (!sgnchg) goto LABEL_250;
// Sign change between X0 and X2, so replace X1 with X2. ----------------
    x1 = x2;
    DCOPY(ng, gx, 1, g1, 1);
    last = 1;
    xroot = false;
    goto LABEL_270;
LABEL_250:
    if (!zroot) goto LABEL_260;
// Zero value at X2 and no sign change in (X0,X2), so X2 is a root. -----
    x1 = x2;
    DCOPY(ng, gx, 1, g1, 1);
    xroot = true;
    goto LABEL_270;
// No sign change between X0 and X2.  Replace X0 with X2. ---------------
LABEL_260:
    DCOPY(ng, gx, 1, g0, 1);
    x0 = x2;
    last = 0;
    xroot = false;
LABEL_270:
    if (std::abs(x1 - x0) <= hmin) xroot = true;
    goto LABEL_150;
//
// Return with X1 as the root.  Set JROOT.  Set X = X1 and GX = G1. -----
LABEL_300:
    jflag = 2;
    x = x1;
    DCOPY(ng, g1, 1, gx, 1);
    for (i = 1; i <= ng; ++i) {
        ARRAY1D(jroot, i) = 0;
        if (std::abs(ARRAY1D(g1, i)) > zero) goto LABEL_310;
        ARRAY1D(jroot, i) = 1;
        continue;
LABEL_310:
        if (SIGN(ARRAY1D(g0, i)) != SIGN(ARRAY1D(g1, i))) ARRAY1D(jroot, i) = 1;
    }
    return;
//
// No sign change in the interval.  Check for zero at right endpoint. ---
LABEL_400:
    if (!zroot) goto LABEL_420;
//
// Zero value at X1 and no sign change in (X0,X1).  Return JFLAG = 3. ---
    x = x1;
    DCOPY(ng, g1, 1, gx, 1);
    for (i = 1; i <= ng; ++i) {
        ARRAY1D(jroot, i) = 0;
        if (std::abs(ARRAY1D(g1, i)) <= zero) ARRAY1D(jroot, i) = 1;
    }
    jflag = 3;
    return;
//
// No sign changes in this interval.  Set X = X1, return JFLAG = 4. -----
LABEL_420: 
    DCOPY(ng, g1, 1, gx, 1);
    x = x1;
    jflag = 4;
    return;
}


/**
 * @fn DSRCAR
 * 
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLSA01, DLSR01, which are used
C internally by one or more ODEPACK solvers.
C
C RSAV = real array of length 245 or more.
C ISAV = integer array of length 55 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
 */
void Odepack::DSRCAR(
    odepack_cpp_real *rsav, int *isav, int job
)
{
    int i, ioff;
// DLS001
    odepack_cpp_real *rls = dls1_.rls;
    int *ils = dls1_.ils;
// DLSA01
    odepack_cpp_real *rlsa = dlsa_.rlsa;
    int *ilsa = dlsa_.ilsa;
// DLSR01
    odepack_cpp_real *rlsr = dlsr_.rlsr;
    int *ilsr = dlsr_.ilsr;
//
    int lenrls = 218;
    int lenils = 37;
    int lenrla = 22;
    int lenila = 9;
    int lenrlr = 5;
    int lenilr = 9;
//
    if (job == 2) goto LABEL_100;
    for (i = 1; i <= lenrls; ++i) {
        ARRAY1D(rsav, i) = ARRAY1D(rls, i);
    }
    for (i = 1; i <= lenrla; ++i) {
        ARRAY1D(rsav, lenrls + i) = ARRAY1D(rlsa, i);
    }
    ioff = lenrls + lenrla;
    for (i = 1; i <= lenrlr; ++i) {
        ARRAY1D(rsav, ioff + i) = ARRAY1D(rlsr, i);
    }
//
    for (i = 1; i <= lenils; ++i) {
        ARRAY1D(isav, i) = ARRAY1D(ils, i);
    }
    for (i = 1; i <= lenila; ++i) {
        ARRAY1D(isav, lenils + i) = ARRAY1D(ilsa, i);
    }
    ioff = lenils + lenila;
    for (i = 1; i <= lenilr; ++i) {
        ARRAY1D(isav, ioff + i) = ARRAY1D(ilsr, i);
    }
//
    return;
//
LABEL_100:
    for (i = 1; i <= lenrls; ++i) {
        ARRAY1D(rls, i) = ARRAY1D(rsav, i);
    }
    for (i = 1; i <= lenrla; ++i) {
        ARRAY1D(rlsa, i) = ARRAY1D(rsav, lenrls + i);
    }
    ioff = lenrls + lenrla;
    for (i = 1; i <= lenrlr; ++i) {
        ARRAY1D(rlsr, i) = ARRAY1D(rsav, ioff + i);
    }
//
    for (i = 1; i <= lenils; ++i) {
        ARRAY1D(ils, i) = ARRAY1D(isav, i);
    }
    for (i = 1; i <= lenila; ++i) {
        ARRAY1D(ilsa, i) = ARRAY1D(isav, lenils + i);
    }
    ioff = lenils + lenila;
    for (i = 1; i <= lenilr; ++i) {
        ARRAY1D(ilsr, i) = ARRAY1D(isav, ioff + i);
    }
//
    return;
}

/**
C-----------------------------------------------------------------------
C DSTODPK performs one step of the integration of an initial value
C problem for a system of Ordinary Differential Equations.
C-----------------------------------------------------------------------
C The following changes were made to generate Subroutine DSTODPK
C from Subroutine DSTODE:
C 1. The array SAVX was added to the call sequence.
C 2. PJAC and SLVS were replaced by PSOL in the call sequence.
C 3. The Common block /DLPK01/ was added for communication.
C 4. The test constant EPCON is loaded into Common below statement
C    numbers 125 and 155, and used below statement 400.
C 5. The Newton iteration counter MNEWT is set below 220 and 400.
C 6. The call to PJAC was replaced with a call to DPKSET (fixed name),
C    with a longer call sequence, called depending on JACFLG.
C 7. The corrector residual is stored in SAVX (not Y) at 360,
C    and the solution vector is in SAVX in the 380 loop.
C 8. SLVS was renamed DSOLPK and includes NEQ, SAVX, EWT, F, and JAC.
C    SAVX was added because DSOLPK now needs Y and SAVF undisturbed.
C 9. The nonlinear convergence failure count NCFN is set at 430.
C-----------------------------------------------------------------------
C Note: DSTODPK is independent of the value of the iteration method
C indicator MITER, when this is .ne. 0, and hence is independent
C of the type of chord method used, or the Jacobian structure.
C Communication with DSTODPK is done with the following variables:
C
C NEQ    = integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls to F and JAC.
C Y      = an array of length .ge. N used as the Y argument in
C          all calls to F and JAC.
C YH     = an NYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C          j-th derivative of y(i), scaled by H**j/factorial(j)
C          (j = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C NYH    = a constant integer .ge. N, the first dimension of YH.
C YH1    = a one-dimensional array occupying the same space as YH.
C EWT    = an array of length N containing multiplicative weights
C          for local error measurements.  Local errors in y(i) are
C          compared to 1.0/EWT(i) in various error tests.
C SAVF   = an array of working storage, of length N.
C          Also used for input of YH(*,MAXORD+2) when JSTART = -1
C          and MAXORD .lt. the current order NQ.
C SAVX   = an array of working storage, of length N.
C ACOR   = a work array of length N, used for the accumulated
C          corrections.  On a successful return, ACOR(i) contains
C          the estimated one-step local error in y(i).
C WM,IWM = real and integer work arrays associated with matrix
C          operations in chord iteration (MITER .ne. 0).
C CCMAX  = maximum relative change in H*EL0 before DPKSET is called.
C H      = the step size to be attempted on the next step.
C          H is altered by the error control algorithm during the
C          problem.  H can be either positive or negative, but its
C          sign must remain constant throughout the problem.
C HMIN   = the minimum absolute value of the step size H to be used.
C HMXI   = inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C          HMIN and HMXI may be changed at any time, but will not
C          take effect until the next change of H is considered.
C TN     = the independent variable. TN is updated on each step taken.
C JSTART = an integer used for input only, with the following
C          values and meanings:
C               0  perform the first step.
C           .gt.0  take a new step continuing from the last.
C              -1  take the next step with a new value of H, MAXORD,
C                    N, METH, MITER, and/or matrix parameters.
C              -2  take the next step with a new value of H,
C                    but with other inputs unchanged.
C          On return, JSTART is set to 1 to facilitate continuation.
C KFLAG  = a completion code with the following meanings:
C               0  the step was succesful.
C              -1  the requested error could not be achieved.
C              -2  corrector convergence could not be achieved.
C              -3  fatal error in DPKSET or DSOLPK.
C          A return with KFLAG = -1 or -2 means either
C          ABS(H) = HMIN or 10 consecutive failures occurred.
C          On a return with KFLAG negative, the values of TN and
C          the YH array are as of the beginning of the last
C          step, and H is the last step size attempted.
C MAXORD = the maximum order of integration method to be allowed.
C MAXCOR = the maximum number of corrector iterations allowed.
C MSBP   = maximum number of steps between DPKSET calls (MITER .gt. 0).
C MXNCF  = maximum number of convergence failures allowed.
C METH/MITER = the method flags.  See description in driver.
C N      = the number of first-order differential equations.
C-----------------------------------------------------------------------
 */
void Odepack::DSTODPK(int neq, odepack_cpp_real *y,odepack_cpp_real *yh, int nyh, odepack_cpp_real *yh1, odepack_cpp_real *ewt, odepack_cpp_real *savf,
    odepack_cpp_real *savx, odepack_cpp_real *acor, odepack_cpp_real *wm, int *iwm, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN3 jac, 
    ODEPACK_PSOL psol, void *user_data)
{
#ifndef YH
#define YH(i, j) ARRAY2D(yh, nyh, i, j)
#endif
#ifndef ELCO
#define ELCO(i, j) ARRAY2D(elco, 13, i, j)
#endif
#ifndef TESCO
#define TESCO(i, j) ARRAY2D(tesco, 3, i, j)
#endif
// DLS001
    odepack_cpp_real &conit = dls1_.conit, &crate = dls1_.crate, *el = dls1_.el, *elco = dls1_.elco,
        &hold = dls1_.hold, &rmax = dls1_.rmax, *tesco = dls1_.tesco,
        &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &ialth = dls1_.ialth, &ipup = dls1_.ipup, &lmax = dls1_.lmax, &meo = dls1_.meo, &nqnyh = dls1_.nqnyh, &nslp = dls1_.nslp,
        &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
// DLPK01
    odepack_cpp_real &delt = dlpk_.delt, &epcon = dlpk_.epcon, &sqrtn = dlpk_.sqrtn, &rsqrtn = dlpk_.rsqrtn;
    int &jpre = dlpk_.jpre, &jacflg = dlpk_.jacflg, &locwp = dlpk_.locwp, &lociwp = dlpk_.lociwp, &lsavx = dlpk_.lsavx, &kmp = dlpk_.kmp, &maxl = dlpk_.maxl, &mnewt = dlpk_.mnewt,
        &nni = dlpk_.nni, &nli = dlpk_.nli, &nps = dlpk_.nps, &ncfn = dlpk_.ncfn, &ncfl = dlpk_.ncfl;
//
    int i, i1, iredo, iret, j, jb, m, ncf, newq;
    odepack_cpp_real dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup,
        r, rh, rhdn, rhsm, rhup, told;
//
    kflag = 0;
    told = tn;
    ncf = 0;
    ierpj = 0;
    iersl = 0;
    jcur = 0;
    icf = 0;
    delp = 0.0;
    if (jstart  >  0) goto LABEL_200;
    if (jstart == -1) goto LABEL_100;
    if (jstart == -2) goto LABEL_160;
//-----------------------------------------------------------------------
// On the first call, the order is set to 1, and other variables are
// initialized.  RMAX is the maximum ratio by which H can be increased
// in a single step.  It is initially 1.E4 to compensate for the small
// initial H, but then is normally equal to 10.  If a failure
// occurs (in corrector convergence or error test), RMAX is set at 2
// for the next increase.
//-----------------------------------------------------------------------
    lmax = maxord + 1;
    nq = 1;
    l = 2;
    ialth = 2;
    rmax = 10000.0;
    rc = 0.0;
    el0 = 1.0;
    crate = 0.7;
    hold = h;
    meo = meth;
    nslp = 0;
    ipup = miter;
    iret = 3;
    goto LABEL_140;
//-----------------------------------------------------------------------
// The following block handles preliminaries needed when JSTART = -1.
// IPUP is set to MITER to force a matrix update.
// If an order increase is about to be considered (IALTH = 1),
// IALTH is reset to 2 to postpone consideration one more step.
// If the caller has changed METH, DCFODE is called to reset
// the coefficients of the method.
// If the caller has changed MAXORD to a value less than the current
// order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
// If H is to be changed, YH must be rescaled.
// If H or METH is being changed, IALTH is reset to L = NQ + 1
// to prevent further changes in H for that many steps.
//-----------------------------------------------------------------------
LABEL_100:
    ipup = miter;
    lmax = maxord + 1;
    if (ialth == 1) ialth = 2;
    if (meth == meo) goto LABEL_110;
    DCFODE (meth, elco, tesco);
    meo = meth;
    if (nq > maxord) goto LABEL_120;
    ialth = l;
    iret = 1;
    goto LABEL_150;
LABEL_110:
    if (nq <= maxord) goto LABEL_160;
LABEL_120:
    nq = maxord;
    l = lmax;
    for (i = 1; i <= l; ++i) {
        ARRAY1D(el, i) = ELCO(i, nq);
    }
    nqnyh = nq * nyh;
    rc = rc * ARRAY1D(el, 1) / el0;
    el0 = ARRAY1D(el, 1);
    conit = 0.5 / static_cast<odepack_cpp_real>(nq + 2);
    epcon = conit * TESCO(2, nq);
    ddn = DVNORM(n, savf, ewt) / TESCO(1, l);
    exdn = 1.0 / static_cast<odepack_cpp_real>(l);
    rhdn = 1.0 / (1.3 * std::pow(ddn, exdn) + 0.0000013);
    rh = std::min(rhdn, 1.0);
    iredo = 3;
    if (h == hold) goto LABEL_170;
    rh = std::min(rh,std::abs(h / hold));
    h = hold;
    goto LABEL_175;
//-----------------------------------------------------------------------
// DCFODE is called to get all the integration coefficients for the
// current METH.  Then the EL vector and related constants are reset
// whenever the order NQ is changed, or at the start of the problem.
//-----------------------------------------------------------------------
LABEL_140:
    DCFODE(meth, elco, tesco);
LABEL_150:
    for (i = 1; i <= l; ++i) {
        ARRAY1D(el, i) = ELCO(i, nq);
    }
    nqnyh = nq * nyh;
    rc = rc * ARRAY1D(el, 1) / el0;
    el0 = ARRAY1D(el, 1);
    conit = 0.5 / static_cast<odepack_cpp_real>(nq + 2);
    epcon = conit * TESCO(2, nq);
    if (iret == 1) {
        goto LABEL_160;
    } else if (iret == 2) {
        goto LABEL_170;
    } else if (iret == 3) {
        goto LABEL_200;
    }
//-----------------------------------------------------------------------
// If H is being changed, the H ratio RH is checked against
// RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
// L = NQ + 1 to prevent a change of H for that many steps, unless
// forced by a convergence or error test failure.
//-----------------------------------------------------------------------
LABEL_160:
    if (h == hold) goto LABEL_200;
    rh = h / hold;
    h = hold;
    iredo = 3;
    goto LABEL_175;
LABEL_170:
    rh = std::max(rh, hmin / std::abs(h));
LABEL_175:
    rh = std::min(rh,rmax);
    rh = rh / std::max(1.0, std::abs(h) * hmxi * rh);
    r  = 1.0;
    for (j = 2; j <= l; ++j) {
        r = r * rh;
        for (i = 1; i <= n; ++i) {
            YH(i, j) = YH(i, j) * r;
        }
    }
    h *= rh;
    rc *= rh;
    ialth = l;
    if (iredo == 0) goto LABEL_690;
//-----------------------------------------------------------------------
// This section computes the predicted values by effectively
// multiplying the YH array by the Pascal triangle matrix.
// The flag IPUP is set according to whether matrix data is involved
// (JACFLG .ne. 0) or not (JACFLG = 0), to trigger a call to DPKSET.
// IPUP is set to MITER when RC differs from 1 by more than CCMAX,
// and at least every MSBP steps, when JACFLG = 1.
// RC is the ratio of new to old values of the coefficient  H*EL(1).
//-----------------------------------------------------------------------
LABEL_200:
    if (jacflg != 0) goto LABEL_202;
    ipup = 0;
    crate = 0.7;
    goto LABEL_205;
LABEL_202:
    if (std::abs(rc-1.0) > ccmax) ipup = miter;
    if (nst >= nslp + msbp) ipup = miter;
LABEL_205:
    tn = tn + h;
    i1 = nqnyh + 1;
    for (jb = 1; jb <= nq; ++jb) {
        i1 = i1 - nyh;
        for (i = i1; i <= nqnyh; ++i) {
            ARRAY1D(yh1, i) = ARRAY1D(yh, i) + ARRAY1D(yh1, i + nyh);
        }
    }
//-----------------------------------------------------------------------
// Up to MAXCOR corrector iterations are taken.  A convergence test is
// made on the RMS-norm of each correction, weighted by the error
// weight vector EWT.  The sum of the corrections is accumulated in the
// vector ACOR(i).  The YH array is not altered in the corrector loop.
//-----------------------------------------------------------------------
LABEL_220:  
    m = 0;
    mnewt = 0;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = YH(i, 1);
    }
    (*f)(neq, tn, y, savf, user_data);
    nfe++;
    if (ipup <= 0) goto LABEL_250;
//-----------------------------------------------------------------------
// If indicated, DPKSET is called to update any matrix data needed,
// before starting the corrector iteration.
// IPUP is set to 0 as an indicator that this has been done.
//-----------------------------------------------------------------------
    DPKSET(neq, y, yh1, ewt, acor, savf, wm, iwm, f, jac, user_data);
    ipup = 0;
    rc = 1.0;
    nslp = nst;
    crate = 0.7;
    if (ierpj != 0) goto LABEL_430;
LABEL_250:
    for (i = 1; i <= n; ++i) {
        ARRAY1D(acor, i) = 0.0;
    }
LABEL_270:
    if (miter != 0) goto LABEL_350;
//-----------------------------------------------------------------------
// In the case of functional iteration, update Y directly from
// the result of the last function evaluation.
//-----------------------------------------------------------------------
    for (i = 1; i <= n; ++i) {
        ARRAY1D(savf, i) = h * ARRAY1D(savf, i) - YH(i, 2);
        ARRAY1D(y, i) = ARRAY1D(savf, i) - ARRAY1D(acor, i);
    }
    del = DVNORM(n, y, ewt);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = YH(i, 1) + ARRAY1D(el, 1) * ARRAY1D(savf, i);
        ARRAY1D(acor, i) = ARRAY1D(savf, i);
    }
    goto LABEL_400;
//-----------------------------------------------------------------------
// In the case of the chord method, compute the corrector error,
// and solve the linear system with that as right-hand side and
// P as coefficient matrix.
//-----------------------------------------------------------------------
LABEL_350:  
    for (i = 1; i <= n; ++i) {
        ARRAY1D(savx, i) = h * ARRAY1D(savf, i) - (YH(i, 2) + ARRAY1D(acor, i));
    }
    DSOLPK(neq, y, savf, savx, ewt, wm, iwm, f, psol, user_data);
    if (iersl < 0) goto LABEL_430;
    if (iersl > 0) goto LABEL_410;
    del = DVNORM (n, savx, ewt);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(acor, i) = ARRAY1D(acor, i) + ARRAY1D(savx, i);
        ARRAY1D(y, i) = YH(i, 1) + ARRAY1D(el, 1) * ARRAY1D(acor, i);
    }
//-----------------------------------------------------------------------
// Test for convergence.  If M .gt. 0, an estimate of the convergence
// rate constant is stored in CRATE, and this is used in the test.
//-----------------------------------------------------------------------
LABEL_400:  
    if (m != 0) crate = std::max(0.2 * crate, del / delp);
    dcon = del * std::min(1.0, 1.5 * crate) / epcon;
    if (dcon <= 1.0) goto LABEL_450;
    m++;
    if (m == maxcor) goto LABEL_410;
    if (m >= 2 && del > 2.0 * delp) goto LABEL_410;
    mnewt = m;
    delp = del;
    (*f)(neq, tn, y, savf, user_data);
    nfe++;
    goto LABEL_270;
//-----------------------------------------------------------------------
// The corrector iteration failed to converge.
// If MITER .ne. 0 and the Jacobian is out of date, DPKSET is called for
// the next try.  Otherwise the YH array is retracted to its values
// before prediction, and H is reduced, if possible.  If H cannot be
// reduced or MXNCF failures have occurred, exit with KFLAG = -2.
//-----------------------------------------------------------------------
LABEL_410:  
    if (miter == 0 || jcur == 1 || jacflg == 0) goto LABEL_430;
    icf = 1;
    ipup = miter;
    goto LABEL_220;
LABEL_430:  
    icf = 2;
    ncf++;
    ncfn++;
    rmax = 2.0;
    tn = told;
    i1 = nqnyh + 1;
    for (jb = 1; jb <= nq; ++jb) {
        i1 = i1 - nyh;
        for (i = i1; i <= nqnyh; ++i) {
            ARRAY1D(yh, i) = ARRAY1D(yh, i) - ARRAY1D(yh, i + nyh);
        }
    }
    if (ierpj < 0 || iersl < 0) goto LABEL_680;
    if (std::abs(h) <= hmin * 1.00001) goto LABEL_670;
    if (ncf == mxncf) goto LABEL_670;
    rh = 0.5;
    ipup = miter;
    iredo = 1;
    goto LABEL_170;
//-----------------------------------------------------------------------
// The corrector has converged.  JCUR is set to 0
// to signal that the Jacobian involved may need updating later.
// The local error test is made and control passes to statement 500
// if it fails.
//-----------------------------------------------------------------------
LABEL_450:  
    jcur = 0;
    if (m == 0) dsm = del / TESCO(2, nq);
    if (m  > 0) dsm = DVNORM(n, acor, ewt) / TESCO(2, nq);
    if (dsm > 1.0) goto LABEL_500;
//-----------------------------------------------------------------------
// After a successful step, update the YH array.
// Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
// If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
// use in a possible order increase on the next step.
// If a change in H is considered, an increase or decrease in order
// by one is considered also.  A change in H is made only if it is by a
// factor of at least 1.1.  If not, IALTH is set to 3 to prevent
// testing for that many steps.
//-----------------------------------------------------------------------
    kflag = 0;
    iredo = 0;
    nst++;
    hu = h;
    nqu = nq;
    for (j = 1; j <= l; ++j) {
        for (i = 1; i <= n; ++i) {
            YH(i, j) = YH(i, j) + ARRAY1D(el, j) * ARRAY1D(acor, i);
        }
    }
    ialth--;
    if (ialth == 0) goto LABEL_520;
    if (ialth  > 1) goto LABEL_700;
    if (l == lmax) goto LABEL_700;
    for (i = 1; i <= n; ++i) {
        YH(i, lmax) = ARRAY1D(acor, i);
    }
    goto LABEL_700;
//-----------------------------------------------------------------------
// The error test failed.  KFLAG keeps track of multiple failures.
// Restore TN and the YH array to their previous values, and prepare
// to try the step again.  Compute the optimum step size for this or
// one lower order.  After 2 or more failures, H is forced to decrease
// by a factor of 0.2 or less.
//-----------------------------------------------------------------------
LABEL_500:  
    kflag--;
    tn = told;
    i1 = nqnyh + 1;
    for (jb = 1; jb <= nq; ++jb) {
        i1 = i1 - nyh;
        for (i = i1; i <= nqnyh; ++i) {
            ARRAY1D(yh, i) = ARRAY1D(yh, i) - ARRAY1D(yh, i + nyh);
        }
    }
    rmax = 2.0;
    if (std::abs(h) <= hmin * 1.00001) goto LABEL_660;
    if (kflag <= -3) goto LABEL_640;
    iredo = 2;
    rhup = 0.0;
    goto LABEL_540;
//-----------------------------------------------------------------------
// Regardless of the success or failure of the step, factors
// RHDN, RHSM, and RHUP are computed, by which H could be multiplied
// at order NQ - 1, order NQ, or order NQ + 1, respectively.
// In the case of failure, RHUP = 0.0 to avoid an order increase.
// the largest of these is determined and the new order chosen
// accordingly.  If the order is to be increased, we compute one
// additional scaled derivative.
//-----------------------------------------------------------------------
LABEL_520:  
    rhup = 0.0;
    if (l == lmax) goto LABEL_540;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(savf, i) = ARRAY1D(acor, i) - YH(i, lmax);
    }
    dup = DVNORM(n, savf, ewt) / TESCO(3, nq);
    exup = 1.0 / (static_cast<odepack_cpp_real>(l + 1));
    rhup = 1.0 / (1.4 * std::pow(dup, exup) + 0.0000014);
LABEL_540:
    exsm = 1.0 / static_cast<odepack_cpp_real>(l);
    rhsm = 1.0 / (1.2 * std::pow(dsm, exsm) + 0.0000012);
    rhdn = 0.0;
    if (nq == 1) goto LABEL_560;
    ddn  = DVNORM(n, &YH(1, l), ewt ) /TESCO(1, nq);
    exdn = 1.0 / static_cast<odepack_cpp_real>(nq);
    rhdn = 1.0 / (1.3 * std::pow(ddn, exdn) + 0.0000013);
LABEL_560:
    if (rhsm >= rhup) goto LABEL_570;
    if (rhup  > rhdn) goto LABEL_590;
    goto LABEL_580;
LABEL_570:
    if (rhsm < rhdn) goto LABEL_580;
    newq = nq;
    rh  = rhsm;
    goto LABEL_620;
LABEL_580:
    newq = nq - 1;
    rh = rhdn;
    if (kflag < 0 && rh > 1.0) rh = 1.0;
    goto LABEL_620;
LABEL_590:
    newq = l;
    rh   = rhup;
    if (rh < 1.1) goto LABEL_610;
    r = ARRAY1D(el, l) / static_cast<odepack_cpp_real>(l);
    for (i = 1; i <= n; ++i) {
        YH(i, newq + 1) = ARRAY1D(acor, i) * r;
    }
    goto LABEL_630;
LABEL_610:  
    ialth = 3;
    goto LABEL_700;
LABEL_620:
    if ((kflag == 0) && (rh < 1.1)) goto LABEL_610;
    if (kflag <= -2) rh = std::min(rh, 0.2);
//-----------------------------------------------------------------------
// If there is a change of order, reset NQ, L, and the coefficients.
// In any case H is reset according to RH and the YH array is rescaled.
// Then exit from 690 if the step was OK, or redo the step otherwise.
//-----------------------------------------------------------------------
    if (newq == nq) goto LABEL_170;
LABEL_630:
    nq = newq;
    l = nq + 1;
    iret = 2;
    goto LABEL_150;
//-----------------------------------------------------------------------
// Control reaches this section if 3 or more failures have occured.
// If 10 failures have occurred, exit with KFLAG = -1.
// It is assumed that the derivatives that have accumulated in the
// YH array have errors of the wrong order.  Hence the first
// derivative is recomputed, and the order is set to 1.  Then
// H is reduced by a factor of 10, and the step is retried,
// until it succeeds or H reaches HMIN.
//-----------------------------------------------------------------------
LABEL_640:  
    if (kflag == -10) goto LABEL_660;
    rh = 0.1;
    rh = std::max(hmin/std::abs(h), rh);
    h *= rh;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = YH(i, 1);
    }
    (*f)(neq, tn, y, savf, user_data);
    nfe = nfe + 1;
    for (i = 1; i <= n; ++i) {
        YH(i, 2) = h * ARRAY1D(savf, i);
    }
    ipup = miter;
    ialth = 5;
    if (nq == 1) goto LABEL_200;
    nq = 1;
    l = 2;
    iret = 3;
    goto LABEL_150;
//-----------------------------------------------------------------------
// All returns are made through this section.  H is saved in HOLD
// to allow the caller to change H on the next step.
//-----------------------------------------------------------------------
LABEL_660:
    kflag = -1;
    goto LABEL_720;
LABEL_670:
    kflag = -2;
    goto LABEL_720;
LABEL_680:
    kflag = -3;
    goto LABEL_720;
LABEL_690:
    rmax = 10.0;
LABEL_700:
    r = 1.0 / TESCO(2, nqu);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(acor, i) *= r;
    }
LABEL_720:  
    hold = h;
    jstart = 1;
    return;
//
#ifdef YH
#undef YH
#endif
#ifdef ELCO
#undef ELCO
#endif
#ifdef TESCO
#undef TESCO
#endif
}


/**
C DPKSET is called by DSTODPK to interface with the user-supplied
C routine JAC, to compute and process relevant parts of
C the matrix P = I - H*EL(1)*J , where J is the Jacobian df/dy,
C as need for preconditioning matrix operations later.
C
C In addition to variables described previously, communication
C with DPKSET uses the following:
C Y     = array containing predicted values on entry.
C YSV   = array containing predicted y, to be saved (YH1 in DSTODPK).
C FTEM  = work array of length N (ACOR in DSTODPK).
C SAVF  = array containing f evaluated at predicted y.
C WM    = real work space for matrices.
C         Space for preconditioning data starts at WM(LOCWP).
C IWM   = integer work space.
C         Space for preconditioning data starts at IWM(LOCIWP).
C IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
C         JAC returned an error flag.
C JCUR  = output flag = 1 to indicate that the Jacobian matrix
C         (or approximation) is now current.
C This routine also uses Common variables EL0, H, TN, IERPJ, JCUR, NJE.
 */
void Odepack::DPKSET(int neq, odepack_cpp_real *y, odepack_cpp_real *ysv, odepack_cpp_real *ewt, odepack_cpp_real *ftem, odepack_cpp_real *savf, odepack_cpp_real *wm, int *iwm,
    ODEPACK_FUNCTION f, ODEPACK_JACOBIAN3 jac, void *user_data)
{
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
// DLPK01
    odepack_cpp_real &delt = dlpk_.delt, &epcon = dlpk_.epcon, &sqrtn = dlpk_.sqrtn, &rsqrtn = dlpk_.rsqrtn;
    int &jpre = dlpk_.jpre, &jacflg = dlpk_.jacflg, &locwp = dlpk_.locwp, &lociwp = dlpk_.lociwp, &lsavx = dlpk_.lsavx, &kmp = dlpk_.kmp, &maxl = dlpk_.maxl, &mnewt = dlpk_.mnewt,
        &nni = dlpk_.nni, &nli = dlpk_.nli, &nps = dlpk_.nps, &ncfn = dlpk_.ncfn, &ncfl = dlpk_.ncfl;
//
    int ier;
    odepack_cpp_real hl0;
//
    ierpj = 0;
    jcur = 1;
    hl0 = el0 * h;
    (*jac)(f, neq, tn, y, ysv, ewt, savf, ftem, hl0, &ARRAY1D(wm, locwp), 
        &ARRAY1D(iwm, lociwp), ier, user_data);
    nje++;
    if (ier == 0) return;
    ierpj = 1;
    return;
}


/**
 * @fn DSOLPK
 * 
C This routine interfaces to one of DSPIOM, DSPIGMR, DPCG, DPCGS, or
C DUSOL, for the solution of the linear system arising from a Newton
C iteration.  It is called if MITER .ne. 0.
C In addition to variables described elsewhere,
C communication with DSOLPK uses the following variables:
C WM    = real work space containing data for the algorithm
C         (Krylov basis vectors, Hessenberg matrix, etc.)
C IWM   = integer work space containing data for the algorithm
C X     = the right-hand side vector on input, and the solution vector
C         on output, of length N.
C IERSL = output flag (in Common):
C         IERSL =  0 means no trouble occurred.
C         IERSL =  1 means the iterative method failed to converge.
C                    If the preconditioner is out of date, the step
C                    is repeated with a new preconditioner.
C                    Otherwise, the stepsize is reduced (forcing a
C                    new evaluation of the preconditioner) and the
C                    step is repeated.
C         IERSL = -1 means there was a nonrecoverable error in the
C                    iterative solver, and an error exit occurs.
C This routine also uses the Common variables TN, EL0, H, N, MITER,
C DELT, EPCON, SQRTN, RSQRTN, MAXL, KMP, MNEWT, NNI, NLI, NPS, NCFL,
C LOCWP, LOCIWP.
 */
void Odepack::DSOLPK(int neq, odepack_cpp_real *y, odepack_cpp_real *savf, odepack_cpp_real *x, odepack_cpp_real *ewt, odepack_cpp_real *wm, int *iwm,
        ODEPACK_FUNCTION f, ODEPACK_PSOL psol, void *user_data)
{
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
// DLPK01
    odepack_cpp_real &delt = dlpk_.delt, &epcon = dlpk_.epcon, &sqrtn = dlpk_.sqrtn, &rsqrtn = dlpk_.rsqrtn;
    int &jpre = dlpk_.jpre, &jacflg = dlpk_.jacflg, &locwp = dlpk_.locwp, &lociwp = dlpk_.lociwp, &lsavx = dlpk_.lsavx, &kmp = dlpk_.kmp, &maxl = dlpk_.maxl, &mnewt = dlpk_.mnewt,
        &nni = dlpk_.nni, &nli = dlpk_.nli, &nps = dlpk_.nps, &ncfn = dlpk_.ncfn, &ncfl = dlpk_.ncfl;
//
    int iflag, lb, ldl, lhes, liom, lgmr, lpcg, lp, lq, lr, lv, lw, lwk, lz, maxlp1, npsl;
    odepack_cpp_real delta, hl0;
//
    iersl = 0;
    hl0 = h * el0;
    delta = delt * epcon;
    if (miter == 1) {
        goto LABEL_100;
    } else if (miter == 2) {
        goto LABEL_200;
    } else if (miter == 3) {
        goto LABEL_300;
    } else if (miter == 4) {
        goto LABEL_400;
    } else {
        goto LABEL_900;
    }
//-----------------------------------------------------------------------
// Use the SPIOM algorithm to solve the linear system P*x = -f.
//-----------------------------------------------------------------------
LABEL_100:
    lv = 1;
    lb = lv + n * maxl;
    lhes = lb + n;
    lwk = lhes + maxl * maxl;
    DCOPY(n, x, 1, &ARRAY1D(wm, lb), 1);
    DSCAL(n, rsqrtn, ewt, 1);
    DSPIOM (neq, tn, y, savf, &ARRAY1D(wm, lb), ewt, n, maxl, kmp, delta,
        hl0, jpre, mnewt, f, psol, npsl, x, &ARRAY1D(wm, lv), &ARRAY1D(wm, lhes), iwm,
        liom, &ARRAY1D(wm, locwp), &ARRAY1D(iwm, lociwp), &ARRAY1D(wm, lwk), iflag, user_data);
    nni = nni + 1;
    nli = nli + liom;
    nps = nps + npsl;
    DSCAL(n, sqrtn, ewt, 1);
    if (iflag != 0) ncfl = ncfl + 1;
    if (iflag >= 2) iersl = 1;
    if (iflag <  0) iersl = -1;
    return;
//-----------------------------------------------------------------------
// Use the SPIGMR algorithm to solve the linear system P*x = -f.
//-----------------------------------------------------------------------
LABEL_200:
    maxlp1 = maxl + 1;
    lv = 1;
    lb = lv + n * maxl;
    lhes = lb + n + 1;
    lq = lhes + maxl * maxlp1;
    lwk = lq + 2 * maxl;
    ldl = lwk + std::min(1, maxl - kmp) * n;
    DCOPY(n, x, 1, &ARRAY1D(wm, lb), 1);
    DSCAL(n, rsqrtn, ewt, 1);
    DSPIGMR(neq, tn, y, savf, &ARRAY1D(wm, lb), ewt, n, maxl, maxlp1, kmp,
        delta, hl0, jpre, mnewt, f, psol, npsl, x, &ARRAY1D(wm, lv), &ARRAY1D(wm, lhes),
        &ARRAY1D(wm, lq), lgmr, &ARRAY1D(wm, locwp), &ARRAY1D(iwm, lociwp), &ARRAY1D(wm, lwk), 
        &ARRAY1D(wm, ldl), iflag, user_data);
    nni = nni + 1;
    nli = nli + lgmr;
    nps = nps + npsl;
    DSCAL(n, sqrtn, ewt, 1);
    if (iflag != 0) ncfl  = ncfl + 1;
    if (iflag >= 2) iersl = 1;
    if (iflag <  0) iersl = -1;
    return;
//-----------------------------------------------------------------------
//Use DPCG to solve the linear system P*x = -f
//-----------------------------------------------------------------------
LABEL_300:
    lr = 1;
    lp = lr + n;
    lw = lp + n;
    lz = lw + n;
    lwk = lz + n;
    DCOPY(n, x, 1, &ARRAY1D(wm, lr), 1);
    DPCG(neq, tn, y, savf, &ARRAY1D(wm, lr), ewt, n, maxl, delta, hl0,
        jpre, mnewt, f, psol, npsl, x, &ARRAY1D(wm, lp), &ARRAY1D(wm, lw), &ARRAY1D(wm, lz),
        lpcg, &ARRAY1D(wm, locwp), &ARRAY1D(iwm, lociwp), &ARRAY1D(wm, lwk), iflag, user_data);
    nni = nni + 1;
    nli = nli + lpcg;
    nps = nps + npsl;
    if (iflag != 0) ncfl  = ncfl + 1;
    if (iflag >= 2) iersl = 1;
    if (iflag <  0) iersl = -1;
    return;
//-----------------------------------------------------------------------
// Use DPCGS to solve the linear system P*x = -f
//-----------------------------------------------------------------------
LABEL_400:
    lr = 1;
    lp = lr + n;
    lw = lp + n;
    lz = lw + n;
    lwk = lz + n;
    DCOPY(n, x, 1, &ARRAY1D(wm, lr), 1);
    DPCGS(neq, tn, y, savf, &ARRAY1D(wm, lr), ewt, n, maxl,delta, hl0, 
        jpre, mnewt, f, psol, npsl, x, &ARRAY1D(wm, lp), &ARRAY1D(wm, lw), &ARRAY1D(wm, lz),
        lpcg, &ARRAY1D(wm, locwp), &ARRAY1D(iwm, lociwp), &ARRAY1D(wm, lwk), iflag, user_data);
    nni = nni + 1;
    nli = nli + lpcg;
    nps = nps + npsl;
    if (iflag != 0) ncfl  = ncfl + 1;
    if (iflag >= 2) iersl = 1;
    if (iflag <  0) iersl = -1;
    return;
//-----------------------------------------------------------------------
// Use DUSOL, which interfaces to psol, to solve the linear system
// (no Krylov iteration).
//-----------------------------------------------------------------------
LABEL_900:
    lb = 1;
    lwk = lb + n;
    DCOPY(n, x, 1, &ARRAY1D(wm, lb), 1);
    DUSOL(neq, tn, y, savf, &ARRAY1D(wm, lb), ewt, n, delta, hl0, mnewt,
        psol, npsl, x, &ARRAY1D(wm, locwp), &ARRAY1D(iwm, lociwp), &ARRAY1D(wm, lwk), iflag, user_data);
    nni = nni + 1;
    nps = nps + npsl;
    if (iflag != 0) ncfl  = ncfl + 1;
    if (iflag == 3) iersl = 1;
    if (iflag <  0) iersl = -1;
    return;
}


/**
 * @fn DSPIOM
 * 
C-----------------------------------------------------------------------
C This routine solves the linear system A * x = b using a scaled
C preconditioned version of the Incomplete Orthogonalization Method.
C An initial guess of x = 0 is assumed.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C           TN = current value of t.
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C         B    = the right hand side of the system A*x = b.
C                B is also used as work space when computing the
C                final approximation.
C                (B is the same as V(*,MAXL+1) in the call to DSPIOM.)
C
C         WGHT = array of length N containing scale factors.
C                1/WGHT(i) are the diagonal elements of the diagonal
C                scaling matrix D.
C
C         N    = the order of the matrix A, and the lengths
C                of the vectors Y, SAVF, B, WGHT, and X.
C
C         MAXL = the maximum allowable order of the matrix HES.
C
C          KMP = the number of previous vectors the new vector VNEW
C                must be made orthogonal to.  KMP .le. MAXL.
C
C        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
C
C          HL0 = current value of (step size h) * (coefficient l0).
C
C         JPRE = preconditioner type flag.
C
C        MNEWT = Newton iteration counter (.ge. 0).
C
C           WK = real work array of length N used by DATV and PSOL.
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C      On return
C
C         X    = the final computed approximation to the solution
C                of the system A*x = b.
C
C         V    = the N by (LIOM+1) array containing the LIOM
C                orthogonal vectors V(*,1) to V(*,LIOM).
C
C         HES  = the LU factorization of the LIOM by LIOM upper
C                Hessenberg matrix whose entries are the
C                scaled inner products of A*V(*,k) and V(*,i).
C
C         IPVT = an integer array containg pivoting information.
C                It is loaded in DHEFA and used in DHESL.
C
C         LIOM = the number of iterations performed, and current
C                order of the upper Hessenberg matrix HES.
C
C         NPSL = the number of calls to PSOL.
C
C        IFLAG = integer error flag:
C                0 means convergence in LIOM iterations, LIOM.le.MAXL.
C                1 means the convergence test did not pass in MAXL
C                  iterations, but the residual norm is .lt. 1,
C                  or .lt. norm(b) if MNEWT = 0, and so X is computed.
C                2 means the convergence test did not pass in MAXL
C                  iterations, residual .gt. 1, and X is undefined.
C                3 means there was a recoverable error in PSOL
C                  caused by the preconditioner being out of date.
C               -1 means there was a nonrecoverable error in PSOL.
C
C-----------------------------------------------------------------------
 */
void Odepack::DSPIOM(int neq, odepack_cpp_real tn, odepack_cpp_real *y, odepack_cpp_real *savf, odepack_cpp_real *b, odepack_cpp_real *wght, int n, int maxl, int kmp,
    odepack_cpp_real &delta, odepack_cpp_real hl0, int jpre, int &mnewt, ODEPACK_FUNCTION f, ODEPACK_PSOL psol, int &npsl, odepack_cpp_real *x,
    odepack_cpp_real *v, odepack_cpp_real *hes, int *ipvt, int &liom, odepack_cpp_real *wp, int *iwp, odepack_cpp_real *wk, int &iflag, void *user_data)
{
#ifndef V
#define V(i, j) ARRAY2D(v, n, i, j)
#endif
#ifndef HES
#define HES(i, j) ARRAY2D(hes, maxl, i, j)
#endif
//
    int i, ier, info, j, k, ll, lm1;
    odepack_cpp_real bnrm, bnrm0, prod, rho, snormw, tem;
//
    iflag = 0;
    liom = 0;
    npsl = 0;
//-----------------------------------------------------------------------
// The initial residual is the vector b.  Apply scaling to b, and test
// for an immediate return with X = 0 or X = b.
//-----------------------------------------------------------------------
    for (i = 1; i <= n; ++i) {
        V(i, 1) = ARRAY1D(b, i) * ARRAY1D(wght, i);
    }
    bnrm0 = DNRM2(n, v, 1);
    bnrm = bnrm0;
    if (bnrm0 > delta) goto LABEL_30;
    if (mnewt > 0) goto LABEL_20;
    DCOPY(n, b, 1, x, 1);
    return;
LABEL_20:
    for (i = 1; i <= n; ++i) {
        ARRAY1D(x, i) = 0.0;
    }
    return;
LABEL_30:
// Apply inverse of left preconditioner to vector b. --------------------
    ier = 0;
    if (jpre == 0 || jpre == 2) goto LABEL_55;
    (*psol)(neq, tn, y, savf, wk, hl0, wp, iwp, b, 1, ier, user_data);
    npsl = 1;
    if (ier != 0) goto LABEL_300;
// Calculate norm of scaled vector V(*,1) and normalize it. -------------
    for (i = 1; i <= n; ++i) {
        V(i, 1) = ARRAY1D(b, i) * ARRAY1D(wght, i);
    }
    bnrm = DNRM2(n, v, 1);
    delta = delta * (bnrm / bnrm0);
LABEL_55:
    tem = 1.0 / bnrm;
    DSCAL(n, tem, &V(1, 1), 1);
// Zero out the HES array. ----------------------------------------------
    for (j = 1; j <= maxl; ++j) {
        for (i = 1; i <= maxl; ++i) {
            HES(i, j) = 0.0;
        }
    }
//-----------------------------------------------------------------------
// Main loop on LL = l to compute the vectors V(*,2) to V(*,MAXL).
// The running product PROD is needed for the convergence test.
//-----------------------------------------------------------------------
    prod = 1.0;
    for (ll = 1; ll <= maxl; ++ll) {
        liom = ll;
//-----------------------------------------------------------------------
// Call routine DATV to compute VNEW = Abar*v(l), where Abar is
// the matrix A with scaling and inverse preconditioner factors applied.
// Call routine DORTHOG to orthogonalize the new vector vnew = V(*,l+1).
// Call routine DHEFA to update the factors of HES.
//-----------------------------------------------------------------------
        DATV(neq, y, savf, &V(1, ll), wght, x, f, psol, &V(1, ll+1),
            wk, wp, iwp, hl0, jpre, ier, npsl, user_data);
        if (ier != 0) goto LABEL_300;
        DORTHOG(&V(1, ll+1), v, hes, n, ll, maxl, kmp, snormw);
        DHEFA(hes, maxl, ll, ipvt, info, ll);
        lm1 = ll - 1;
        if (ll > 1 && ARRAY1D(ipvt, lm1) == lm1) prod = prod * HES(ll, lm1);
        if (info != ll) goto LABEL_70;
//-----------------------------------------------------------------------
// The last pivot in HES was found to be zero.
// If vnew = 0 or l = MAXL, take an error return with IFLAG = 2.
// otherwise, continue the iteration without a convergence test.
//-----------------------------------------------------------------------
        if (snormw == 0.0) goto LABEL_120;
        if (ll == maxl) goto LABEL_120;
        goto LABEL_80;
//-----------------------------------------------------------------------
// Update RHO, the estimate of the norm of the residual b - A*x(l).
// test for convergence.  If passed, compute approximation x(l).
// If failed and l .lt. MAXL, then continue iterating.
//-----------------------------------------------------------------------
LABEL_70:
        rho = bnrm * snormw * std::abs(prod / HES(ll, ll));
        if (rho <= delta) goto LABEL_200;
        if (ll == maxl) goto LABEL_100;
// If l .lt. MAXL, store HES(l+1,l) and normalize the vector v(*,l+1).
LABEL_80:
        HES(ll+1, ll) = snormw;
        tem = 1.0 / snormw;
        DSCAL(n, tem, &V(1, ll+1), 1);
    }
//-----------------------------------------------------------------------
// l has reached MAXL without passing the convergence test:
// If RHO is not too large, compute a solution anyway and return with
// IFLAG = 1.  Otherwise return with IFLAG = 2.
//-----------------------------------------------------------------------
LABEL_100:
    if (rho <= 1.0) goto LABEL_150;
    if (rho <= bnrm && mnewt == 0) goto LABEL_150;
LABEL_120:
    iflag = 2;
    return;
LABEL_150:
    iflag = 1;
//-----------------------------------------------------------------------
// Compute the approximation x(l) to the solution.
// Since the vector X was used as work space, and the initial guess
// of the Newton correction is zero, X must be reset to zero.
//-----------------------------------------------------------------------
LABEL_200:
    ll = liom;
    for (k = 1; k <= ll; ++k) {
        ARRAY1D(b, k) = 0.0;
    }
    ARRAY1D(b, 1) = bnrm;
    DHESL(hes, maxl, ll, ipvt, b);
    for (k = 1; k <= n; ++k) {
        ARRAY1D(x, k) = 0.0;
    }
    for (i = 1; i <= ll; ++i) {
        DAXPY(n, ARRAY1D(b, i), &V(1, i), 1, x, 1);
    }
    for (i = 1; i <= n; ++i) {
        ARRAY1D(x, i) = ARRAY1D(x, i) / ARRAY1D(wght, i);
    }
    if (jpre <= 1) return;
    (*psol)(neq, tn, y, savf, wk, hl0, wp, iwp, x, 2, ier, user_data);
    npsl = npsl + 1;
    if (ier != 0) goto LABEL_300;
    return;
//-----------------------------------------------------------------------
// This block handles error returns forced by routine PSOL.
//-----------------------------------------------------------------------
LABEL_300:
    if (ier < 0) iflag = -1;
    if (ier > 0) iflag = 3;
//
#ifdef V
#undef V
#endif
#ifdef HES
#undef HES
#endif
}


/**
 * @fn DATV
C-----------------------------------------------------------------------
C This routine computes the product
C
C   (D-inverse)*(P1-inverse)*(I - hl0*df/dy)*(P2-inverse)*(D*v),
C
C where D is a diagonal scaling matrix, and P1 and P2 are the
C left and right preconditioning matrices, respectively.
C v is assumed to have WRMS norm equal to 1.
C The product is stored in z.  This is computed by a
C difference quotient, a call to F, and two calls to PSOL.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            V = real array of length N (can be the same array as Z).
C
C         WGHT = array of length N containing scale factors.
C                1/WGHT(i) are the diagonal elements of the matrix D.
C
C         FTEM = work array of length N.
C
C         VTEM = work array of length N used to store the
C                unscaled version of V.
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C          HL0 = current value of (step size h) * (coefficient l0).
C
C         JPRE = preconditioner type flag.
C
C
C      On return
C
C            Z = array of length N containing desired scaled
C                matrix-vector product.
C
C          IER = error flag from PSOL.
C
C         NPSL = the number of calls to PSOL.
C
C In addition, this routine uses the Common variables TN, N, NFE.
C-----------------------------------------------------------------------
 */
void Odepack::DATV(int neq, odepack_cpp_real *y, odepack_cpp_real *savf, odepack_cpp_real *v, odepack_cpp_real *wght, odepack_cpp_real *ftem, ODEPACK_FUNCTION f,
    ODEPACK_PSOL psol, odepack_cpp_real *z, odepack_cpp_real *vtem, odepack_cpp_real *wp, int *iwp, odepack_cpp_real hl0, int &jpre, int &ier, int &npsl,
    void *user_data)
{
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
//
    int i;
    odepack_cpp_real fac, rnorm, tempn;
//
// Set VTEM = D * V.
    for (i = 1; i <= n; ++i) {
        ARRAY1D(vtem, i) = ARRAY1D(v, i) / ARRAY1D(wght, i);
    }
    ier = 0;
    if (jpre >= 2) goto LABEL_30;
//
// JPRE = 0 or 1.  Save Y in Z and increment Y by VTEM.
    DCOPY(n, y, 1, z, 1);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = ARRAY1D(z, i) + ARRAY1D(vtem, i);
    }
    fac = hl0;
    goto LABEL_60;
//
// JPRE = 2 or 3.  Apply inverse of right preconditioner to VTEM.
LABEL_30:
    (*psol)(neq, tn, y, savf, ftem, hl0, wp, iwp, vtem, 2, ier, user_data);
    npsl = npsl + 1;
    if (ier != 0) return;
// Calculate L-2 norm of (D-inverse) * VTEM.
    for (i = 1; i <= n; ++i) {
        ARRAY1D(z, i) = ARRAY1D(vtem, i) * ARRAY1D(wght,i);
    }
    tempn = DNRM2(n, z, 1);
    rnorm = 1.0 / tempn;
// Save Y in Z and increment Y by VTEM/norm.
    DCOPY(n, y, 1, z, 1);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = ARRAY1D(z, i) + ARRAY1D(vtem, i) * rnorm;
    }
    fac = hl0 * tempn;
//
// For all JPRE, call F with incremented Y argument, and restore Y.
LABEL_60:
    (*f)(neq, tn, y, ftem, user_data);
    nfe = nfe + 1;
    DCOPY(n, z, 1, y, 1);
// Set Z = (identity - hl0*Jacobian) * VTEM, using difference quotient.
    for (i = 1; i <= n; ++i) {
        ARRAY1D(z, i) = ARRAY1D(ftem, i) - ARRAY1D(savf, i);
    }
    for (i = 1; i <= n; ++i) {
        ARRAY1D(z, i) = ARRAY1D(vtem, i) - fac * ARRAY1D(z, i);
    }
// Apply inverse of left preconditioner to Z, if nontrivial.
    if (jpre == 0 || jpre == 2) goto LABEL_85;
    (*psol)(neq, n, y, savf, ftem, hl0, wp, iwp, z, 1, ier, user_data);
    npsl = npsl + 1;
    if (ier != 0) return;
LABEL_85:
// Apply D-inverse to Z and return.
    for (i = 1; i <= n; ++i) {
        ARRAY1D(z, i) = ARRAY1D(z, i) * ARRAY1D(wght, i);
    }
    return;
}


/**
 * @fn DORTHOG
 * 
C-----------------------------------------------------------------------
C This routine orthogonalizes the vector VNEW against the previous
C KMP vectors in the V array.  It uses a modified Gram-Schmidt
C orthogonalization procedure with conditional reorthogonalization.
C This is the version of 28 may 1986.
C-----------------------------------------------------------------------
C
C      On entry
C
C         VNEW = the vector of length N containing a scaled product
C                of the Jacobian and the vector V(*,LL).
C
C         V    = the N x l array containing the previous LL
C                orthogonal vectors v(*,1) to v(*,LL).
C
C         HES  = an LL x LL upper Hessenberg matrix containing,
C                in HES(i,k), k.lt.LL, scaled inner products of
C                A*V(*,k) and V(*,i).
C
C        LDHES = the leading dimension of the HES array.
C
C         N    = the order of the matrix A, and the length of VNEW.
C
C         LL   = the current order of the matrix HES.
C
C          KMP = the number of previous vectors the new vector VNEW
C                must be made orthogonal to (KMP .le. MAXL).
C
C
C      On return
C
C         VNEW = the new vector orthogonal to V(*,i0) to V(*,LL),
C                where i0 = MAX(1, LL-KMP+1).
C
C         HES  = upper Hessenberg matrix with column LL filled in with
C                scaled inner products of A*V(*,LL) and V(*,i).
C
C       SNORMW = L-2 norm of VNEW.
C
C-----------------------------------------------------------------------
 */
void Odepack::DORTHOG(odepack_cpp_real *vnew, odepack_cpp_real *v, odepack_cpp_real *hes, int n, int ll, int ldhes, int kmp, odepack_cpp_real &snormw)
{
#ifndef V
#define V(i, j) ARRAY2D(v, n, i, j)
#endif
#ifndef HES
#define HES(i, j) ARRAY2D(hes, ldhes, i, j)
#endif
//
    int i, i0;
    odepack_cpp_real arg, sumdsq, tem, vnrm;
//
    vnrm = DNRM2(n, vnew, 1);
//-----------------------------------------------------------------------
// Do modified Gram-Schmidt on VNEW = A*v(LL).
// Scaled inner products give new column of HES.
// Projections of earlier vectors are subtracted from VNEW.
//-----------------------------------------------------------------------
    i0 = std::max(1, ll - kmp + 1);
    for (i = i0; i <= ll; ++i) {
        HES(i, ll) = DDOT(n, &V(1, i), 1, vnew, 1);
        tem = - HES(i, ll);
        DAXPY(n, tem, &V(1, i), 1, vnew, 1);
    }
//-----------------------------------------------------------------------
// Compute SNORMW = norm of VNEW.
// If VNEW is small compared to its input value (in norm), then
// reorthogonalize VNEW to V(*,1) through V(*,LL).
// Correct if relative correction exceeds 1000*(unit roundoff).
// finally, correct SNORMW using the dot products involved.
//-----------------------------------------------------------------------
    snormw = DNRM2(n, vnew, 1);
    if ((vnrm + 0.001 * snormw) != vnrm) return;
    sumdsq = 0.0;
    for (i = i0; i <= ll; ++i) {
        tem = - DDOT(n, &V(1, i), 1, vnew, 1);
        if ((HES(i, ll) + 0.001 * tem) == HES(i, ll)) continue;
        HES(i, ll) = HES(i, ll) - tem;
        DAXPY(n, tem, &V(1, i), 1, vnew, 1);
        sumdsq = sumdsq + tem * tem;
    }
    if (sumdsq == 0.0) return;
    arg = std::max(0.0, snormw * snormw - sumdsq);
    snormw = std::sqrt(arg);
//
    return;
//
#ifdef V
#undef V
#endif
#ifdef HES
#undef HES
#endif
}


/**
 * @fn DSPIGMR
 * 
C-----------------------------------------------------------------------
C This routine solves the linear system A * x = b using a scaled
C preconditioned version of the Generalized Minimal Residual method.
C An initial guess of x = 0 is assumed.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C           TN = current value of t.
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            B = the right hand side of the system A*x = b.
C                B is also used as work space when computing
C                the final approximation.
C                (B is the same as V(*,MAXL+1) in the call to DSPIGMR.)
C
C         WGHT = the vector of length N containing the nonzero
C                elements of the diagonal scaling matrix.
C
C            N = the order of the matrix A, and the lengths
C                of the vectors WGHT, B and X.
C
C         MAXL = the maximum allowable order of the matrix HES.
C
C       MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES.
C
C          KMP = the number of previous vectors the new vector VNEW
C                must be made orthogonal to.  KMP .le. MAXL.
C
C        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
C
C          HL0 = current value of (step size h) * (coefficient l0).
C
C         JPRE = preconditioner type flag.
C
C        MNEWT = Newton iteration counter (.ge. 0).
C
C           WK = real work array used by routine DATV and PSOL.
C
C           DL = real work array used for calculation of the residual
C                norm RHO when the method is incomplete (KMP .lt. MAXL).
C                Not needed or referenced in complete case (KMP = MAXL).
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C      On return
C
C         X    = the final computed approximation to the solution
C                of the system A*x = b.
C
C         LGMR = the number of iterations performed and
C                the current order of the upper Hessenberg
C                matrix HES.
C
C         NPSL = the number of calls to PSOL.
C
C         V    = the N by (LGMR+1) array containing the LGMR
C                orthogonal vectors V(*,1) to V(*,LGMR).
C
C         HES  = the upper triangular factor of the QR decomposition
C                of the (LGMR+1) by lgmr upper Hessenberg matrix whose
C                entries are the scaled inner-products of A*V(*,i)
C                and V(*,k).
C
C         Q    = real array of length 2*MAXL containing the components
C                of the Givens rotations used in the QR decomposition
C                of HES.  It is loaded in DHEQR and used in DHELS.
C
C        IFLAG = integer error flag:
C                0 means convergence in LGMR iterations, LGMR .le. MAXL.
C                1 means the convergence test did not pass in MAXL
C                  iterations, but the residual norm is .lt. 1,
C                  or .lt. norm(b) if MNEWT = 0, and so x is computed.
C                2 means the convergence test did not pass in MAXL
C                  iterations, residual .gt. 1, and X is undefined.
C                3 means there was a recoverable error in PSOL
C                  caused by the preconditioner being out of date.
C               -1 means there was a nonrecoverable error in PSOL.
C
C-----------------------------------------------------------------------
 */
void Odepack::DSPIGMR(
    int neq, odepack_cpp_real tn, odepack_cpp_real *y, odepack_cpp_real *savf, odepack_cpp_real *b, odepack_cpp_real *wght, int n, int maxl, int maxlp1,
    int kmp, odepack_cpp_real &delta, odepack_cpp_real hl0, int jpre, int mnewt, ODEPACK_FUNCTION f, ODEPACK_PSOL psol, int &npsl, odepack_cpp_real *x, odepack_cpp_real *v, odepack_cpp_real *hes, odepack_cpp_real *q,
    int &lgmr, odepack_cpp_real *wp, int *iwp, odepack_cpp_real *wk, odepack_cpp_real *dl, int &iflag, void *user_data)
{
#ifndef V
#define V(i, j) ARRAY2D(v, n, i, j)
#endif
#ifndef HES
#define HES(i, j) ARRAY2D(hes, maxlp1, i, j)
#endif
//
    int i, ier, info, ip1, i2, j, k, ll, llp1;
    odepack_cpp_real bnrm, bnrm0, c, dlnrm, prod, rho, s, snormw, tem;
//
    iflag = 0;
    lgmr = 0;
    npsl = 0;
//-----------------------------------------------------------------------
// The initial residual is the vector b.  Apply scaling to b, and test
// for an immediate return with X = 0 or X = b.
//-----------------------------------------------------------------------
    for (i = 1; i <= n; ++i) {
        V(i, 1) = ARRAY1D(b, i) * ARRAY1D(wght, i);
    }
    bnrm0 = DNRM2(n, v, 1);
    bnrm  = bnrm0;
    if (bnrm0 > delta) goto LABEL_30;
    if (mnewt > 0) goto LABEL_20;
    DCOPY(n, b, 1, x, 1);
    return;
LABEL_20:
    for (i = 1; i <= n; ++i) {
        ARRAY1D(x, i) = 0.0;
    }
    return;
LABEL_30:
// Apply inverse of left preconditioner to vector b. --------------------
    ier = 0;
    if (jpre == 0 || jpre == 2) goto LABEL_55;
    (*psol)(neq, tn, y, savf, wk, hl0, wp, iwp, b, 1, ier, user_data);
    npsl = 1;
    if (ier != 0) goto LABEL_300;
// Calculate norm of scaled vector V(*,1) and normalize it. -------------
    for (i = 1; i <= n; ++i) {
        V(i, 1) = ARRAY1D(b, i) * ARRAY1D(wght, i);
    }
    bnrm = DNRM2(n, v, 1);
    delta = delta * (bnrm / bnrm0);
LABEL_55:
    tem = 1.0 / bnrm;
    DSCAL(n, tem, &V(1, 1), 1);
// Zero out the HES array. ----------------------------------------------
    for (j = 1; j <= maxl; ++j) {
        for (i = 1; i <= maxlp1; ++i) {
            HES(i, j) = 0.0;
        }
    }
//-----------------------------------------------------------------------
// Main loop to compute the vectors V(*,2) to V(*,MAXL).
// The running product PROD is needed for the convergence test.
//-----------------------------------------------------------------------
    prod = 1.0;
    for (ll = 1; ll <= maxl; ++ll) {
        lgmr = ll;
//-----------------------------------------------------------------------
// Call routine DATV to compute VNEW = Abar*v(ll), where Abar is
// the matrix A with scaling and inverse preconditioner factors applied.
// Call routine DORTHOG to orthogonalize the new vector VNEW = V(*,LL+1).
// Call routine DHEQR to update the factors of HES.
//-----------------------------------------------------------------------
        DATV(neq, y, savf, &V(1, ll), wght, x, f, psol, &V(1, ll+1),
            wk, wp, iwp, hl0, jpre, ier, npsl, user_data);
        if (ier != 0) goto LABEL_300;
        DORTHOG(&V(1, ll+1), v, hes, n, ll, maxlp1, kmp, snormw);
        HES(ll+1, ll) = snormw;
        DHEQR(hes, maxlp1, ll, q, info, ll);
        if (info == ll) goto LABEL_120;
//-----------------------------------------------------------------------
// Update RHO, the estimate of the norm of the residual b - A*xl.
// If KMP .lt. MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
// necessarily orthogonal for LL .gt. KMP.  The vector DL must then
// be computed, and its norm used in the calculation of RHO.
//-----------------------------------------------------------------------
        prod = prod * ARRAY1D(q, 2 * ll);
        rho = std::abs(prod * bnrm);
        if (ll > kmp && kmp < maxl) {
            if (ll == (kmp + 1)) {
                DCOPY(n, &V(1, 1), 1, dl, 1);
                for (i = 1; i <= kmp; ++i) {
                    ip1 = i + 1;
                    i2  = i * 2;
                    s   = ARRAY1D(q, i2);
                    c   = ARRAY1D(q, i2 - 1);
                    for (k = 1; k <= n; ++k) {
                        ARRAY1D(dl, k) = s * ARRAY1D(dl, k) + c * V(k, ip1);
                    }
                }
            }
            s = ARRAY1D(q, 2 * ll);
            c = ARRAY1D(q, 2 * ll - 1) / snormw;
            llp1 = ll + 1;
            for (k = 1; k <= n; ++k) {
                ARRAY1D(dl, k) = s * ARRAY1D(dl, k) + c * V(k, llp1);
            }
            dlnrm = DNRM2(n, dl, 1);
            rho = rho * dlnrm;
        }
//-----------------------------------------------------------------------
// Test for convergence.  If passed, compute approximation xl.
// if failed and LL .lt. MAXL, then continue iterating.
//-----------------------------------------------------------------------
        if (rho <= delta) goto LABEL_200;
        if (ll == maxl) goto LABEL_100;
//-----------------------------------------------------------------------
// Rescale so that the norm of V(1,LL+1) is one.
//-----------------------------------------------------------------------
        tem = 1.0 / snormw;
        DSCAL(n, tem, &V(1, ll + 1), 1);
    }
LABEL_100:
    if (rho <= 1.0) goto LABEL_150;
    if (rho <= bnrm && mnewt == 0) goto LABEL_150;
LABEL_120:
    iflag = 2;
    return;
LABEL_150:
    return;
//-----------------------------------------------------------------------
// Compute the approximation xl to the solution.
// Since the vector X was used as work space, and the initial guess
// of the Newton correction is zero, X must be reset to zero.
//-----------------------------------------------------------------------
LABEL_200:
    ll = lgmr;
    llp1 = ll + 1;
    for (k = 1; k <= llp1; ++k) {
        ARRAY1D(b, k) = 0.0;
    }
    ARRAY1D(b, 1) = bnrm;
    DHELS(hes, maxlp1, ll, q, b);
    for (k = 1; k <= n; ++k) {
        ARRAY1D(x, k) = 0.0;
    }
    for (i = 1; i <= ll; ++i) {
        DAXPY(n, ARRAY1D(b, i), &V(1, i), 1, x, 1);
    }
    for (i = 1; i <= n; ++i) {
        ARRAY1D(x, i) = ARRAY1D(x, i) / ARRAY1D(wght, i);
    }
    if (jpre <= 1) return;
    (*psol)(neq, tn, y, savf, wk, hl0, wp, iwp, x, 2, ier, user_data);
    npsl = npsl + 1;
    if (ier != 0) goto LABEL_300;
    return;
//-----------------------------------------------------------------------
// This block handles error returns forced by routine PSOL.
//-----------------------------------------------------------------------
LABEL_300:
    if (ier < 0) iflag = -1;
    if (ier > 0) iflag = 3;
//
    return;
//
#ifdef V
#undef V
#endif
#ifdef HES
#undef HES
#endif
}


/**
 * @fn DPCG
 * 
C-----------------------------------------------------------------------
C This routine computes the solution to the system A*x = b using a
C preconditioned version of the Conjugate Gradient algorithm.
C It is assumed here that the matrix A and the preconditioner
C matrix M are symmetric positive definite or nearly so.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C           TN = current value of t.
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            R = the right hand side of the system A*x = b.
C
C         WGHT = array of length N containing scale factors.
C                1/WGHT(i) are the diagonal elements of the diagonal
C                scaling matrix D.
C
C            N = the order of the matrix A, and the lengths
C                of the vectors Y, SAVF, R, WGHT, P, W, Z, WK, and X.
C
C         MAXL = the maximum allowable number of iterates.
C
C        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
C
C          HL0 = current value of (step size h) * (coefficient l0).
C
C         JPRE = preconditioner type flag.
C
C        MNEWT = Newton iteration counter (.ge. 0).
C
C           WK = real work array used by routine DATP.
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C      On return
C
C         X    = the final computed approximation to the solution
C                of the system A*x = b.
C
C         LPCG = the number of iterations performed, and current
C                order of the upper Hessenberg matrix HES.
C
C         NPSL = the number of calls to PSOL.
C
C        IFLAG = integer error flag:
C                0 means convergence in LPCG iterations, LPCG .le. MAXL.
C                1 means the convergence test did not pass in MAXL
C                  iterations, but the residual norm is .lt. 1,
C                  or .lt. norm(b) if MNEWT = 0, and so X is computed.
C                2 means the convergence test did not pass in MAXL
C                  iterations, residual .gt. 1, and X is undefined.
C                3 means there was a recoverable error in PSOL
C                  caused by the preconditioner being out of date.
C                4 means there was a zero denominator in the algorithm.
C                  The system matrix or preconditioner matrix is not
C                  sufficiently close to being symmetric pos. definite.
C               -1 means there was a nonrecoverable error in PSOL.
C
C-----------------------------------------------------------------------
 */
void Odepack::DPCG(int neq, odepack_cpp_real tn, odepack_cpp_real *y, odepack_cpp_real *savf, odepack_cpp_real *r, odepack_cpp_real *wght, int n, int maxl,
    odepack_cpp_real delta, odepack_cpp_real hl0, int &jpre, int mnewt, ODEPACK_FUNCTION f, ODEPACK_PSOL psol,
    int &npsl, odepack_cpp_real *x, odepack_cpp_real *p, odepack_cpp_real *w, odepack_cpp_real *z, int &lpcg, odepack_cpp_real *wp,int *iwp, odepack_cpp_real *wk, int &iflag,
    void *user_data)
{
    int i, ier;
    odepack_cpp_real alpha, beta, bnrm, ptw, rnrm, ztr, ztr0;
//
    iflag = 0;
    npsl = 0;
    lpcg = 0;
    for (i = 1; i <= n; ++i) {
        x[i] = 0.0;
    }
    bnrm = DVNORM(n, r, wght);
// Test for immediate return with X = 0 or X = b. -----------------------
    if (bnrm > delta) goto LABEL_20;
    if (mnewt > 0) return;
    DCOPY(n, r, 1, x, 1);
    return;
//
LABEL_20:
    ztr = 0.0;
// Loop point for PCG iterations. ---------------------------------------
LABEL_30:
    lpcg++;
    DCOPY(n, r, 1, z, 1);
    ier = 0;
    if (jpre == 0) goto LABEL_40;
    (*psol)(neq, tn, y, savf, wk, hl0, wp, iwp, z, 3, ier, user_data);
    npsl++;
    if (ier != 0) goto LABEL_100;
LABEL_40:
    ztr0 = ztr;
    ztr = DDOT(n, z, 1, r, 1);
    if (lpcg != 1) goto LABEL_50;
    DCOPY(n, z, 1, p, 1);
    goto LABEL_70;
LABEL_50:
    if (ztr0 == 0.0) goto LABEL_200;
    beta = ztr / ztr0;
    for (i = 1; i <= n; ++i) {
        p[i] = z[i] + beta * p[i];
    }
LABEL_70:
//-----------------------------------------------------------------------
//  Call DATP to compute A*p and return the answer in W.
//-----------------------------------------------------------------------
    DATP(neq, y, savf, p, wght, hl0, wk, f, w, user_data);
//
    ptw = DDOT(n, p, 1, w, 1);
    if (ptw == 0.0) goto LABEL_200;
    alpha = ztr / ztr0;
    DAXPY(n, alpha, p, 1, x, 1);
    alpha *= -1.0;
    DAXPY(n, alpha, w, 1, r, 1);
    rnrm = DVNORM(n, r, wght);
    if (rnrm <= delta) return;
    if (lpcg < maxl) goto LABEL_30;
    iflag = 2;
    if (rnrm <= 1.0) iflag = 1;
    if (rnrm <= bnrm && mnewt == 0) iflag = 1;
    return;
//-----------------------------------------------------------------------
// This block handles error returns from PSOL.
//-----------------------------------------------------------------------
LABEL_100:
    if (ier < 0) iflag = -1;
    if (ier > 0) iflag = 3;
    return;
//-----------------------------------------------------------------------
// This block handles division by zero errors.
//-----------------------------------------------------------------------
LABEL_200:
    iflag = 4;
    return;
}


/**
 * @fn DPCGS
 * 
C-----------------------------------------------------------------------
C This routine computes the solution to the system A*x = b using a
C scaled preconditioned version of the Conjugate Gradient algorithm.
C It is assumed here that the scaled matrix D**-1 * A * D and the
C scaled preconditioner D**-1 * M * D are close to being
C symmetric positive definite.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C           TN = current value of t.
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            R = the right hand side of the system A*x = b.
C
C         WGHT = array of length N containing scale factors.
C                1/WGHT(i) are the diagonal elements of the diagonal
C                scaling matrix D.
C
C            N = the order of the matrix A, and the lengths
C                of the vectors Y, SAVF, R, WGHT, P, W, Z, WK, and X.
C
C         MAXL = the maximum allowable number of iterates.
C
C        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
C
C          HL0 = current value of (step size h) * (coefficient l0).
C
C         JPRE = preconditioner type flag.
C
C        MNEWT = Newton iteration counter (.ge. 0).
C
C           WK = real work array used by routine DATP.
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C      On return
C
C         X    = the final computed approximation to the solution
C                of the system A*x = b.
C
C         LPCG = the number of iterations performed, and current
C                order of the upper Hessenberg matrix HES.
C
C         NPSL = the number of calls to PSOL.
C
C        IFLAG = integer error flag:
C                0 means convergence in LPCG iterations, LPCG .le. MAXL.
C                1 means the convergence test did not pass in MAXL
C                  iterations, but the residual norm is .lt. 1,
C                  or .lt. norm(b) if MNEWT = 0, and so X is computed.
C                2 means the convergence test did not pass in MAXL
C                  iterations, residual .gt. 1, and X is undefined.
C                3 means there was a recoverable error in PSOL
C                  caused by the preconditioner being out of date.
C                4 means there was a zero denominator in the algorithm.
C                  the scaled matrix or scaled preconditioner is not
C                  sufficiently close to being symmetric pos. definite.
C               -1 means there was a nonrecoverable error in PSOL.
C
C-----------------------------------------------------------------------
 */
void Odepack::DPCGS(int neq, odepack_cpp_real tn, odepack_cpp_real *y, odepack_cpp_real *savf, odepack_cpp_real *r, odepack_cpp_real *wght, int n, int maxl,
    odepack_cpp_real delta, odepack_cpp_real hl0, int jpre, int mnewt, ODEPACK_FUNCTION f, ODEPACK_PSOL psol,
    int &npsl, odepack_cpp_real *x, odepack_cpp_real *p, odepack_cpp_real *w, odepack_cpp_real *z, int &lpcg, odepack_cpp_real *wp, int *iwp, odepack_cpp_real *wk, int &iflag,
    void *user_data)
{
    int i, ier;
    odepack_cpp_real alpha, beta, bnrm, ptw, rnrm, ztr, ztr0;
//
    iflag = 0;
    npsl = 0;
    lpcg = 0;
    for (i = 1; i <=  n; ++i) {
        x[i] = 0.0;
    }
    bnrm = DVNORM(n, r, wght);
// Test for immediate return with X = 0 or X = b. -----------------------
    if (bnrm > delta) goto LABEL_20;
    if (mnewt > 0) return;
    DCOPY(n, r, 1, x, 1);
    return;
//
LABEL_20:
    ztr = 0.0;
// Loop point for PCG iterations. ---------------------------------------
LABEL_30:
    lpcg++;
    DCOPY(n, r, 1, z, 1);
    ier = 0;
    if (jpre == 0) goto LABEL_40;
    (*psol)(neq, tn, y, savf, wk, hl0, wp, iwp, z, 3, ier, user_data);
    npsl++;
    if (ier != 0) goto LABEL_100;
LABEL_40:
    ztr0 = ztr;
    ztr = 0.0;
    for (i = 1; i <= n; ++i) {
        ztr += z[i] * r[i] * wght[i] * wght[i];
    }
    if (lpcg != 1) goto LABEL_50;
    DCOPY(n, z, 1, p, 1);
    goto LABEL_70;
LABEL_50:
    if (ztr0 == 0.0) goto LABEL_200;
    beta = ztr / ztr0;
    for (i = 1; i <= n; ++i) {
        p[i] = z[i] + beta * p[i];
    }
LABEL_70:
//-----------------------------------------------------------------------
//  Call DATP to compute A*p and return the answer in W.
//-----------------------------------------------------------------------
    DATP(neq, y, savf, p, wght, hl0, wk, f, w, user_data);
//
    ptw = 0.0;
    for (i = 1; i <= n; ++i) {
        ptw += p[i] * w[i] * wght[i] * wght[i];
    }
    if (ptw == 0.0) goto LABEL_200;
    alpha = ztr / ptw;
    DAXPY(n, alpha, p, 1, x, 1);
    alpha *= -1.0;
    DAXPY(n, alpha, w, 1, r, 1);
    rnrm = DVNORM(n, r, wght);
    if (rnrm <= delta) return;
    if (lpcg < maxl) goto LABEL_30;
    iflag = 2;
    if (rnrm <= 1.0) iflag = 1;
    if (rnrm <= bnrm && mnewt == 0) iflag = 1;
    return;
//-----------------------------------------------------------------------
// This block handles error returns from PSOL.
//-----------------------------------------------------------------------
LABEL_100:
    if (ier < 0) iflag = -1;
    if (ier > 0) iflag = 3;
    return;
//-----------------------------------------------------------------------
// This block handles division by zero errors.
//-----------------------------------------------------------------------
LABEL_200:
    iflag = 4;
    return;
}


/**
 * @fn DATP
 * 
C-----------------------------------------------------------------------
C This routine computes the product
C
C              w = (I - hl0*df/dy)*p
C
C This is computed by a call to F and a difference quotient.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            P = real array of length N.
C
C         WGHT = array of length N containing scale factors.
C                1/WGHT(i) are the diagonal elements of the matrix D.
C
C           WK = work array of length N.
C
C      On return
C
C
C            W = array of length N containing desired
C                matrix-vector product.
C
C In addition, this routine uses the Common variables TN, N, NFE.
C-----------------------------------------------------------------------
 */
void Odepack::DATP(int neq, odepack_cpp_real *y, odepack_cpp_real *savf, odepack_cpp_real *p, odepack_cpp_real *wght, odepack_cpp_real hl0, odepack_cpp_real *wk,
    ODEPACK_FUNCTION f, odepack_cpp_real *w, void *user_data)
{
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
//
    int i;
    odepack_cpp_real fac, pnrm, rpnrm;
//
    pnrm = DVNORM(n, p, wght);
    rpnrm = 1.0 / pnrm;
    DCOPY(n, y, 1, w, 1);
    for (i = 1; i <= n; ++i) {
        y[i] = w[i] + p[i] * rpnrm;
    }
    (*f)(neq, tn, y, wk, user_data);
    nfe++;
    DCOPY(n, w, 1, y, 1);
    fac = hl0 * pnrm;
    for (i = 1; i <= n; ++i) {
        w[i] = p[i] - fac * (wk[i] - savf[i]);
    }
    return;
}


/**
 * @fn DUSOL
 * 
C-----------------------------------------------------------------------
C This routine solves the linear system A * x = b using only a call
C to the user-supplied routine PSOL (no Krylov iteration).
C If the norm of the right-hand side vector b is smaller than DELTA,
C the vector X returned is X = b (if MNEWT = 0) or X = 0 otherwise.
C PSOL is called with an LR argument of 0.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C           TN = current value of t.
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            B = the right hand side of the system A*x = b.
C
C         WGHT = the vector of length N containing the nonzero
C                elements of the diagonal scaling matrix.
C
C            N = the order of the matrix A, and the lengths
C                of the vectors WGHT, B and X.
C
C        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
C
C          HL0 = current value of (step size h) * (coefficient l0).
C
C        MNEWT = Newton iteration counter (.ge. 0).
C
C           WK = real work array used by PSOL.
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C      On return
C
C         X    = the final computed approximation to the solution
C                of the system A*x = b.
C
C         NPSL = the number of calls to PSOL.
C
C        IFLAG = integer error flag:
C                0 means no trouble occurred.
C                3 means there was a recoverable error in PSOL
C                  caused by the preconditioner being out of date.
C               -1 means there was a nonrecoverable error in PSOL.
C
C-----------------------------------------------------------------------
 */
void Odepack::DUSOL(int neq, odepack_cpp_real tn, odepack_cpp_real *y, odepack_cpp_real *savf, odepack_cpp_real *b, odepack_cpp_real *wght, int n, odepack_cpp_real delta,
    odepack_cpp_real hl0, int mnewt, ODEPACK_PSOL psol, int &npsl, odepack_cpp_real *x, odepack_cpp_real *wp, int *iwp, odepack_cpp_real *wk, int &iflag,
    void *user_data)
{
    int i, ier;
    odepack_cpp_real bnrm;
//
    iflag = 0;
    npsl = 0;
//-----------------------------------------------------------------------
// Test for an immediate return with X = 0 or X = b.
//-----------------------------------------------------------------------
    bnrm = DVNORM(n, b, wght);
    if (bnrm > delta) goto LABEL_30;
    if (mnewt > 0) goto LABEL_10;
    DCOPY(n, b, 1, x, 1);
    return;
LABEL_10:
    for (i = 1; i <= n; ++i) {
        x[i] = 0.0;
    }
    return;
// Make call to PSOL and copy result from B to X. -----------------------
LABEL_30:
    ier = 0;
    (*psol)(neq, tn, y, savf, wk, hl0, wp, iwp, b, 0, ier, user_data);
    npsl = 1;
    if (ier != 0) goto LABEL_100;
    DCOPY(n, b, 1, x, 1);
    return;
//-----------------------------------------------------------------------
// This block handles error returns forced by routine PSOL.
//-----------------------------------------------------------------------
LABEL_100:
    if (ier < 0) iflag = -1;
    if (ier > 0) iflag = 3;
    return;
}


/**
 * @fn DSRCPK
 * 
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLPK01, which are used
C internally by the DLSODPK solver.
C
C RSAV = real array of length 222 or more.
C ISAV = integer array of length 50 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
 */
void Odepack::DSRCPK(odepack_cpp_real *rsav, int *isav, int job)
{
    int i;
// DLS001
    odepack_cpp_real *rls = dls1_.rls;
    int *ils = dls1_.ils;
// DLPK01
    odepack_cpp_real *rlsp = dlpk_.rlsp;
    int *ilsp = dlpk_.ilsp;
//
    int lenrls = 218;
    int lenils = 37;
    int lenrlp = 4;
    int lenilp = 13;
//
    if (job == 2) goto LABEL_100;
    DCOPY(lenrls, rls, 1, rsav, 1);
    DCOPY(lenrlp, rlsp, 1, &ARRAY1D(rsav, lenrls+1), 1);
    for (i = 1; i <= lenils; ++i) {
        ARRAY1D(isav, i) = ARRAY1D(ils, i);
    }
    for (i = 1; i <= lenilp; ++i) {
        ARRAY1D(isav, lenils + i) = ARRAY1D(ilsp, i);
    }
    return;
//
LABEL_100:
    DCOPY(lenrls, rsav, 1, rls, 1);
    DCOPY(lenrlp, &ARRAY1D(rsav, lenrls + 1), 1, rlsp, 1);
    for (i = 1; i <= lenils; ++i) {
        ARRAY1D(ils, i) = ARRAY1D(isav, i);
    }
    for (i = 1; i <= lenilp; ++i) {
        ARRAY1D(ilsp, i) = ARRAY1D(isav, lenils + i);
    }
    return;
}


/**
 * @fn DHEFA
 * 
C-----------------------------------------------------------------------
C     This routine is a modification of the LINPACK routine DGEFA and
C     performs an LU decomposition of an upper Hessenberg matrix A.
C     There are two options available:
C
C          (1)  performing a fresh factorization
C          (2)  updating the LU factors by adding a row and a
C               column to the matrix A.
C-----------------------------------------------------------------------
C     DHEFA factors an upper Hessenberg matrix by elimination.
C
C     On entry
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
C        JOB     INTEGER
C                JOB = 1    means that a fresh factorization of the
C                           matrix A is desired.
C                JOB .ge. 2 means that the current factorization of A
C                           will be updated by the addition of a row
C                           and a column.
C
C     On return
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
C                = k  if  U(k,k) .eq. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DHESL will divide by zero if called.
C
C     Modification of LINPACK, by Peter Brown, LLNL.
C     Written 7/20/83.  This version dated 6/20/01.
C    
C     BLAS called: DAXPY, IDAMAX
C-----------------------------------------------------------------------
 */
void Odepack::DHEFA(odepack_cpp_real *a, int lda, int n, int *ipvt, int &info, int job)
{
#ifndef MATA
#define MATA(i, j) ARRAY2D(a, lda, i, j)
#endif
//
    int j, k, km1, kp1, l, nm1;
    odepack_cpp_real t;
//
    if (job > 1) goto LABEL_80;
//
// A new facorization is desired.  This is essentially the LINPACK
// code with the exception that we know there is only one nonzero
// element below the main diagonal.
//
//     Gaussian elimination with partial pivoting
//
    info = 0;
    nm1 = n - 1;
    if (nm1 < 1) goto LABEL_70;
    for (k = 1; k <= nm1; ++k) {
        kp1 = k + 1;
//
//        Find L = pivot index
//
        l = IDAMAX(2, &MATA(k, k), 1) + k - 1;
        ARRAY1D(ipvt, k) = l;
//
//        Zero pivot implies this column already triangularized
//
        if (l == k) goto LABEL_10;
        t = MATA(l, k);
        MATA(l, k) = MATA(k, k);
        MATA(k, k) = t;
LABEL_10:
//
//           Compute multipliers
//
        t = -1.0 / MATA(k, k);
        MATA(k+1, k) *= t;
//
//           Row elimination with column indexing
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
// The old factorization of A will be updated.  A row and a column
// has been added to the matrix A.
// N-1 is now the old order of the matrix.
//
LABEL_80:
    nm1 = n - 1;
//
// Perform row interchanges on the elements of the new column, and
// perform elimination operations on the elements using the multipliers.
//
    if (nm1 <= 1) goto LABEL_105;
    for (k = 2; k <= nm1; ++k) {
        km1 = k -1;
        l = ARRAY1D(ipvt, km1);
        t = MATA(l, n);
        if (l == km1) goto LABEL_90;
        MATA(l, n) = MATA(km1, n);
        MATA(km1, n) = t;
LABEL_90:
        MATA(k, n) += MATA(k, km1) * t;
    }
LABEL_105:
//
// Complete update of factorization by decomposing last 2x2 block.
//
    info = 0;
//
//        Find L = pivot index
//
    l = IDAMAX(2, &MATA(nm1, nm1), 1) + nm1 - 1;
    ARRAY1D(ipvt, nm1) = l;
//
//        Zero pivot implies this column already triangularized
//
    if (MATA(l, nm1) == 0.0) goto LABEL_140;
//
//           Interchange if necessary
//
    if (l == nm1) goto LABEL_110;
    t = MATA(l, nm1);
    MATA(l, nm1) = MATA(nm1, nm1);
    MATA(nm1, nm1) = t;
LABEL_110:
//
//           Compute multipliers
//
    t = -1.0 / MATA(nm1, nm1);
    MATA(n, nm1) *= t;
//
//           Row elimination with column indexing
//
    t = MATA(l, n);
    if (l == nm1) goto LABEL_120;
    MATA(l, n) = MATA(nm1, n);
    MATA(nm1, n) = t;
LABEL_120:
    MATA(n, n) += t * MATA(n, nm1);
    goto LABEL_150;
LABEL_140:
    info = nm1;
LABEL_150:
    ARRAY1D(ipvt, n) = n;
    if (MATA(n, n) == 0.0) info = n;
    return;
//
#ifdef MATA
#undef MATA
#endif
}


/**
 * @fn DHESL
 * 
C-----------------------------------------------------------------------
C This is essentially the LINPACK routine DGESL except for changes
C due to the fact that A is an upper Hessenberg matrix.
C-----------------------------------------------------------------------
C     DHESL solves the real system A * x = b
C     using the factors computed by DHEFA.
C
C     On entry
C
C        A       odepack_cpp_real PRECISION(LDA, N)
C                the output from DHEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DHEFA.
C
C        B       odepack_cpp_real PRECISION(N)
C                the right hand side vector.
C
C     On return
C
C        B       the solution vector  x .
C
C     Modification of LINPACK, by Peter Brown, LLNL.
C     Written 7/20/83.  This version dated 6/20/01.
C
C     BLAS called: DAXPY
C-----------------------------------------------------------------------
 */
void Odepack::DHESL(odepack_cpp_real *a, int lda, int n, int *ipvt, odepack_cpp_real *b)
{
#ifndef MATA
#define MATA(i, j) ARRAY2D(a, lda, i, j)
#endif
//
    int k, kb, l, nm1;
    odepack_cpp_real t;
//
    nm1 = n - 1;
//
//        Solve  A * x = b
//        First solve  L*y = b
//
    if (nm1 < 1) goto LABEL_30;
    for (k = 1; k <= nm1; ++k) {
        l = ARRAY1D(ipvt, k);
        t = ARRAY1D(b, l);
        if (l == k) goto LABEL_10;
        ARRAY1D(b, l) = ARRAY1D(b, k);
        ARRAY1D(b, k) = t;
LABEL_10:
        ARRAY1D(b, k+1) += t * MATA(k+1, k);
    }
LABEL_30:
//
//        Now solve  U*x = y
//
    for (kb = 1; kb <= n; ++kb) {
        k = n + 1 - kb;
        ARRAY1D(b, k) /= MATA(k, k);
        t = - ARRAY1D(b, k);
        DAXPY(k-1, t, &MATA(1, k), 1, &ARRAY1D(b, 1), 1);
    }
    return;
//
#ifdef MATA
#undef MATA
#endif
}


/**
 * @fn DHEQR
 * 
C-----------------------------------------------------------------------
C     This routine performs a QR decomposition of an upper
C     Hessenberg matrix A.  There are two options available:
C
C          (1)  performing a fresh decomposition
C          (2)  updating the QR factors by adding a row and a
C               column to the matrix A.
C-----------------------------------------------------------------------
C     DHEQR decomposes an upper Hessenberg matrix by using Givens
C     rotations.
C
C     On entry
C
C        A       odepack_cpp_real PRECISION(LDA, N)
C                the matrix to be decomposed.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                A is an (N+1) by N Hessenberg matrix.
C
C        IJOB    INTEGER
C                = 1     means that a fresh decomposition of the
C                        matrix A is desired.
C                .ge. 2  means that the current decomposition of A
C                        will be updated by the addition of a row
C                        and a column.
C     On return
C
C        A       the upper triangular matrix R.
C                The factorization can be written Q*A = R, where
C                Q is a product of Givens rotations and R is upper
C                triangular.
C
C        Q       odepack_cpp_real PRECISION(2*N)
C                the factors c and s of each Givens rotation used
C                in decomposing A.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = k  if  A(k,k) .eq. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DHELS will divide by zero
C                     if called.
C
C     Modification of LINPACK, by Peter Brown, LLNL.
C     Written 1/13/86.  This version dated 6/20/01.
C-----------------------------------------------------------------------
 */
void Odepack::DHEQR(odepack_cpp_real *a, int lda, int n, odepack_cpp_real *q, int &info, int ijob)
{
#ifndef MATA
#define MATA(i, j) ARRAY2D(a, lda, i, j)
#endif
//
    int i, iq, j, k, km1, kp1, nm1;
    odepack_cpp_real c, s, t, t1, t2;
//
    if (ijob > 1) goto LABEL_70;
//
// A new facorization is desired.
//
//     QR decomposition without pivoting
//
    info = 0;
    for (k = 1; k <= n; ++k) {
        km1 = k - 1;
        kp1 = k + 1;
//
//           Compute kth column of R.
//           First, multiply the kth column of A by the previous
//           k-1 Givens rotations.
//
        if (km1 < 1) goto LABEL_20;
        for (j = 1; j <= km1; ++j) {
            i = 2 * (j - 1) + 1;
            t1 = MATA(j, k);
            t2 = MATA(j+1, k);
            c = ARRAY1D(q, i);
            s = ARRAY1D(q, i+1);
            MATA(j, k) = c * t1 - s * t2;
            MATA(j+1, k) = s * t1 + c * t2;
        }
//
//           Compute Givens components c and s
//
LABEL_20:
        iq = 2 * km1 + 1;
        t1 = MATA(k, k);
        t2 = MATA(kp1, k);
        if (t2 != 0.0) goto LABEL_30;
        c = 1.0;
        s = 0.0;
        goto LABEL_50;
LABEL_30:
        if (std::abs(t2) < std::abs(t1)) goto LABEL_40;
        t = t1 / t2;
        s = -1.0 / std::sqrt(1.0 + t * t);
        c = -s * t;
        goto LABEL_50;
LABEL_40:
        t = t2 / t1;
        c = -1.0 / std::sqrt(1.0 + t * t);
        s = -c * t;
LABEL_50:
        ARRAY1D(q, iq) = c;
        ARRAY1D(q, iq+1) = s;
        MATA(k, k) = c * t1 - s * t2;
        if (MATA(k, k) == 0.0) info = k;
    }
    return;
//
// The old factorization of A will be updated.  A row and a column
// has been added to the matrix A.
// N by N-1 is now the old size of the matrix.
//
LABEL_70:
    nm1 = n - 1;
//
// Multiply the new column by the N previous Givens rotations.
//    
    for (k = 1; k <= nm1; ++k) {
        i = 2 * (k - 1) + 1;
        t2 = MATA(k, n);
        t2 = MATA(k+1, n);
        c = ARRAY1D(q, i);
        s = ARRAY1D(q, i+1);
        MATA(k, n) = c * t1 - s * t2;
        MATA(k+1, n) = s * t1 + c * t2;
    }
//
// Complete update of decomposition by forming last Givens rotation,
// and multiplying it times the column vector (A(N,N), A(N+1,N)).
//
    info = 0;
    t1 = MATA(n, n);
    t2 = MATA(n+1, n);
    if (t2 != 0.0) goto LABEL_110;
    c = 1.0;
    s = 0.0;
    goto LABEL_130;
LABEL_110:
    if (std::abs(t2) < std::abs(t1)) goto LABEL_120;
    t = t1 / t2;
    s = -1.0 /  std::sqrt(1.0 + t * t);
    c = -s * t;
    goto LABEL_130;
LABEL_120:
    t = t2 / t1;
    c = 1.0 / std::sqrt(1.0 + t * t);
    s = -c * t;
LABEL_130:
    iq = 2 * n - 1;
    ARRAY1D(q, iq) = c;
    ARRAY1D(q, iq+1) = s;
    MATA(n, n) = c * t1 - s * t2;
    if (MATA(n, n) == 0.0) info = n;
    return;
//
#ifdef MATA
#undef MATA
#endif
}


/**
 * @fn DHELS
 * 
C-----------------------------------------------------------------------
C This is part of the LINPACK routine DGESL with changes
C due to the fact that A is an upper Hessenberg matrix.
C-----------------------------------------------------------------------
C     DHELS solves the least squares problem
C
C           min (b-A*x, b-A*x)
C
C     using the factors computed by DHEQR.
C
C     On entry
C
C        A       odepack_cpp_real PRECISION(LDA, N)
C                the output from DHEQR which contains the upper
C                triangular factor R in the QR decomposition of A.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                A is originally an (N+1) by N matrix.
C
C        Q       odepack_cpp_real PRECISION(2*N)
C                The coefficients of the N givens rotations
C                used in the QR factorization of A.
C
C        B       odepack_cpp_real PRECISION(N+1)
C                the right hand side vector.
C
C     On return
C
C        B       the solution vector  x .
C
C     Modification of LINPACK, by Peter Brown, LLNL.
C     Written 1/13/86.  This version dated 6/20/01.
C
C     BLAS called: DAXPY
C-----------------------------------------------------------------------
 */
void Odepack::DHELS(odepack_cpp_real *a, int lda, int n, odepack_cpp_real *q, odepack_cpp_real *b)
{
#ifndef MATA
#define MATA(i, j) ARRAY2D(a, lda, i, j)
#endif
//
    int iq, k, kb, kp1;
    odepack_cpp_real c, s, t, t1, t2;
//
//  Minimize (b-A*x, b-A*x)
//  First form Q*b.
//
    for (k = 1; k <= n; ++k) {
        kp1 = k + 1;
        iq = 2 * (k - 1) + 1;
        c = ARRAY1D(q, iq);
        s = ARRAY1D(q, iq+1);
        t1 = ARRAY1D(b, k);
        t2 = ARRAY1D(b, kp1);
        ARRAY1D(b, k) = c * t1 - s * t2;
        ARRAY1D(b, kp1) = s * t1 + c * t2;
    }
//
//  Now solve  R*x = Q*b.
//
    for (kb = 1; kb <= n; ++kb) {
        k = n + 1 - kb;
        ARRAY1D(b, k) /= MATA(k, k);
        t = -ARRAY1D(b, k);
        DAXPY(k-1, t, &MATA(1, k), 1, &ARRAY1D(b, 1), 1);
    }
    return;
//
#ifdef MATA
#undef MATA
#endif
}


/**
 * @fn DLHIN
 * 
C-----------------------------------------------------------------------
C Call sequence input -- NEQ, N, T0, Y0, YDOT, F, TOUT, UROUND,
C                        EWT, ITOL, ATOL, Y, TEMP
C Call sequence output -- H0, NITER, IER
C Common block variables accessed -- None
C
C Subroutines called by DLHIN: F, DCOPY
C Function routines called by DLHIN: DVNORM
C-----------------------------------------------------------------------
C This routine computes the step size, H0, to be attempted on the
C first step, when the user has not supplied a value for this.
C
C First we check that TOUT - T0 differs significantly from zero.  Then
C an iteration is done to approximate the initial second derivative
C and this is used to define H from WRMS-norm(H**2 * yddot / 2) = 1.
C A bias factor of 1/2 is applied to the resulting h.
C The sign of H0 is inferred from the initial values of TOUT and T0.
C
C Communication with DLHIN is done with the following variables:
C
C NEQ    = NEQ array of solver, passed to F.
C N      = size of ODE system, input.
C T0     = initial value of independent variable, input.
C Y0     = vector of initial conditions, input.
C YDOT   = vector of initial first derivatives, input.
C F      = name of subroutine for right-hand side f(t,y), input.
C TOUT   = first output value of independent variable
C UROUND = machine unit roundoff
C EWT, ITOL, ATOL = error weights and tolerance parameters
C                   as described in the driver routine, input.
C Y, TEMP = work arrays of length N.
C H0     = step size to be attempted, output.
C NITER  = number of iterations (and of f evaluations) to compute H0,
C          output.
C IER    = the error flag, returned with the value
C          IER = 0  if no trouble occurred, or
C          IER = -1 if TOUT and t0 are considered too close to proceed.
C-----------------------------------------------------------------------
 */
void Odepack::DLHIN(int neq, int n, odepack_cpp_real t0, odepack_cpp_real *y0, odepack_cpp_real *ydot, ODEPACK_FUNCTION f, odepack_cpp_real tout,
    odepack_cpp_real uround, odepack_cpp_real *ewt, int itol, odepack_cpp_real *atol, odepack_cpp_real *y, odepack_cpp_real *temp, odepack_cpp_real &h0, 
    int &niter, int &ier, void *user_data)
{
    odepack_cpp_real afi, atoli, delyi, hg, hlb, hnew, hrat, hub, t1, tdist, tround, yddnrm;
    int i, iter;
    odepack_cpp_real half = 0.5;
    odepack_cpp_real hun = 100.0;
    odepack_cpp_real pt1 = 0.1;
    odepack_cpp_real two = 2.0;
//
    niter = 0;
    tdist = std::abs(tout - t0);
    tround = uround * std::max(std::abs(t0), std::abs(tout));
    if (tdist < two * tround) goto LABEL_100;
//
// Set a lower bound on H based on the roundoff level in T0 and TOUT. ---
    hlb = hun * tround;
// Set an upper bound on H based on TOUT-T0 and the initial Y and YDOT. -
    hub = pt1 * tdist;
    atoli = ARRAY1D(atol, 1);
    for (i = 1; i <= n; ++i) {
        if (itol == 2 || itol == 4) atoli = ARRAY1D(atol, i);
        delyi = pt1 * std::abs(ARRAY1D(y0, i)) + atoli;
        afi = std::abs(ARRAY1D(ydot, i));
        if (afi * hub > delyi) hub = delyi / afi;
    }
//
// Set initial guess for H as geometric mean of upper and lower bounds. -
    iter = 0;
    hg = std::sqrt(hlb * hub);
// If the bounds have crossed, exit with the mean value. ----------------
    if (hub < hlb) {
        h0 = hg;
        goto LABEL_90;
    }
//
// Looping point for iteration. -----------------------------------------
LABEL_50:
// Estimate the second derivative as a difference quotient in f. --------
    t1 = t0 + hg;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = ARRAY1D(y0, i) + hg * ARRAY1D(ydot, i);
    }
    (*f)(neq, t1, y, temp, user_data);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(temp, i) = (ARRAY1D(temp, i) - ARRAY1D(ydot, i)) / hg;
    }
    yddnrm = DVNORM(n, temp, ewt);
// Get the corresponding new value of H. --------------------------------
    if (yddnrm * hub * hub > two) {
        hnew = std::sqrt(two / yddnrm);
    } else {
        hnew = std::sqrt(hg * hub);
    }
    iter++;
//-----------------------------------------------------------------------
// Test the stopping conditions.
// Stop if the new and previous H values differ by a factor of .lt. 2.
// Stop if four iterations have been done.  Also, stop with previous H
// if hnew/hg .gt. 2 after first iteration, as this probably means that
// the second derivative value is bad because of cancellation error.
//-----------------------------------------------------------------------
    if (iter >= 4) goto LABEL_80;
    hrat = hnew / hg;
    if ((hrat > half) && (hrat < two)) goto LABEL_80;
    if ((iter >= 2) && (hnew > two * hg)) {
        hnew = hg;
        goto LABEL_80;
    }
    hg = hnew;
    goto LABEL_50;
//
// Iteration done.  Apply bounds, bias factor, and sign. ----------------
LABEL_80:
    h0 = hnew * half;
    if (h0 < hlb) h0 = hlb;
    if (h0 > hub) h0 = hub;
LABEL_90:
    h0 = std::copysign(h0, tout - t0);
// Restore Y array from Y0, then exit. ----------------------------------
    DCOPY(n, y0, 1, y, 1);
    niter = iter;
    ier = 0;
    return;
// Error return for TOUT - T0 too small. --------------------------------
LABEL_100:
    ier = -1;
    return;
}


/**
 * @fn DSTOKA
 * 
C-----------------------------------------------------------------------
C DSTOKA performs one step of the integration of an initial value
C problem for a system of Ordinary Differential Equations.
C
C This routine was derived from Subroutine DSTODPK in the DLSODPK
C package by the addition of automatic functional/Newton iteration
C switching and logic for re-use of Jacobian data.
C-----------------------------------------------------------------------
C Note: DSTOKA is independent of the value of the iteration method
C indicator MITER, when this is .ne. 0, and hence is independent
C of the type of chord method used, or the Jacobian structure.
C Communication with DSTOKA is done with the following variables:
C
C NEQ    = integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls to F and JAC.
C Y      = an array of length .ge. N used as the Y argument in
C          all calls to F and JAC.
C YH     = an NYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C          j-th derivative of y(i), scaled by H**j/factorial(j)
C          (j = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C NYH    = a constant integer .ge. N, the first dimension of YH.
C YH1    = a one-dimensional array occupying the same space as YH.
C EWT    = an array of length N containing multiplicative weights
C          for local error measurements.  Local errors in y(i) are
C          compared to 1.0/EWT(i) in various error tests.
C SAVF   = an array of working storage, of length N.
C          Also used for input of YH(*,MAXORD+2) when JSTART = -1
C          and MAXORD .lt. the current order NQ.
C SAVX   = an array of working storage, of length N.
C ACOR   = a work array of length N, used for the accumulated
C          corrections.  On a successful return, ACOR(i) contains
C          the estimated one-step local error in y(i).
C WM,IWM = real and integer work arrays associated with matrix
C          operations in chord iteration (MITER .ne. 0).
C CCMAX  = maximum relative change in H*EL0 before DSETPK is called.
C H      = the step size to be attempted on the next step.
C          H is altered by the error control algorithm during the
C          problem.  H can be either positive or negative, but its
C          sign must remain constant throughout the problem.
C HMIN   = the minimum absolute value of the step size H to be used.
C HMXI   = inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C          HMIN and HMXI may be changed at any time, but will not
C          take effect until the next change of H is considered.
C TN     = the independent variable. TN is updated on each step taken.
C JSTART = an integer used for input only, with the following
C          values and meanings:
C               0  perform the first step.
C           .gt.0  take a new step continuing from the last.
C              -1  take the next step with a new value of H, MAXORD,
C                    N, METH, MITER, and/or matrix parameters.
C              -2  take the next step with a new value of H,
C                    but with other inputs unchanged.
C          On return, JSTART is set to 1 to facilitate continuation.
C KFLAG  = a completion code with the following meanings:
C               0  the step was succesful.
C              -1  the requested error could not be achieved.
C              -2  corrector convergence could not be achieved.
C              -3  fatal error in DSETPK or DSOLPK.
C          A return with KFLAG = -1 or -2 means either
C          ABS(H) = HMIN or 10 consecutive failures occurred.
C          On a return with KFLAG negative, the values of TN and
C          the YH array are as of the beginning of the last
C          step, and H is the last step size attempted.
C MAXORD = the maximum order of integration method to be allowed.
C MAXCOR = the maximum number of corrector iterations allowed.
C MSBP   = maximum number of steps between DSETPK calls (MITER .gt. 0).
C MXNCF  = maximum number of convergence failures allowed.
C METH/MITER = the method flags.  See description in driver.
C N      = the number of first-order differential equations.
C-----------------------------------------------------------------------
 */
void Odepack::DSTOKA(int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *yh1, odepack_cpp_real *ewt, odepack_cpp_real *savf, odepack_cpp_real *savx,
    odepack_cpp_real *acor, odepack_cpp_real *wm, int *iwm, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN3 jac, ODEPACK_PSOL psol, void *user_data)
{
#ifndef YH
#define YH(i, j) ARRAY2D(yh, nyh, i, j)
#endif
#ifndef ELCO
#define ELCO(i, j) ARRAY2D(elco, 13, i, j)
#endif
#ifndef TESCO
#define TESCO(i, j) ARRAY2D(tesco, 3, i, j)
#endif
// DLS001
    odepack_cpp_real &conit = dls1_.conit, &crate = dls1_.crate, *el = dls1_.el, *elco = dls1_.elco,
        &hold = dls1_.hold, &rmax = dls1_.rmax, *tesco = dls1_.tesco,
        &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &ialth = dls1_.ialth, &ipup = dls1_.ipup, &lmax = dls1_.lmax, &meo = dls1_.meo, &nqnyh = dls1_.nqnyh, &nslp = dls1_.nslp,
        &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
// DLS002
    odepack_cpp_real &stifr = dls2_.stifr;
    int &newt = dls2_.newt, &nsfi = dls2_.nsfi, &nslj = dls2_.nslj, &njev = dls2_.njev;
// DLPK01
    odepack_cpp_real &delt = dlpk_.delt, &epcon = dlpk_.epcon, &sqrtn = dlpk_.sqrtn, &rsqrtn = dlpk_.rsqrtn;
    int &jpre = dlpk_.jpre, &jacflg = dlpk_.jacflg, &locwp = dlpk_.locwp, &lociwp = dlpk_.lociwp, &lsavx = dlpk_.lsavx, &kmp = dlpk_.kmp, &maxl = dlpk_.maxl, &mnewt = dlpk_.mnewt,
        &nni = dlpk_.nni, &nli = dlpk_.nli, &nps = dlpk_.nps, &ncfn = dlpk_.ncfn, &ncfl = dlpk_.ncfl;
//
    int i, i1, iredo, iret, j, jb, jok, m, ncf, newq, nslow;
    odepack_cpp_real dcon, ddn, del, delp, drc, dsm, dup, exdn, exsm, exup, r, rh, rhdn, rhsm, rhup, roc, stiff, told, dfnorm;
//
    kflag = 0;
    told = tn;
    ncf = 0;
    ierpj = 0;
    iersl = 0;
    jcur = 0;
    icf = 0;
    delp = 0.0;
    if (jstart > 0) goto LABEL_200;
    if (jstart == -1) goto LABEL_100;
    if (jstart == -2) goto LABEL_160;
//-----------------------------------------------------------------------
// On the first call, the order is set to 1, and other variables are
// initialized.  RMAX is the maximum ratio by which H can be increased
// in a single step.  It is initially 1.E4 to compensate for the small
// initial H, but then is normally equal to 10.  If a failure
// occurs (in corrector convergence or error test), RMAX is set at 2
// for the next increase.
//-----------------------------------------------------------------------
    lmax = maxord + 1;
    nq = 1;
    l = 2;
    ialth = 2;
    rmax = 10000.0;
    rc = 0.0;
    el0 = 1.0;
    crate = 0.7;
    hold = h;
    meo = meth;
    nslp = 0;
    nslj = 0;
    ipup = 0;
    iret = 3;
    newt = 0;
    stifr = 0.0;
    goto LABEL_140;
//-----------------------------------------------------------------------
// The following block handles preliminaries needed when JSTART = -1.
// IPUP is set to MITER to force a matrix update.
// If an order increase is about to be considered (IALTH = 1),
// IALTH is reset to 2 to postpone consideration one more step.
// If the caller has changed METH, DCFODE is called to reset
// the coefficients of the method.
// If the caller has changed MAXORD to a value less than the current
// order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
// If H is to be changed, YH must be rescaled.
// If H or METH is being changed, IALTH is reset to L = NQ + 1
// to prevent further changes in H for that many steps.
//-----------------------------------------------------------------------
LABEL_100:
    ipup = miter;
    lmax = maxord + 1;
    if (ialth == 1) ialth = 2;
    if (meth == meo) goto LABEL_110;
    DCFODE(meth, elco, tesco);
    meo = meth;
    if (nq > maxord) goto LABEL_120;
    ialth = l;
    iret = 1;
    goto LABEL_150;
LABEL_110:
    if (nq <= maxord) goto LABEL_160;
LABEL_120:
    nq = maxord;
    l = lmax;
    for (i = 1; i <= l; ++i) {
        ARRAY1D(el, i) = ELCO(i, nq);
    }
    nqnyh = nq * nyh;
    rc *= ARRAY1D(el, 1) / el0;
    el0 = ARRAY1D(el, 1);
    conit = 0.5 / static_cast<odepack_cpp_real>(nq + 2);
    epcon = conit * TESCO(2, nq);
    ddn = DVNORM(n, savf, ewt) / TESCO(1, l);
    exdn = 1.0 / static_cast<odepack_cpp_real>(l);
    rhdn = 1.0 / (1.3 * std::pow(ddn, exdn) + 0.0000013);
    rh = std::min(rhdn, 1.0);
    h = hold;
    goto LABEL_175;
//-----------------------------------------------------------------------
// DCFODE is called to get all the integration coefficients for the
// current METH.  Then the EL vector and related constants are reset
// whenever the order NQ is changed, or at the start of the problem.
//-----------------------------------------------------------------------
LABEL_140:
    DCFODE(meth, elco, tesco);
LABEL_150:
    for (i = 1; i <= l; ++i) {
        ARRAY1D(el, i) = ELCO(i, nq);
    }
    nqnyh = nq * nyh;
    rc *= ARRAY1D(el, 1) / el0;
    el0 = ARRAY1D(el, 1);
    conit = 0.5 / static_cast<odepack_cpp_real>(nq + 2);
    epcon = conit * TESCO(2, nq);
    if (iret == 1) {
        goto LABEL_160;
    } else if (iret == 2) {
        goto LABEL_170;
    } else if (iret == 3) {
        goto LABEL_200;
    }
//-----------------------------------------------------------------------
// If H is being changed, the H ratio RH is checked against
// RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
// L = NQ + 1 to prevent a change of H for that many steps, unless
// forced by a convergence or error test failure.
//-----------------------------------------------------------------------
LABEL_160:
    if (h == hold) goto LABEL_200;
    rh = h / hold;
    h = hold;
    iredo = 3;
    goto LABEL_175;
LABEL_170:
    rh = std::max(rh, hmin / std::abs(h));
LABEL_175:
    rh /= std::max(1.0, std::abs(h) * hmxi * rh);
    r = 1.0;
    for (j = 2; j <= l; ++j) {
        r *= rh;
        for (i = 1; i <= n; ++i) {
            YH(i, j) *= r;
        }
    }
    h *= rh;
    rc *= rh;
    ialth = l;
    if (iredo == 0) goto LABEL_690;
//-----------------------------------------------------------------------
// This section computes the predicted values by effectively
// multiplying the YH array by the Pascal triangle matrix.
// The flag IPUP is set according to whether matrix data is involved
// (NEWT .gt. 0 .and. JACFLG .ne. 0) or not, to trigger a call to DSETPK.
// IPUP is set to MITER when RC differs from 1 by more than CCMAX,
// and at least every MSBP steps, when JACFLG = 1.
// RC is the ratio of new to old values of the coefficient  H*EL(1).
//-----------------------------------------------------------------------
LABEL_200:
    if (newt == 0 || jacflg == 0) {
        drc = 0.0;
        ipup = 0;
        crate = 0.7;
    } else {
        drc = std::abs(rc - 1.0);
        if (drc > ccmax) ipup = miter;
        if (nst >= nslp + msbp) ipup = miter;
    }
    tn += h;
    i1 = nqnyh + 1;
    for (jb = 1; jb <= nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= nqnyh; ++i) {
            ARRAY1D(yh1, i) += ARRAY1D(yh1, i + nyh);
        }
    }
//-----------------------------------------------------------------------
// Up to MAXCOR corrector iterations are taken.  A convergence test is
// made on the RMS-norm of each correction, weighted by the error
// weight vector EWT.  The sum of the corrections is accumulated in the
// vector ACOR(i).  The YH array is not altered in the corrector loop.
// Within the corrector loop, an estimated rate of convergence (ROC)
// and a stiffness ratio estimate (STIFF) are kept.  Corresponding
// global estimates are kept as CRATE and stifr.
//-----------------------------------------------------------------------
LABEL_220:
    m = 0;
    mnewt = 0;
    stiff = 0.0;
    roc = 0.05;
    nslow = 0;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = YH(i, 1);
    }
    (*f)(neq, tn, y, savf, user_data);
    nfe++;
    if (newt == 0 || ipup <= 0) goto LABEL_250;
//-----------------------------------------------------------------------
// If indicated, DSETPK is called to update any matrix data needed,
// before starting the corrector iteration.
// JOK is set to indicate if the matrix data need not be recomputed.
// IPUP is set to 0 as an indicator that the matrix data is up to date.
//-----------------------------------------------------------------------
    jok = 1;
    if (nst == 0 || nst > nslj + 50) jok = -1;
    if (icf == 1 && drc < 0.2) jok = -1;
    if (icf == 2) jok = -1;
    if (jok == -1) {
        nslj = nst;
        njev++;
    }
    DSETPK(neq, y, yh1, ewt, acor, savf, jok, wm, iwm, f, jac, user_data);
    ipup = 0;
    rc = 1.0;
    drc = 0.0;
    nslp = nst;
    crate = 0.7;
    if (ierpj != 0) goto LABEL_430;
LABEL_250:
    for (i = 1; i <= n; ++i) {
        ARRAY1D(acor, i) = 0.0;
    }
LABEL_270:
    if (newt != 0) goto LABEL_350;
//-----------------------------------------------------------------------
// In the case of functional iteration, update Y directly from
// the result of the last function evaluation, and STIFF is set to 1.0.
//-----------------------------------------------------------------------
    for (i = 1; i <= n; ++i) {
        ARRAY1D(savf, i) = h * ARRAY1D(savf, i) - YH(i, 2);
        ARRAY1D(y, i) = ARRAY1D(savf, i) - ARRAY1D(acor, i);
    }
    del = DVNORM(n, y, ewt);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = YH(i, 1) + ARRAY1D(el, i) * ARRAY1D(savf, i);
        ARRAY1D(acor, i) = ARRAY1D(savf, i);
    }
    stiff = 1.0;
    goto LABEL_400;
//-----------------------------------------------------------------------
// In the case of the chord method, compute the corrector error,
// and solve the linear system with that as right-hand side and
// P as coefficient matrix.  STIFF is set to the ratio of the norms
// of the residual and the correction vector.
//-----------------------------------------------------------------------
LABEL_350:
    for (i = 1; i <= n; ++i) {
        ARRAY1D(savx, i) = h * ARRAY1D(savf, i) - (YH(i, 2) + ARRAY1D(acor, i));
    }
    dfnorm = DVNORM(n, savx, ewt);
    DSOLPK(neq, y, savf, savx, ewt, wm, iwm, f, psol, user_data);
    if (iersl < 0) goto LABEL_430;
    if (iersl > 0) goto LABEL_410;
    del = DVNORM(n, savx, ewt);
    if (del > 1.0e-8) stiff = std::max(stiff, dfnorm / del);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(acor, i) += ARRAY1D(savx, i);
        ARRAY1D(y, i) = YH(i, 1) + ARRAY1D(el, i) * ARRAY1D(acor, i);
    }
//-----------------------------------------------------------------------
// Test for convergence.  If M .gt. 0, an estimate of the convergence
// rate constant is made for the iteration switch, and is also used
// in the convergence test.   If the iteration seems to be diverging or
// converging at a slow rate (.gt. 0.8 more than once), it is stopped.
//-----------------------------------------------------------------------
LABEL_400:
    if (m != 0) {
        roc = std::max(0.05, del / delp);
        crate = std::max(0.2 * crate, roc);
    }
    dcon = del * std::min(1.0, 1.5 * crate) / epcon;
    if (dcon <= 1.0) goto LABEL_450;
    m++;
    if (m == maxord) goto LABEL_410;
    if (m >= 2 && del > 2.0 * delp)  goto LABEL_410;
    if (roc > 10.0) goto LABEL_410;
    if (roc > 0.8) nslow++;
    if (nslow >= 2) goto LABEL_410;
    mnewt = m;
    delp = del;
    (*f)(neq, tn, y, savf, user_data);
    nfe++;
    goto LABEL_270;
//-----------------------------------------------------------------------
// The corrector iteration failed to converge.
// If functional iteration is being done (NEWT = 0) and MITER .gt. 0
// (and this is not the first step), then switch to Newton
// (NEWT = MITER), and retry the step.  (Setting STIFR = 1023 insures
// that a switch back will not occur for 10 step attempts.)
// If Newton iteration is being done, but using a preconditioner that
// is out of date (JACFLG .ne. 0 .and. JCUR = 0), then signal for a
// re-evalutation of the preconditioner, and retry the step.
// In all other cases, the YH array is retracted to its values
// before prediction, and H is reduced, if possible.  If H cannot be
// reduced or MXNCF failures have occurred, exit with KFLAG = -2.
//-----------------------------------------------------------------------
LABEL_410:
    icf = 1;
    if (newt == 0) {
        if (nst == 0) goto LABEL_430;
        if (miter == 0) goto LABEL_430;
        newt = miter;
        stifr = 1023.0;
        ipup = miter;
        goto LABEL_220;
    }
    if (jcur == 1 || jacflg == 0) goto LABEL_430;
    ipup = miter;
    goto LABEL_220;
LABEL_430:
    icf = 2;
    ncf++;
    ncfn++;
    rmax = 2.0;
    tn = told;
    i1 = nqnyh + 1;
    for (jb = 1; jb <= nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= nqnyh; ++i) {
            ARRAY1D(yh1, i) -= ARRAY1D(yh1, i + nyh);
        }
    }
    if (ierpj < 0 || iersl < 0) goto LABEL_680;
    if (std::abs(h) <= hmin * 1.00001) goto LABEL_670;
    if (ncf == mxncf) goto LABEL_670;
    rh = 0.5;
    ipup = miter;
    iredo = 1;
    goto LABEL_170;
//-----------------------------------------------------------------------
// The corrector has converged.  JCUR is set to 0 to signal that the
// preconditioner involved may need updating later.
// The stiffness ratio STIFR is updated using the latest STIFF value.
// The local error test is made and control passes to statement 500
// if it fails.
//-----------------------------------------------------------------------
LABEL_450:
    jcur = 0;
    if (newt > 0) stifr = 0.5 * (stifr + stiff);
    if (m == 0) dsm = del / TESCO(2, nq);
    if (m > 0) dsm = DVNORM(n, acor, ewt) / TESCO(2, nq);
    if (dsm > 1.0) goto LABEL_500;
//-----------------------------------------------------------------------
// After a successful step, update the YH array.
// If Newton iteration is being done and STIFR is less than 1.5,
// then switch to functional iteration.
// Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
// If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
// use in a possible order increase on the next step.
// If a change in H is considered, an increase or decrease in order
// by one is considered also.  A change in H is made only if it is by a
// factor of at least 1.1.  If not, IALTH is set to 3 to prevent
// testing for that many steps.
//-----------------------------------------------------------------------
    kflag = 0;
    iredo = 0;
    nst++;
    if (newt == 0) nsfi++;
    if (newt > 0 && stifr < 1.5) newt = 0;
    hu = h;
    nqu = nq;
    for (j = 1; j <= l; ++j) {
        for (i = 1; i <= n; ++i) {
            YH(i, j) += ARRAY1D(el, j) * ARRAY1D(acor, i);
        }
    }
    ialth--;
    if (ialth == 0) goto LABEL_520;
    if (ialth > 1) goto LABEL_700;
    if (l == lmax) goto LABEL_700;
    for (i = 1; i <= n; ++i) {
        YH(i, lmax) = ARRAY1D(acor, i);
    }
    goto LABEL_700;
//-----------------------------------------------------------------------
// The error test failed.  KFLAG keeps track of multiple failures.
// Restore TN and the YH array to their previous values, and prepare
// to try the step again.  Compute the optimum step size for this or
// one lower order.  After 2 or more failures, H is forced to decrease
// by a factor of 0.2 or less.
//-----------------------------------------------------------------------
LABEL_500:
    kflag--;
    tn = told;
    i1 = nqnyh + 1;
    for (jb = 1; jb <= nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= nqnyh; ++i) {
            ARRAY1D(yh1, i) -= ARRAY1D(yh1, i + nyh);
        }
    }
    rmax = 2.0;
    if (std::abs(h) <= hmin * 1.00001) goto LABEL_660;
    if (kflag <= -3) goto LABEL_640;
    iredo = 2;
    rhup = 0.0;
    goto LABEL_540;
//-----------------------------------------------------------------------
// Regardless of the success or failure of the step, factors
// RHDN, RHSM, and RHUP are computed, by which H could be multiplied
// at order NQ - 1, order NQ, or order NQ + 1, respectively.
// in the case of failure, RHUP = 0.0 to avoid an order increase.
// the largest of these is determined and the new order chosen
// accordingly.  If the order is to be increased, we compute one
// additional scaled derivative.
//-----------------------------------------------------------------------
LABEL_520:
    if (l == lmax) goto LABEL_540;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(savf, i) = ARRAY1D(acor, i) - YH(i, lmax);
    }
    dup = DVNORM(n, savf, ewt) / TESCO(3, nq);
    exup = 1.0 / static_cast<odepack_cpp_real>(l + 1);
    rhup = 1.0 / (1.4 * std::pow(dup, exup) + 0.0000014);
LABEL_540:
    exsm = 1.0 / static_cast<odepack_cpp_real>(l);
    rhsm = 1.0 / (1.2 * std::pow(dsm, exsm) + 0.0000012);
    rhdn = 0.0;
    if (nq == 1) goto LABEL_560;
    ddn = DVNORM(n, &YH(1, l), ewt) / TESCO(1, nq);
    exdn = 1.0 / static_cast<odepack_cpp_real>(nq);
    rhdn = 1.0 / (1.3 * std::pow(ddn, exdn) + 0.0000013);
LABEL_560:
    if (rhsm >= rhup) goto LABEL_570;
    if (rhup > rhdn) goto LABEL_590;
    goto LABEL_580;
LABEL_570:
    if (rhsm < rhdn) goto LABEL_580;
    newq = nq;
    rh = rhsm;
    goto LABEL_620;
LABEL_580:
    newq = nq - 1;
    rh = rhdn;
    if (kflag < 0 && rh > 1.0) rh = 1.0;
    goto LABEL_620;
LABEL_590:
    newq = l;
    rh = rhup;
    if (rh < 1.1) goto LABEL_610;
    r = ARRAY1D(el, l) / static_cast<odepack_cpp_real>(l);
    for (i = 1; i <= n; ++i) {
        YH(i, newq+1) = ARRAY1D(acor, i) * r;
    }
    goto LABEL_630;
LABEL_610:
    ialth = 3;
    goto LABEL_700;
LABEL_620:
    if ((kflag == 0) && (rh < 1.1)) goto LABEL_610;
    if (kflag <= -2) rh = std::min(rh, 0.2);
//-----------------------------------------------------------------------
// If there is a change of order, reset NQ, L, and the coefficients.
// In any case H is reset according to RH and the YH array is rescaled.
// Then exit from 690 if the step was OK, or redo the step otherwise.
//-----------------------------------------------------------------------
    if (newq == nq) goto LABEL_170;
LABEL_630:
    nq = newq;
    l = nq + 1;
    iret = 2;
    goto LABEL_150;
//-----------------------------------------------------------------------
// Control reaches this section if 3 or more failures have occured.
// If 10 failures have occurred, exit with KFLAG = -1.
// It is assumed that the derivatives that have accumulated in the
// YH array have errors of the wrong order.  Hence the first
// derivative is recomputed, and the order is set to 1.  Then
// H is reduced by a factor of 10, and the step is retried,
// until it succeeds or H reaches HMIN.
//-----------------------------------------------------------------------
LABEL_640:
    if (kflag == -10) goto LABEL_660;
    rh = 0.1;
    rh = std::max(hmin / std::abs(h), rh);
    h *= rh;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = YH(i, 1);
    }
    (*f)(neq, tn, y, savf, user_data);
    nfe++;
    for (i = 1; i <= n; ++i) {
        YH(i, 2) = h * ARRAY1D(savf, i);
    }
    ipup = miter;
    ialth = 5;
    if (nq == 1) goto LABEL_200;
    nq = 1;
    l = 2;
    iret = 3;
    goto LABEL_150;
//-----------------------------------------------------------------------
// All returns are made through this section.  H is saved in HOLD
// to allow the caller to change H on the next step.
//-----------------------------------------------------------------------
LABEL_660:
    kflag = -1;
    goto LABEL_720;
LABEL_670:
    kflag = -2;
    goto LABEL_720;
LABEL_680:
    kflag = -3;
    goto LABEL_720;
LABEL_690:
    rmax = 10.0;
LABEL_700:
    r = 1.0 / TESCO(2, nqu);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(acor, i) *= r;
    }
LABEL_720:
    hold = h;
    jstart = 1;
    return;
//
#ifdef YH
#undef YH
#endif
#ifdef ELCO
#undef ELCO
#endif
#ifdef TESCO
#undef TESCO
#endif
}


/**
 * @brief DSETPK
 * 
 * DSETPK is called by DSTOKA to interface with the user-supplied
 * routine JAC, to compute and process relevant parts of
 * the matrix P = I - H*EL(1)*J , where J is the Jacobian df/dy,
 * as need for preconditioning matrix operations later.
 */
void Odepack::DSETPK(int neq, odepack_cpp_real *y, odepack_cpp_real *ysv, odepack_cpp_real *ewt, odepack_cpp_real *ftem, odepack_cpp_real *savf, int jok, odepack_cpp_real *wm,int *iwm, 
    ODEPACK_FUNCTION f, ODEPACK_JACOBIAN3 jac, void *user_data)
{
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &init = dls1_.init, &mxstep = dls1_.mxstep, &mxhnil = dls1_.mxhnil, &nhnil = dls1_.nhnil, &nslast = dls1_.nslast, &nyh = dls1_.nyh, 
        &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
// DLPK01
    odepack_cpp_real &delt = dlpk_.delt, &epcon = dlpk_.epcon, &sqrtn = dlpk_.sqrtn, &rsqrtn = dlpk_.rsqrtn;
    int &jpre = dlpk_.jpre, &jacflg = dlpk_.jacflg, &locwp = dlpk_.locwp, &lociwp = dlpk_.lociwp, &lsavx = dlpk_.lsavx, &kmp = dlpk_.kmp, &maxl = dlpk_.maxl, &mnewt = dlpk_.mnewt,
        &nni = dlpk_.nni, &nli = dlpk_.nli, &nps = dlpk_.nps, &ncfn = dlpk_.ncfn, &ncfl = dlpk_.ncfl;
//
    int ier;
    odepack_cpp_real hl0;
//
    ierpj = 0;
    jcur = 0;
    if (jok == -1) jcur = 1;
    hl0 = el0 * h;
    // (*jac)(f, neq, tn, y, ysv, ewt, savf, ftem, hl0, jok, 
    //     &ARRAY1D(wm, locwp), &ARRAY1D(iwm, lociwp), ier, user_data);
    (*jac)(f, neq, tn, y, ysv, ewt, savf, ftem, hl0, 
        &ARRAY1D(wm, locwp), &ARRAY1D(iwm, lociwp), ier, user_data);
    nje++;
    if (ier == 0) return;
    ierpj = 1;
    return;
}


/**
 * @fn DSRCKR
 * 
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLS002, DLSR01, DLPK01, which
C are used internally by the DLSODKR solver.
C
C RSAV = real array of length 228 or more.
C ISAV = integer array of length 63 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
 */
void Odepack::DSRCKR(odepack_cpp_real *rsav, int *isav, int job)
{
// DLS001
    odepack_cpp_real *rls = dls1_.rls;
    int *ils = dls1_.ils;
// DLS002
    odepack_cpp_real &rls2 = dls2_.rls2;
    int *ils2 = dls2_.ils2;
// DLSR01
    odepack_cpp_real *rlsr = dlsr_.rlsr;
    int *ilsr = dlsr_.ilsr;
// DLPK01
    odepack_cpp_real *rlsp = dlpk_.rlsp;
    int *ilsp = dlpk_.ilsp;
//
    int lenrls = 218; // lenrls;
    int lenils = 37;  // lenils;
    int lenrlp = 4;
    int lenilp = 13;
    int lenrlr = 5;
    int lenilr = 9;
    int i, ioff;
//
    if (job == 2) goto LABEL_100;
    DCOPY(lenrls, rls, 1, rsav, 1);
    ARRAY1D(rsav, lenrls + 1) = rls2;
    DCOPY(lenrlr, rlsr, 1, &ARRAY1D(rsav, lenrls + 2), 1);
    DCOPY(lenrlp, rlsp, 1, &ARRAY1D(rsav, lenrls + lenrlr + 2), 1);
    for (i = 1; i <= lenils; ++i) {
        ARRAY1D(isav, i) = ARRAY1D(ils, i);
    }
    ARRAY1D(isav, lenils + 1) = ARRAY1D(ils2, 1);
    ARRAY1D(isav, lenils + 2) = ARRAY1D(ils2, 2);
    ARRAY1D(isav, lenils + 3) = ARRAY1D(ils2, 3);
    ARRAY1D(isav, lenils + 4) = ARRAY1D(ils2, 4);
    ioff = lenils + 2;
    for (i = 1; i <= lenilr; ++i) {
        ARRAY1D(isav, ioff + i) = ARRAY1D(ilsr, i); 
    }
    ioff += lenilr;
    for (i = 1; i <= lenilp; ++i) {
        ARRAY1D(isav, ioff + i) = ARRAY1D(ilsp, i);
    }
    return;
//
LABEL_100:
    DCOPY(lenrls, rsav, 1, rls, 1);
    rls2 = ARRAY1D(rsav, lenrls + 1);
    DCOPY(lenrlr, &ARRAY1D(rsav, lenrls + 2), 1, rlsr, 1);
    DCOPY(lenrlp, &ARRAY1D(rsav, lenrls + lenrlr + 2), 1, rlsp, 1);
    for (i = 1; i <= lenils; ++i) {
        ARRAY1D(ils, i) = ARRAY1D(isav, i);
    }
    ARRAY1D(ils2, 1) = ARRAY1D(isav, lenils + 1);
    ARRAY1D(ils2, 2) = ARRAY1D(isav, lenils + 2);
    ARRAY1D(ils2, 3) = ARRAY1D(isav, lenils + 3);
    ARRAY1D(ils2, 4) = ARRAY1D(isav, lenils + 4);
    ioff = lenils + 2;
    for (i = 1; i <= lenilr; ++i) {
        ARRAY1D(ilsr, i) = ARRAY1D(isav, ioff + i);
    }
    ioff += lenilr;
    for (i = 1; i <= lenilp; ++i) {
        ARRAY1D(ilsp, i) = ARRAY1D(isav, ioff + i);
    }
    return;
}


/**
 * @fn DAINVG
 * 
C-----------------------------------------------------------------------
C This subroutine computes the initial value
C of the vector YDOT satisfying
C     A * YDOT = g(t,y)
C when A is nonsingular.  It is called by DLSODI for
C initialization only, when ISTATE = 0 .
C DAINVG returns an error flag IER:
C   IER  =  0  means DAINVG was successful.
C   IER .ge. 2 means RES returned an error flag IRES = IER.
C   IER .lt. 0 means the a-matrix was found to be singular.
C-----------------------------------------------------------------------
 */
void Odepack::DAINVG(ODEPACK_RESIDUAL res, ODEPACK_ADDA1 adda, int neq, odepack_cpp_real t, odepack_cpp_real *y, odepack_cpp_real *ydot, int &miter, int &ml,
    int &mu, odepack_cpp_real *pw, int *ipvt, int &ier, void *user_data)
{
    int i, lenpw, mlp1, nrowpw;
//
    if (miter >= 4) goto LABEL_100;
//
// Full matrix case -----------------------------------------------------
//
    lenpw = neq * neq;
    for (i = 1; i <= lenpw; ++i) {
        ARRAY1D(pw, i) = 0.0;
    }
//
    ier = 1;
    (*res)(neq, t, y, pw, ydot, ier, user_data);
    if (ier > 1) return;
//
    (*adda)(neq, t, y, 0, 0, pw, neq, user_data);
    DGEFA(pw, neq, neq, ipvt, ier);
    if (ier == 0) goto LABEL_20;
    ier *= -1;
    return;
LABEL_20:
    DGESL(pw, neq, neq, ipvt, ydot, 0);
    return;
//
// Band matrix case -----------------------------------------------------
//
LABEL_100:
    nrowpw = 2 * ml + mu + 1;
    lenpw = neq * nrowpw;
    for (i = 1; i <= lenpw; ++i) {
        ARRAY1D(pw, i) = 0.0;
    }
//
    ier = 1;
    (*res)(neq, t, y, pw, ydot, ier, user_data);
    if (ier > 1) return;
//
    mlp1 = ml + 1;
    (*adda)(neq, t, y, ml, mu, &ARRAY1D(pw, mlp1), nrowpw, user_data);
    DGBSL(pw, nrowpw, neq, ml, mu, ipvt, ydot, 0);
    if (ier == 0) goto LABEL_120;
    ier *= -1;
    return;
LABEL_120:
    DGBSL(pw, nrowpw, neq, ml, mu, ipvt, ydot, 0);
    return;
}


/**
 * @fn DSTODI
 * 
C-----------------------------------------------------------------------
C DSTODI performs one step of the integration of an initial value
C problem for a system of Ordinary Differential Equations.
C Note: DSTODI is independent of the value of the iteration method
C indicator MITER, and hence is independent
C of the type of chord method used, or the Jacobian structure.
C Communication with DSTODI is done with the following variables:
C
C NEQ    = integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls to RES, ADDA,
C          and JAC.
C Y      = an array of length .ge. N used as the Y argument in
C          all calls to RES, JAC, and ADDA.
C NEQ    = integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls tO RES, G, ADDA,
C          and JAC.
C YH     = an NYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C          j-th derivative of y(i), scaled by H**j/factorial(j)
C          (j = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C NYH    = a constant integer .ge. N, the first dimension of YH.
C YH1    = a one-dimensional array occupying the same space as YH.
C EWT    = an array of length N containing multiplicative weights
C          for local error measurements.  Local errors in y(i) are
C          compared to 1.0/EWT(i) in various error tests.
C SAVF   = an array of working storage, of length N. also used for
C          input of YH(*,MAXORD+2) when JSTART = -1 and MAXORD is less
C          than the current order NQ.
C          Same as YDOTI in the driver.
C SAVR   = an array of working storage, of length N.
C ACOR   = a work array of length N used for the accumulated
C          corrections. On a succesful return, ACOR(i) contains
C          the estimated one-step local error in y(i).
C WM,IWM = real and integer work arrays associated with matrix
C          operations in chord iteration.
C PJAC   = name of routine to evaluate and preprocess Jacobian matrix.
C SLVS   = name of routine to solve linear system in chord iteration.
C CCMAX  = maximum relative change in H*EL0 before PJAC is called.
C H      = the step size to be attempted on the next step.
C          H is altered by the error control algorithm during the
C          problem.  H can be either positive or negative, but its
C          sign must remain constant throughout the problem.
C HMIN   = the minimum absolute value of the step size H to be used.
C HMXI   = inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C          HMIN and HMXI may be changed at any time, but will not
C          take effect until the next change of H is considered.
C TN     = the independent variable. TN is updated on each step taken.
C JSTART = an integer used for input only, with the following
C          values and meanings:
C               0  perform the first step.
C           .gt.0  take a new step continuing from the last.
C              -1  take the next step with a new value of H, MAXORD,
C                    N, METH, MITER, and/or matrix parameters.
C              -2  take the next step with a new value of H,
C                    but with other inputs unchanged.
C          On return, JSTART is set to 1 to facilitate continuation.
C KFLAG  = a completion code with the following meanings:
C               0  the step was succesful.
C              -1  the requested error could not be achieved.
C              -2  corrector convergence could not be achieved.
C              -3  RES ordered immediate return.
C              -4  error condition from RES could not be avoided.
C              -5  fatal error in PJAC or SLVS.
C          A return with KFLAG = -1, -2, or -4 means either
C          ABS(H) = HMIN or 10 consecutive failures occurred.
C          On a return with KFLAG negative, the values of TN and
C          the YH array are as of the beginning of the last
C          step, and H is the last step size attempted.
C MAXORD = the maximum order of integration method to be allowed.
C MAXCOR = the maximum number of corrector iterations allowed.
C MSBP   = maximum number of steps between PJAC calls.
C MXNCF  = maximum number of convergence failures allowed.
C METH/MITER = the method flags.  See description in driver.
C N      = the number of first-order differential equations.
C-----------------------------------------------------------------------
 */
template <typename ODEPACK_ADDA, typename ODEPACK_JACOBIAN>
void Odepack::DSTODI(
    int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *yh1, odepack_cpp_real *ewt, odepack_cpp_real *savf, odepack_cpp_real *savr,
    odepack_cpp_real *acor, odepack_cpp_real *wm, void *iwm_in, ODEPACK_RESIDUAL res, ODEPACK_ADDA adda, ODEPACK_JACOBIAN jac, 
    FUNC_PJAC2<ODEPACK_JACOBIAN, ODEPACK_ADDA> pjac, FUNC_SLVS slvs, void *user_data)
{
#ifndef YH
#define YH(i, j) ARRAY2D(yh, nyh, i, j)
#endif
#ifndef ELCO
#define ELCO(i, j) ARRAY2D(elco, 13, i, j)
#endif
#ifndef TESCO
#define TESCO(i, j) ARRAY2D(tesco, 3, i, j)
    #endif
// DLS001
    odepack_cpp_real &conit = dls1_.conit, &crate = dls1_.crate, *el = dls1_.el, *elco = dls1_.elco, 
        &hold = dls1_.hold, &rmax = dls1_.rmax, *tesco = dls1_.tesco,
        &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &ialth = dls1_.ialth, &ipup = dls1_.ipup, &lmax = dls1_.lmax, &meo = dls1_.meo, &nqnyh = dls1_.nqnyh, &nslp = dls1_.nslp,
        &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
//
    int i, i1, iredo, ires, iret, j, jb, kgo, m, ncf, newq;
    odepack_cpp_real dcon, ddn, del, delp, dsm, dup,
        eljh, el1h, exdn, exsm, exup,
        r, rh, rhdn, rhsm, rhup, told;
//
    int *iwm = static_cast<int*>(iwm_in);
//
    kflag = 0;
    told = tn;
    ncf = 0;
    ierpj = 0;
    iersl = 0;
    jcur = 0;
    icf = 0;
    delp = 0.0;
    if (jstart > 0) goto LABEL_200;
    if (jstart == -1) goto LABEL_100;
    if (jstart == -2) goto LABEL_160;
//-----------------------------------------------------------------------
// On the first call, the order is set to 1, and other variables are
// initialized.  RMAX is the maximum ratio by which H can be increased
// in a single step.  It is initially 1.E4 to compensate for the small
// initial H, but then is normally equal to 10.  If a failure
// occurs (in corrector convergence or error test), RMAX is set at 2
// for the next increase.
//-----------------------------------------------------------------------
    lmax = maxord + 1;
    nq = 1;
    l = 2;
    ialth = 2;
    rmax = 10000.0;
    rc = 0.0;
    el0 = 1.0;
    crate = 0.7;
    hold = h;
    meo = meth;
    nslp = 0;
    ipup = miter;
    iret = 3;
    goto LABEL_140;
//-----------------------------------------------------------------------
// The following block handles preliminaries needed when JSTART = -1.
// IPUP is set to MITER to force a matrix update.
// If an order increase is about to be considered (IALTH = 1),
// IALTH is reset to 2 to postpone consideration one more step.
// If the caller has changed METH, DCFODE is called to reset
// the coefficients of the method.
// If the caller has changed MAXORD to a value less than the current
// order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
// If H is to be changed, YH must be rescaled.
// If H or METH is being changed, IALTH is reset to L = NQ + 1
// to prevent further changes in H for that many steps.
//-----------------------------------------------------------------------
LABEL_100:
    ipup = miter;
    lmax = maxord + 1;
    if (ialth == 1) ialth = 2;
    if (meth == meo) goto LABEL_110;
    DCFODE(meth, elco, tesco);
    meo = meth;
    if (nq > maxord) goto LABEL_120;
    ialth = l;
    iret = 1;
    goto LABEL_150;
LABEL_110:
    if (nq <= maxord) goto LABEL_160;
LABEL_120:
    nq = maxord;
    l = lmax;
    for (i = 1; i <= l; ++i) {
        ARRAY1D(el, i) = ELCO(i, nq);
    }
    nqnyh = nq * nyh;
    rc *= ARRAY1D(el, 1) / el0;
    el0 = ARRAY1D(el, 1);
    conit = 0.5 / static_cast<odepack_cpp_real>(nq + 2);
    ddn = DVNORM(n, savf, ewt) / TESCO(1, l);
    exdn = 1.0 / static_cast<odepack_cpp_real>(l);
    rhdn = 1.0 / (1.3 * std::pow(ddn, exdn) + 0.0000013);
    rh = std::min(rhdn, 1.0);
    iredo = 3;
    if (h == hold) goto LABEL_170;
    rh = std::min(rh, std::abs(h / hold));
    h = hold;
    goto LABEL_175;
//-----------------------------------------------------------------------
// DCFODE is called to get all the integration coefficients for the
// current METH.  Then the EL vector and related constants are reset
// whenever the order NQ is changed, or at the start of the problem.
//-----------------------------------------------------------------------
LABEL_140:
    DCFODE(meth, elco, tesco);
LABEL_150:
    for (i = 1; i <= l; ++i) {
        ARRAY1D(el, i) = ELCO(i, nq);
    }
    nqnyh = nq * nyh;
    rc *= ARRAY1D(el, 1) / el0;
    el0 = ARRAY1D(el, 1);
    conit = 0.5 / static_cast<odepack_cpp_real>(nq + 2);
    if (iret == 1) {
        goto LABEL_160;
    } else if (iret == 2) {
        goto LABEL_170;
    } else if (iret == 3) {
        goto LABEL_200;
    }
//-----------------------------------------------------------------------
// If H is being changed, the H ratio RH is checked against
// RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
// L = NQ + 1 to prevent a change of H for that many steps, unless
// forced by a convergence or error test failure.
//-----------------------------------------------------------------------
LABEL_160:
    if (h == hold) goto LABEL_200;
    rh = h / hold;
    h = hold;
    iredo = 3;
    goto LABEL_175;
LABEL_170:
    rh = std::max(rh, hmin / std::abs(h));
LABEL_175:
    rh = std::min(rh, rmax);
    rh /= std::max(1.0, std::abs(h) * hmxi * rh);
    r = 1.0;
    for (j = 2; j <= l; ++j) {
        r *= rh;
        for (i = 1; i <= n; ++i) {
            YH(i, j) *= r;
        }
    }
    h *= rh;
    rc *= rh;
    ialth = l;
    if (iredo == 0) goto LABEL_690;
//-----------------------------------------------------------------------
// This section computes the predicted values by effectively
// multiplying the YH array by the Pascal triangle matrix.
// RC is the ratio of new to old values of the coefficient  H*EL(1).
// When RC differs from 1 by more than CCMAX, IPUP is set to MITER
// to force PJAC to be called.
// In any case, PJAC is called at least every MSBP steps.
//-----------------------------------------------------------------------
LABEL_200:
    if (std::abs(rc - 1.0) > ccmax) ipup = miter;
    if (nst >= nslp + msbp) ipup = miter;
    tn += h;
    i1 = nqnyh + 1;
    for (jb = 1; jb <= nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= nqnyh; ++i) {
            ARRAY1D(yh1, i) += ARRAY1D(yh1, i + nyh);
        }
    }
//-----------------------------------------------------------------------
// Up to MAXCOR corrector iterations are taken.  A convergence test is
// made on the RMS-norm of each correction, weighted by H and the
// error weight vector EWT.  The sum of the corrections is accumulated
// in ACOR(i).  The YH array is not altered in the corrector loop.
//-----------------------------------------------------------------------
LABEL_220:
    m = 0;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(savf, i) = YH(i, 2) / h;
        ARRAY1D(y, i) = YH(i, 1);
    }
    if (ipup <= 0) goto LABEL_240;
//-----------------------------------------------------------------------
// If indicated, the matrix P = A - H*EL(1)*dr/dy is reevaluated and
// preprocessed before starting the corrector iteration.  IPUP is set
// to 0 as an indicator that this has been done.
//-----------------------------------------------------------------------
    (this->*pjac)(neq, y, yh, nyh, ewt, acor, savr, savf, wm, iwm, 
        res, jac, adda, user_data);
    ipup = 0;
    rc = 1.0;
    nslp = nst;
    crate = 0.7;
    if (ierpj == 0) goto LABEL_250;
    if (iersl < 0) goto LABEL_435;
    ires = ierpj;
    if (ires == 1) {
        goto LABEL_430;
    } else if (ires == 2) {
        goto LABEL_435;
    } else if (ires == 3) {
        goto LABEL_430;
    }
// Get residual at predicted values, if not already done in PJAC. -------
LABEL_240:
    ires = 1;
    (*res)(neq, tn, y, savf, savr, ires, user_data);
    nfe++;
    kgo = std::abs(ires);
    if (kgo == 1) {
        goto LABEL_250;
    } else if (kgo == 2) {
        goto LABEL_435;
    } else if (kgo == 3) {
        goto LABEL_430;
    }
LABEL_250:
    for (i = 1; i <= n; ++i) {
        ARRAY1D(acor, i) = 0.0;
    }
//-----------------------------------------------------------------------
// Solve the linear system with the current residual as
// right-hand side and P as coefficient matrix.
//-----------------------------------------------------------------------
LABEL_270:
    (this->*slvs)(wm, iwm, savr, savf);
    if (iersl < 0) goto LABEL_430;
    if (iersl > 0) goto LABEL_410;
    el1h = ARRAY1D(el, 1) * h;
    del = DVNORM(n, savr, ewt) * std::abs(h);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(acor, i) += ARRAY1D(savr, i);
        ARRAY1D(savf, i) = ARRAY1D(acor, i) + YH(i, 2) / h;
        ARRAY1D(y, i) = YH(i, 1) + el1h * ARRAY1D(acor, i);
    }
//-----------------------------------------------------------------------
// Test for convergence.  If M .gt. 0, an estimate of the convergence
// rate constant is stored in CRATE, and this is used in the test.
//-----------------------------------------------------------------------
    if (m != 0) crate = std::max(0.2 * crate, del / delp);
    dcon = del * std::min(1.0, 1.5 * crate) / (TESCO(2, nq) * conit);
    if (dcon <= 1.0) goto LABEL_460;
    m++;
    if (m == maxcor) goto LABEL_410;
    if (m >= 2 && del > 2.0 * delp) goto LABEL_410;
    delp = del;
    ires = 1;
    (*res)(neq, tn, y, savf, savr, ires, user_data);
    nfe++;
    kgo = std::abs(ires);
    if (kgo == 1) {
        goto LABEL_270;
    } else if (kgo == 2) {
        goto LABEL_435;
    } else if (kgo == 3) {
        goto LABEL_410;
    }
//-----------------------------------------------------------------------
// The correctors failed to converge, or RES has returned abnormally.
// on a convergence failure, if the Jacobian is out of date, PJAC is
// called for the next try.  Otherwise the YH array is retracted to its
// values before prediction, and H is reduced, if possible.
// take an error exit if IRES = 2, or H cannot be reduced, or MXNCF
// failures have occurred, or a fatal error occurred in PJAC or SLVS.
//-----------------------------------------------------------------------
LABEL_410:
    icf = 1;
    if (jcur == 1) goto LABEL_430;
    ipup = miter;
    goto LABEL_220;
LABEL_430:
    icf = 2;
    ncf++;
    rmax = 2.0;
LABEL_435:
    tn = told;
    i1 = nqnyh + 1;
    for (jb = 1; jb <= nq; ++jb) {
        i1 -= nqnyh;
        for (i = i1; i <= nqnyh; ++i) {
            ARRAY1D(yh1, i) -= ARRAY1D(yh1, i + nyh);
        }
    }
    if (ires == 2) goto LABEL_680;
    if (ierpj < 0 || iersl < 0) goto LABEL_685;
    if (std::abs(h) <= hmin * 1.00001) goto LABEL_450;
    if (ncf == mxncf) goto LABEL_450;
    rh = 0.25;
    ipup = miter;
    iredo = 1;
    goto LABEL_170;
LABEL_450:
    if (ires == 3) goto LABEL_680;
    goto LABEL_670;
//-----------------------------------------------------------------------
// The corrector has converged.  JCUR is set to 0
// to signal that the Jacobian involved may need updating later.
// The local error test is made and control passes to statement 500
// if it fails.
//-----------------------------------------------------------------------
LABEL_460:
    jcur = 0;
    if (m == 0) dsm = del / TESCO(2, nq);
    if (m > 0) dsm = std::abs(h) * DVNORM(n, acor, ewt) / TESCO(2, nq);
    if (dsm > 1.0) goto LABEL_500;
//-----------------------------------------------------------------------
// After a successful step, update the YH array.
// Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
// If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
// use in a possible order increase on the next step.
// If a change in H is considered, an increase or decrease in order
// by one is considered also.  A change in H is made only if it is by a
// factor of at least 1.1.  If not, IALTH is set to 3 to prevent
// testing for that many steps.
//-----------------------------------------------------------------------
    kflag = 0;
    iredo = 0;
    nst++;
    hu = h;
    nqu = nq;
    for (j = 1; j <= l; ++j) {
        eljh = ARRAY1D(el, j) * h;
        for (i = 1; i <= n; ++i) {
            YH(i, j) += eljh * ARRAY1D(acor, i);
        }
    }
    ialth--;
    if (ialth == 0) goto LABEL_520;
    if (ialth > 1) goto LABEL_700;
    if (l == lmax) goto LABEL_700;
    for (i = 1; i <= n; ++i) {
        YH(i, lmax) = ARRAY1D(acor, i);
    }
    goto LABEL_700;
//-----------------------------------------------------------------------
// The error test failed.  KFLAG keeps track of multiple failures.
// restore TN and the YH array to their previous values, and prepare
// to try the step again.  Compute the optimum step size for this or
// one lower order.  After 2 or more failures, H is forced to decrease
// by a factor of 0.1 or less.
//-----------------------------------------------------------------------
LABEL_500:
    kflag--;
    tn = told;
    i1 = nqnyh + 1;
    for (jb = 1; jb <= nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= nqnyh; ++i) {
            ARRAY1D(yh1, i) -= ARRAY1D(yh1, i + nyh);
        }
    }
    rmax = 2.0;
    if (std::abs(h) <= hmin * 1.00001) goto LABEL_660;
    if (kflag <= -7) goto LABEL_660;
    iredo = 2;
    rhup = 0.0;
    goto LABEL_540;
//-----------------------------------------------------------------------
// Regardless of the success or failure of the step, factors
// RHDN, RHSM, and RHUP are computed, by which H could be multiplied
// at order NQ - 1, order NQ, or order NQ + 1, respectively.
// In the case of failure, RHUP = 0.0 to avoid an order increase.
// The largest of these is determined and the new order chosen
// accordingly.  If the order is to be increased, we compute one
// additional scaled derivative.
//-----------------------------------------------------------------------
LABEL_520:
    rhup = 0.0;
    if (l == lmax) goto LABEL_540;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(savf, i) = ARRAY1D(acor, i) - YH(i, lmax);
    }
    dup = std::abs(h) * DVNORM(n, savf, ewt) / TESCO(3, nq);
    exup = 1.0 / static_cast<odepack_cpp_real>(l + 1);
    rhup = 1.0 / (1.4 * std::pow(dup, exup) + 0.0000014);
LABEL_540:
    exsm = 1.0 / static_cast<odepack_cpp_real>(l);
    rhsm = 1.0 / (1.2 * std::pow(dsm, exsm) + 0.0000012);
    rhdn = 0.0;
    if (nq == 1) goto LABEL_560;
    ddn = DVNORM(n, &YH(1, l), ewt) / TESCO(1, nq);
    exdn = 1.0 / static_cast<odepack_cpp_real>(nq);
    rhdn = 1.0 / (1.3 * std::pow(ddn, exdn) + 0.0000013);
LABEL_560:
    if (rhsm >= rhup) goto LABEL_570;
    if (rhup > rhdn) goto LABEL_590;
    goto LABEL_580;
LABEL_570:
    if (rhsm < rhdn) goto LABEL_580;
    newq = nq;
    rh = rhsm;
    goto LABEL_620;
LABEL_580:
    newq = nq - 1;
    rh = rhdn;
    if (kflag < 0 && rh > 1.0) rh = 1.0;
    goto LABEL_620;
LABEL_590:
    newq = l;
    rh = rhup;
    if (rh < 1.1) goto LABEL_610;
    r = h * ARRAY1D(el, l) / static_cast<odepack_cpp_real>(l);
    for (i = 1; i <= n; ++i) {
        YH(i, newq + 1) = ARRAY1D(acor, i) * r;
    }
    goto LABEL_630;
LABEL_610:
    ialth = 3;
    goto LABEL_700;
LABEL_620:
    if (kflag == 0 && rh < 1.1) goto LABEL_610;
    if (kflag <= -2) rh = std::min(rh, 0.1);
//-----------------------------------------------------------------------
// If there is a change of order, reset NQ, L, and the coefficients.
// In any case H is reset according to RH and the YH array is rescaled.
// Then exit from 690 if the step was OK, or redo the step otherwise.
//-----------------------------------------------------------------------
    if (newq == nq) goto LABEL_170;
LABEL_630:
    nq = newq;
    l = nq + 1;
    iret = 2;
    goto LABEL_150;
//-----------------------------------------------------------------------
// All returns are made through this section.  H is saved in HOLD
// to allow the caller to change H on the next step.
//-----------------------------------------------------------------------
LABEL_660:
    kflag = -1;
    goto LABEL_720;
LABEL_670:
    kflag = -2;
    goto LABEL_720;
LABEL_680:
    kflag = -1 - ires;
    goto LABEL_720;
LABEL_685:
    kflag = -5;
    goto LABEL_720;
LABEL_690:
    rmax = 10.0;
LABEL_700:
    r = h / TESCO(2, nqu);
    for (i = 1; i <= n; ++i) {
        ARRAY1D(acor, i) *= r;
    }
LABEL_720:
    hold = h;
    jstart = 1;
    return;
//
#ifdef YH
#undef YH
#endif
#ifdef ELCO
#undef ELCO
#endif
#ifdef TESCO
#undef TESCO
#endif
}

// Explicit instantiation
template void Odepack::DSTODI<ODEPACK_ADDA1, ODEPACK_JACOBIAN4>(
    int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *yh1, odepack_cpp_real *ewt, odepack_cpp_real *savf, odepack_cpp_real *savr,
    odepack_cpp_real *acor, odepack_cpp_real *wm, void *iwm_in, ODEPACK_RESIDUAL res, ODEPACK_ADDA1 adda, ODEPACK_JACOBIAN4 jac, 
    FUNC_PJAC2<ODEPACK_JACOBIAN4, ODEPACK_ADDA1> pjac, FUNC_SLVS slvs, void *user_data
);

template void Odepack::DSTODI<ODEPACK_ADDA2, ODEPACK_JACOBIAN5>(
    int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *yh1, odepack_cpp_real *ewt, odepack_cpp_real *savf, odepack_cpp_real *savr,
    odepack_cpp_real *acor, odepack_cpp_real *wm, void *iwm_in, ODEPACK_RESIDUAL res, ODEPACK_ADDA2 adda, ODEPACK_JACOBIAN5 jac, 
    FUNC_PJAC2<ODEPACK_JACOBIAN5, ODEPACK_ADDA2> pjac, FUNC_SLVS slvs, void *user_data
);

template void Odepack::DSTODI<ODEPACK_ADDA3, ODEPACK_JACOBIAN6>(
    int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *yh1, odepack_cpp_real *ewt, odepack_cpp_real *savf, odepack_cpp_real *savr,
    odepack_cpp_real *acor, odepack_cpp_real *wm, void *iwm_in, ODEPACK_RESIDUAL res, ODEPACK_ADDA3 adda, ODEPACK_JACOBIAN6 jac, 
    FUNC_PJAC2<ODEPACK_JACOBIAN6, ODEPACK_ADDA3> pjac, FUNC_SLVS slvs, void *user_data
);


/**
 * @fn DPREPJI
 * 
C-----------------------------------------------------------------------
C DPREPJI is called by DSTODI to compute and process the matrix
C P = A - H*EL(1)*J , where J is an approximation to the Jacobian dr/dy,
C where r = g(t,y) - A(t,y)*s.  Here J is computed by the user-supplied
C routine JAC if MITER = 1 or 4, or by finite differencing if MITER =
C 2 or 5.  J is stored in WM, rescaled, and ADDA is called to generate
C P. P is then subjected to LU decomposition in preparation
C for later solution of linear systems with P as coefficient
C matrix.  This is done by DGEFA if MITER = 1 or 2, and by
C DGBFA if MITER = 4 or 5.
C
C In addition to variables described previously, communication
C with DPREPJI uses the following:
C Y     = array containing predicted values on entry.
C RTEM  = work array of length N (ACOR in DSTODI).
C SAVR  = array used for output only.  On output it contains the
C         residual evaluated at current values of t and y.
C S     = array containing predicted values of dy/dt (SAVF in DSTODI).
C WM    = real work space for matrices.  On output it contains the
C         LU decomposition of P.
C         Storage of matrix elements starts at WM(3).
C         WM also contains the following matrix-related data:
C         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
C IWM   = integer work space containing pivot information, starting at
C         IWM(21).  IWM also contains the band parameters
C         ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C EL0   = el(1) (input).
C IERPJ = output error flag.
C         = 0 if no trouble occurred,
C         = 1 if the P matrix was found to be singular,
C         = IRES (= 2 or 3) if RES returned IRES = 2 or 3.
C JCUR  = output flag = 1 to indicate that the Jacobian matrix
C         (or approximation) is now current.
C This routine also uses the Common variables EL0, H, TN, UROUND,
C MITER, N, NFE, and NJE.
C-----------------------------------------------------------------------
 */
// NEQ, Y, YH, NYH, EWT, RTEM, SAVR, S, WM, IWM, RES, JAC, ADDA
void Odepack::DPREPJI(
    int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *ewt, odepack_cpp_real *rtem, odepack_cpp_real *savr, odepack_cpp_real *s,
    odepack_cpp_real *wm, int *iwm, ODEPACK_RESIDUAL res, ODEPACK_JACOBIAN4 jac, ODEPACK_ADDA1 adda, void *user_data)
{
#ifndef YH
#define YH(i, j) ARRAY2D(yh, nyh, i, j)
#endif
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
//
    int i, i1, i2, ier, ii, ires, j, j1, jj, lenp,
        mba, mband, meb1, meband, ml, ml3, mu;
    odepack_cpp_real con, fac, hl0, r, srur, yi, yj, yjj;
//
    nje++;
    hl0 = h * el0;
    ierpj = 0;
    jcur = 1;
    if (miter == 1) {
        goto LABEL_100;
    } else if (miter == 2) {
        goto LABEL_200;
    } else if (miter == 3) {
        goto LABEL_300;
    } else if (miter == 4) {
        goto LABEL_400;
    } else if (miter == 5) {
        goto LABEL_500;
    }
// If MITER = 1, call RES, then JAC, and multiply by scalar. ------------
LABEL_100:
    ires = 1;
    (*res)(neq, tn, y, s, savr, ires, user_data);
    nfe++;
    if (ires > 1) goto LABEL_600;
    lenp = n * n;
    for (i = 1; i <= lenp; ++i) {
        ARRAY1D(wm, i + 2) = 0.0;
    }
    (*jac)(neq, tn, y, s, 0, 0, &ARRAY1D(wm, 3), n, user_data);
    con = -hl0;
    for (i = 1; i <= lenp; ++i) {
        ARRAY1D(wm, i + 2) *= con;
    }
    goto LABEL_240;
// If MITER = 2, make N + 1 calls to RES to approximate J. --------------
LABEL_200:
    ires = -1;
    (*res)(neq, tn, y, s, savr, ires, user_data);
    nfe++;
    if (ires > 1) goto LABEL_600;
    srur = ARRAY1D(wm, 1);
    j1 = 2;
    for (j = 1; j <= n; ++j) {
        yj = ARRAY1D(y, j);
        r = std::max(srur * std::abs(yj), 0.01 / ARRAY1D(ewt, j));
        ARRAY1D(y, j) += r;
        fac = - hl0 / r;
        (*res)(neq, tn, y, s, rtem, ires, user_data);
        nfe++;
        if (ires > 1) goto LABEL_600;
        for (i = 1; i <= n; ++i) {
            ARRAY1D(wm, i + j1) = (ARRAY1D(rtem, i) - ARRAY1D(savr, i)) * fac;
        }
        ARRAY1D(y, j) = yj;
        j1 += n;
    }
    ires = 1;
    (*res)(neq, tn, y, s, savr, ires, user_data);
    nfe++;
    if (ires > 1) goto LABEL_600;
// Add matrix A. --------------------------------------------------------
LABEL_240:
    (*adda)(neq, tn, y, 0, 0, &ARRAY1D(wm, 3), n, user_data);
// Do LU decomposition on P. --------------------------------------------
    DGEFA(&ARRAY1D(wm, 3), n, n, &ARRAY1D(iwm, 21), ier);
    if (ier != 0) ierpj = 1;
    return;
// Dummy section for MITER = 3
LABEL_300:
    return;
// If MITER = 4, call RES, then JAC, and multiply by scalar. ------------
LABEL_400:
    ires = 1;
    (*res)(neq, tn, y, s, savr, ires, user_data);
    nfe++;
    if (ires > 1) goto LABEL_600;
    ml = ARRAY1D(iwm, 1);
    mu = ARRAY1D(iwm, 2);
    ml3 = ml + 3;
    mband = ml + mu + 1;
    meband = mband + ml;
    lenp = meband * n;
    for (i = 1; i <= lenp; ++i) {
        ARRAY1D(wm, i + 2) = 0.0;
    }
    (*jac)(neq, tn, y, s, ml, mu, &ARRAY1D(wm, ml3), meband, user_data);
    con = -hl0;
    for (i = 1; i <= lenp; ++i) {
        ARRAY1D(wm, i + 2) *= con;
    }
    goto LABEL_570;
// If MITER = 5, make ML + MU + 2 calls to RES to approximate J. --------
LABEL_500:
    ires = -1;
    (*res)(neq, tn, y, s, savr, ires, user_data);
    nfe++;
    if (ires > 1) goto LABEL_600;
    ml = ARRAY1D(iwm, 1);
    mu = ARRAY1D(iwm, 2);
    ml3 = ml + 3;
    mband = ml + mu + 1;
    mba = std::min(mband, n);
    meband = mband + ml;
    meb1 = meband - 1;
    srur = ARRAY1D(wm, 1);
    for (j = 1; j <= mba; ++j) {
        for (i = j; i <= n; i += mband) {
            yi = ARRAY1D(y, i);
            r = std::max(srur * std::abs(yi), 0.01 / ARRAY1D(ewt, i));
            ARRAY1D(y, i) += r;
        }
        (*res)(neq, tn, y, s, rtem, ires, user_data);
        nfe++;
        if (ires > 1) goto LABEL_600;
        for (jj = j; jj <= n; jj += mband) {
            ARRAY1D(y, jj) = YH(jj, 1);
            yjj = ARRAY1D(y, jj);
            r = std::max(srur * std::abs(yjj), 0.01 / ARRAY1D(ewt, jj));
            fac = -hl0 / r;
            i1 = std::max(jj - mu, 1);
            i2 = std::min(jj + ml, n);
            ii = jj * meb1 - ml + 2;
            for (i = i1; i <= i2; ++i) {
                ARRAY1D(wm, ii + i) = (ARRAY1D(rtem, i) - ARRAY1D(savr, i)) * fac;
            }
        }
    }
    ires = 1;
    (*res)(neq, tn, y, s, savr, ires, user_data);
    nfe++;
    if (ires > 1) goto LABEL_600;
// Add matrix A. --------------------------------------------------------
LABEL_570:
    (*adda)(neq, tn, y, ml, mu, &ARRAY1D(wm, ml3), meband, user_data);
// Do LU decomposition of P. --------------------------------------------
    DGBFA(&ARRAY1D(wm, 3), meband, n, ml, mu, &ARRAY1D(iwm, 21), ier);
    if (ier != 0) ierpj = 1;
    return;
// Error return for IRES = 2 or IRES = 3 return from RES. ---------------
LABEL_600:
    ierpj = ires;
    return;
//
#ifdef YH
#undef YH
#endif
}


/**
 * @fn DAIGBT
 * 
C-----------------------------------------------------------------------
C This subroutine computes the initial value
C of the vector YDOT satisfying
C     A * YDOT = g(t,y)
C when A is nonsingular.  It is called by DLSOIBT for
C initialization only, when ISTATE = 0 .
C DAIGBT returns an error flag IER:
C   IER  =  0  means DAIGBT was successful.
C   IER .ge. 2 means RES returned an error flag IRES = IER.
C   IER .lt. 0 means the A matrix was found to have a singular
C              diagonal block (hence YDOT could not be solved for).
C-----------------------------------------------------------------------
 */
void Odepack::DAIGBT(ODEPACK_RESIDUAL res, ODEPACK_ADDA2 adda, int neq, odepack_cpp_real t, odepack_cpp_real *y, odepack_cpp_real *ydot, int &mb, int &nb,
    odepack_cpp_real *pw, int *ipvt, int &ier, void *user_data)
{
    int i, lenpw, lblox, lpb, lpc;
//
    lblox = mb * mb * nb;
    lpb = 1 + lblox;
    lpc = lpb + lblox;
    lenpw = 3 * lblox;
    for (i = 1; i <= lenpw; ++i) {
        ARRAY1D(pw, i) = 0.0;
    }
    ier = 1;
    (*res)(neq, t, y, pw, ydot, ier, user_data);
    if (ier > 1) return;
    (*adda)(neq, t, y, mb, nb, &ARRAY1D(pw, 1), &ARRAY1D(pw, lpb), &ARRAY1D(pw, lpc), user_data);
    DDECBT(mb, nb, pw, &ARRAY1D(pw, lpb), &ARRAY1D(pw, lpc), ipvt, ier);
    if (ier == 0) goto LABEL_20;
    ier *= -1;
    return;
LABEL_20:
    DSOLBT(mb, nb, pw, &ARRAY1D(pw, lpb), &ARRAY1D(pw, lpc), ydot, ipvt);
    return;
}


/**
 * @fn DPJIBT
 * 
C-----------------------------------------------------------------------
C DPJIBT is called by DSTODI to compute and process the matrix
C P = A - H*EL(1)*J , where J is an approximation to the Jacobian dr/dy,
C and r = g(t,y) - A(t,y)*s.  Here J is computed by the user-supplied
C routine JAC if MITER = 1, or by finite differencing if MITER = 2.
C J is stored in WM, rescaled, and ADDA is called to generate P.
C P is then subjected to LU decomposition by DDECBT in preparation
C for later solution of linear systems with P as coefficient matrix.
C
C In addition to variables described previously, communication
C with DPJIBT uses the following:
C Y     = array containing predicted values on entry.
C RTEM  = work array of length N (ACOR in DSTODI).
C SAVR  = array used for output only.  On output it contains the
C         residual evaluated at current values of t and y.
C S     = array containing predicted values of dy/dt (SAVF in DSTODI).
C WM    = real work space for matrices.  On output it contains the
C         LU decomposition of P.
C         Storage of matrix elements starts at WM(3).
C         WM also contains the following matrix-related data:
C         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
C IWM   = integer work space containing pivot information, starting at
C         IWM(21).  IWM also contains block structure parameters
C         MB = IWM(1) and NB = IWM(2).
C EL0   = EL(1) (input).
C IERPJ = output error flag.
C         = 0 if no trouble occurred,
C         = 1 if the P matrix was found to be unfactorable,
C         = IRES (= 2 or 3) if RES returned IRES = 2 or 3.
C JCUR  = output flag = 1 to indicate that the Jacobian matrix
C         (or approximation) is now current.
C This routine also uses the Common variables EL0, H, TN, UROUND,
C MITER, N, NFE, and NJE.
C-----------------------------------------------------------------------
 */
void Odepack::DPJIBT(int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *ewt, odepack_cpp_real *rtem, odepack_cpp_real *savr, odepack_cpp_real *s,
    odepack_cpp_real *wm, int *iwm, ODEPACK_RESIDUAL res, ODEPACK_JACOBIAN5 jac, ODEPACK_ADDA2 adda, void *user_data)
{
#ifndef YH
#define YH(i, j) ARRAY2D(yh, nyh, i, j)
#endif
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
//
    int i, ier, iia, iib, iic, ipa, ipb, ipc, ires, j, j1, j2,
        k, k1, lenp, lblox, lpb, lpc, mb, mbsq, mwid, nb;
    odepack_cpp_real con, fac, hl0, r, srur;
//
    nje++;
    hl0 = h * el0;
    ierpj = 0;
    jcur = 1;
    mb = ARRAY1D(iwm, 1);
    nb = ARRAY1D(iwm, 2);
    mbsq = mb * mb;
    lblox = mbsq * nb;
    lpb = 3 + lblox;
    lpc = lpb + lblox;
    lenp = 3 * lblox;
    if (miter == 1) {
        goto LABEL_100;
    } else if (miter == 2) {
        goto LABEL_200;
    }
// If MITER = 1, call RES, then JAC, and multiply by scalar. ------------
LABEL_100:
    ires = 1;
    (*res)(neq, tn, y, s, savr, ires, user_data);
    nfe++;
    if (ires > 1) goto LABEL_600;
    for (i = 1; i <= lenp; ++i) {
        ARRAY1D(wm, i + 2) = 0.0;
    }
    (*jac)(neq, tn, y, s, mb, nb, &ARRAY1D(wm, 3), &ARRAY1D(wm, lpb), &ARRAY1D(wm, lpc), user_data);
    con = -hl0;
    for (i = 1; i <= lenp; ++i) {
        ARRAY1D(wm, i + 2) *= con;
    }
    goto LABEL_260;
//
// If MITER = 2, make 3*MB + 1 calls to RES to approximate J. -----------
LABEL_200:
    ires = -1;
    (*res)(neq, tn, y, s, savr, ires, user_data);
    nfe++;
    if (ires > 1) goto LABEL_600;
    mwid = 3 * mb;
    srur = ARRAY1D(wm, 1);
    for (i = 1; i <= lenp; ++i) {
        ARRAY1D(wm, i + 2) = 0.0;
    }
    for (k = 1; k <= 3; ++k) {
        for (j = 1; j <= mb; ++j) {
//          Increment Y(I) for group of column indices, and call RES. ----
            j1 = j + (k - 1) * mb;
            for (i = j1; i <= n; i += mwid) {
                r = std::max(srur * std::abs(ARRAY1D(y, i)), 0.01 / ARRAY1D(ewt, i));
                ARRAY1D(y, i) += r;
            }
            (*res)(neq, tn, y, s, rtem, ires, user_data);
            nfe++;
            if (ires > 1) goto LABEL_600;
            for (i = 1; i <= n; ++i) {
                ARRAY1D(rtem, i) -= ARRAY1D(savr, i);
            }
            k1 = k;
            for (i = j1; i <= n; i += mwid) {
//              Get Jacobian elements in column I (block-column K1). -------
                ARRAY1D(y, i) = YH(i, 1);
                r = std::max(srur * std::abs(ARRAY1D(y, i)), 0.01 / ARRAY1D(ewt, i));
                fac = - hl0 / r;
//              Compute and load elements PA(*,J,K1). ----------------------
                iia = i - j;
                ipa = 2 + (j - 1) * mb + (k1 - 1) * mbsq;
                for (j2 = 1; j2 <= mb; ++j2) {
                    ARRAY1D(wm, ipa + j2) = ARRAY1D(rtem, iia + j2) * fac;
                }
                if (k1 <= 1) goto LABEL_223;
//              Compute and load elements PB(*,J,K1-1). --------------------
                iib = iia - mb;
                ipb = ipa + lblox - mbsq;
                for (j2 = 1; j2 <= mb; ++j2) {
                    ARRAY1D(wm, ipb + j2) = ARRAY1D(rtem, iib + j2) * fac;
                }
LABEL_223:
                if (k1 >= nb) goto LABEL_225;
//              Compute and load elements PC(*,J,K1+1). --------------------
                iic = iia + mb;
                ipc = ipa + 2 * lblox + mbsq;
                for (j2 = 1; j2 <= mb; ++j2) {
                    ARRAY1D(wm, ipc + j2) = ARRAY1D(rtem, iic + j2) * fac;
                }
LABEL_225:
                if (k1 != 3) goto LABEL_227;
//              Compute and load elements PC(*,J,1). -----------------------
                ipc = ipa - 2 * mbsq + 2 * lblox;
                for (j2 = 1; j2 <= mb; ++j2) {
                    ARRAY1D(wm, ipc + j2) = ARRAY1D(rtem, j2) * fac;
                }
LABEL_227:
                if (k1 != nb - 2) goto LABEL_229;
//              Compute and load elements PB(*,J,NB). ----------------------
                iib = n - mb;
                ipb = ipa + 2 * mbsq + lblox;
                for (j2 = 1; j2 <= mb; ++j2) {
                    ARRAY1D(wm, ipb + j2) = ARRAY1D(rtem, iib + j2) * fac;
                }
LABEL_229:
                k1 += 3;
            }
        }
    }
// RES call for first corrector iteration. ------------------------------
    ires = 1;
    (*res)(neq, tn, y, s, savr, ires, user_data);
    nfe++;
    if (ires > 1) goto LABEL_600;
// Add matrix A. --------------------------------------------------------
LABEL_260:
    (*adda)(neq, tn, y, mb, nb, &ARRAY1D(wm, 3), &ARRAY1D(wm, lpb), &ARRAY1D(wm, lpc), user_data);
// Do LU decomposition on P. --------------------------------------------
    DDECBT(mb, nb, &ARRAY1D(wm, 3), &ARRAY1D(wm, lpb), &ARRAY1D(wm, lpc), &ARRAY1D(iwm, 21), ier);
    if (ier != 0) ierpj = 1;
    return;
// Error return for IRES = 2 or IRES = 3 return from RES. ---------------
LABEL_600:
    ierpj = ires;
    return;
//
#ifdef YH
#undef YH
#endif
}


/**
 * @fn DSLSBT
 * 
C-----------------------------------------------------------------------
C This routine acts as an interface between the core integrator
C routine and the DSOLBT routine for the solution of the linear system
C arising from chord iteration.
C Communication with DSLSBT uses the following variables:
C WM    = real work space containing the LU decomposition,
C         starting at WM(3).
C IWM   = integer work space containing pivot information, starting at
C         IWM(21).  IWM also contains block structure parameters
C         MB = IWM(1) and NB = IWM(2).
C X     = the right-hand side vector on input, and the solution vector
C         on output, of length N.
C TEM   = vector of work space of length N, not used in this version.
C-----------------------------------------------------------------------
 */
void Odepack::DSLSBT(odepack_cpp_real *wm, int *iwm, odepack_cpp_real *x, odepack_cpp_real *tem)
{
    int lblox, lpb, lpc, mb, nb;
//
    mb = ARRAY1D(iwm, 1);
    nb = ARRAY1D(iwm, 2);
    lblox = mb * mb * nb;
    lpb = 3 + lblox;
    lpc = lpb + lblox;
    DSOLBT(mb, nb, &ARRAY1D(wm, 3), &ARRAY1D(wm, lpb), &ARRAY1D(wm, lpc), x, &ARRAY1D(iwm, 21));
    return;
}


/**
 * @fn DDECBT
 * 
C-----------------------------------------------------------------------
C Block-tridiagonal matrix decomposition routine.
C Written by A. C. Hindmarsh.
C Latest revision:  November 10, 1983 (ACH)
C Reference:  UCID-30150
C             Solution of Block-Tridiagonal Systems of Linear
C             Algebraic Equations
C             A.C. Hindmarsh
C             February 1977
C The input matrix contains three blocks of elements in each block-row,
C including blocks in the (1,3) and (N,N-2) block positions.
C DDECBT uses block Gauss elimination and Subroutines DGEFA and DGESL
C for solution of blocks.  Partial pivoting is done within
C block-rows only.
C
C Note: this version uses LINPACK routines DGEFA/DGESL instead of
C of dec/sol for solution of blocks, and it uses the BLAS routine DDOT
C for dot product calculations.
C
C Input:
C     M = order of each block.
C     N = number of blocks in each direction of the matrix.
C         N must be 4 or more.  The complete matrix has order M*N.
C     A = M by M by N array containing diagonal blocks.
C         A(i,j,k) contains the (i,j) element of the k-th block.
C     B = M by M by N array containing the super-diagonal blocks
C         (in B(*,*,k) for k = 1,...,N-1) and the block in the (N,N-2)
C         block position (in B(*,*,N)).
C     C = M by M by N array containing the subdiagonal blocks
C         (in C(*,*,k) for k = 2,3,...,N) and the block in the
C         (1,3) block position (in C(*,*,1)).
C    IP = integer array of length M*N for working storage.
C Output:
C A,B,C = M by M by N arrays containing the block-LU decomposition
C         of the input matrix.
C    IP = M by N array of pivot information.  IP(*,k) contains
C         information for the k-th digonal block.
C   IER = 0  if no trouble occurred, or
C       = -1 if the input value of M or N was illegal, or
C       = k  if a singular matrix was found in the k-th diagonal block.
C Use DSOLBT to solve the associated linear system.
C
C External routines required: DGEFA and DGESL (from LINPACK) and
C DDOT (from the BLAS, or Basic Linear Algebra package).
C-----------------------------------------------------------------------
 */
void Odepack::DDECBT(int m, int n, odepack_cpp_real *a, odepack_cpp_real *b, odepack_cpp_real *c, int *ip, int &ier)
{
#ifndef MAT3DA
#define MAT3DA(i, j, k) ARRAY3D(a, m, m, i, j, k)
#endif
#ifndef MAT3DB
#define MAT3DB(i, j, k) ARRAY3D(b, m, m, i, j, k)
#endif
#ifndef MAT3DC
#define MAT3DC(i, j, k) ARRAY3D(c, m, m, i, j, k)
#endif
#ifndef IP
#define IP(i, j) ARRAY2D(ip, m, i, j)
#endif
//
    int nm1, nm2, km1, i, j, k;
    odepack_cpp_real dp;
    if (m < 1 || n < 4) goto LABEL_210;
    nm1 = n - 1;
    nm2 = n - 2;
// Process the first block-row. -----------------------------------------
    DGEFA(a, m, m, ip, ier);
    k = 1;
    if (ier != 0) goto LABEL_200;
    for (j = 1; j <= m; ++j) {
        DGESL(a, m, m, ip, &MAT3DB(1, j, 1), 0);
        DGESL(a, m, m, ip, &MAT3DC(1, j, 1), 0);
    }
// Adjust B(*,*,2). -----------------------------------------------------
    for (j = 1; j <= m; ++j) {
        for (i = 1; i <= m; ++i) {
            dp = DDOT(m, &MAT3DC(i, 1, 2), m, &MAT3DC(1, j, 1), 1);
            MAT3DB(i, j, 2) -= dp;
        }
    }
// Main loop.  Process block-rows 2 to N-1. -----------------------------
    for (k = 2; k <= nm1; ++k) {
        km1 = k - 1;
        for (j = 1; j <= m; ++j) {
            for (i = 1; i <= m; ++i) {
                dp = DDOT(m, &MAT3DC(i, 1, k), m, &MAT3DB(1, j, km1), 1);
                MAT3DA(i, j, k) -= dp;
            }
        }
        DGEFA(&MAT3DA(1, 1, k), m, m, &IP(1, k), ier);
        if (ier != 0) goto LABEL_200;
        for (j = 1; j <= m; ++j) {
            DGESL(&MAT3DA(1, 1, k), m, m, &IP(1, k), &MAT3DB(1, j, k), 0);
        }
    }
// Process last block-row and return. -----------------------------------
    for (j = 1; j <= m; ++j) {
        for (i = 1; i <= m; ++i) {
            dp = DDOT(m, &MAT3DB(i, 1, n), m, &MAT3DB(1, j, nm2), 1);
            MAT3DC(i, j, n) -= dp;
        }
    }
    for (j = 1; j <= m; ++j) {
        for (i = 1; i <= m; ++i) {
            dp = DDOT(m, &MAT3DC(i, 1, n), m, &MAT3DB(1, j, nm1), 1);
            MAT3DA(i, j, n) -= dp;
        }
    }
    DGEFA(&MAT3DA(1, 1, n), m, m, &IP(1, n), ier);
    k = n;
    if (ier != 0) goto LABEL_200;
    return;
// Error returns. -------------------------------------------------------
LABEL_200:
    ier = k;
    return;
LABEL_210:
    ier = -1;
    return;
//
#ifdef MAT3DA
#undef MAT3DA
#endif
#ifdef MAT3DB
#undef MAT3DB
#endif
#ifdef MAT3DC
#undef MAT3DC
#endif
#ifdef IP
#undef IP
#endif
}


/**
 * @fn DSOLBT
 * 
C-----------------------------------------------------------------------
C Solution of block-tridiagonal linear system.
C Coefficient matrix must have been previously processed by DDECBT.
C M, N, A,B,C, and IP  must not have been changed since call to DDECBT.
C Written by A. C. Hindmarsh.
C Input:
C     M = order of each block.
C     N = number of blocks in each direction of matrix.
C A,B,C = M by M by N arrays containing block LU decomposition
C         of coefficient matrix from DDECBT.
C    IP = M by N integer array of pivot information from DDECBT.
C     Y = array of length M*N containg the right-hand side vector
C         (treated as an M by N array here).
C Output:
C     Y = solution vector, of length M*N.
C
C External routines required: DGESL (LINPACK) and DDOT (BLAS).
C-----------------------------------------------------------------------
 */
void Odepack::DSOLBT(int m, int n, odepack_cpp_real *a, odepack_cpp_real *b, odepack_cpp_real *c, odepack_cpp_real *y, int *ip)
{
#ifndef MAT3DA
#define MAT3DA(i, j, k) ARRAY3D(a, m, m, i, j, k)
#endif
#ifndef MAT3DB
#define MAT3DB(i, j, k) ARRAY3D(b, m, m, i, j, k)
#endif
#ifndef MAT3DC
#define MAT3DC(i, j, k) ARRAY3D(c, m, m, i, j, k)
#endif
#ifndef MATY
#define MATY(i, j) ARRAY2D(y, m, i, j)
#endif
#ifndef IP
#define IP(i, j) ARRAY2D(ip, m, i, j)
#endif
//
    int nm1, nm2, i, k, kb, km1, kp1;
    odepack_cpp_real dp;
    nm1 = n - 1;
    nm2 = n - 2;
// Forward solution sweep. ----------------------------------------------
    DGESL(a, m, m, ip, y, 0);
    for (k = 2; k <= nm1; ++k) {
        km1 = k - 1;
        for (i = 1; i <= m; ++i) {
            dp = DDOT(m, &MAT3DC(i, 1, k), m, &MATY(1, km1), 1);
            MATY(i, k) -= dp;
        }
        DGESL(&MAT3DA(1, 1, k), m, m, &IP(1, k), &MATY(1, k), 0);
    }
    for (i = 1; i <= m; ++i) {
        dp = DDOT(m, &MAT3DC(i, 1, n), m, &MATY(1, nm1), 1)
           + DDOT(m, &MAT3DB(i, 1, n), m, &MATY(1, nm2), 1);
        MATY(i, n) -= dp;
    }
    DGESL(&MAT3DA(1, 1, n), m, m, &IP(1, n), &MATY(1, n), 0);
// Backward solution sweep. ---------------------------------------------
    for (kb = 1; kb <= nm1; ++kb) {
        k = n - kb;
        kp1 = k + 1;
        for (i = 1; i <= m; ++i) {
            dp = DDOT(m, &MAT3DB(i, 1, k), m, &MATY(1, kp1), 1);
            MATY(i, k) -= dp;
        }
    }
    for (i = 1; i <= m; ++i) {
        dp = DDOT(m, &MAT3DC(i, 1, 1), m, &MATY(1, 3), 1);
        MATY(i, 1) -= dp;
    }
    return;
//
#ifdef MAT3DA
#undef MAT3DA
#endif
#ifdef MAT3DB
#undef MAT3DB
#endif
#ifdef MAT3DC
#undef MAT3DC
#endif
#ifdef MATY
#undef MATY
#endif
#ifdef IP
#undef IP
#endif
}


/**
 * @fn DIPREPI
 * 
C-----------------------------------------------------------------------
C This routine serves as an interface between the driver and
C Subroutine DPREPI.  Tasks performed here are:
C  * call DPREPI,
C  * reset the required WM segment length LENWK,
C  * move YH back to its final location (following WM in RWORK),
C  * reset pointers for YH, SAVR, EWT, and ACOR, and
C  * move EWT to its new position if ISTATE = 0 or 1.
C IPFLAG is an output error indication flag.  IPFLAG = 0 if there was
C no trouble, and IPFLAG is the value of the DPREPI error flag IPPER
C if there was trouble in Subroutine DPREPI.
C-----------------------------------------------------------------------
 */
void Odepack::DIPREPI(int neq, odepack_cpp_real *y, odepack_cpp_real *s, odepack_cpp_real *rwork, int *ia, int *ja, int *ic, int *jc, int &ipflag,
    ODEPACK_RESIDUAL res, ODEPACK_JACOBIAN6 jac, ODEPACK_ADDA3 adda, void *user_data)
{
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
// DLSS01
    int &iplost = dlss_.iplost, &iesp = dlss_.iesp, &istatc = dlss_.istatc, &iys = dlss_.iys, &iba = dlss_.iba, &ibian = dlss_.ibian, &ibjan = dlss_.ibjan, &ibjgp = dlss_.ibjgp,
        &ipian = dlss_.ipian, &ipjan = dlss_.ipjan, &ipjgp = dlss_.ipjgp, &ipigp = dlss_.ipigp, &ipr = dlss_.ipr, &ipc = dlss_.ipc, &ipic = dlss_.ipic, &ipisp = dlss_.ipisp, &iprsp = dlss_.iprsp, &ipa = dlss_.ipa,
        &lenyh = dlss_.lenyh, &lenyhm = dlss_.lenyhm, &lenwk = dlss_.lenwk, &lreq = dlss_.lreq, &lrat = dlss_.lrat, &lrest = dlss_.lrest, &lwmin = dlss_.lwmin, &moss = dlss_.moss, &msbj = dlss_.msbj,
        &nslj = dlss_.nslj, &ngp = dlss_.ngp, &nlu = dlss_.nlu, &nnz = dlss_.nnz, &nsp = dlss_.nsp, &nzl = dlss_.nzl, &nzu = dlss_.nzu;
//
    int i, imax, lentn, lyhd, lyhn;
//
    ipflag = 0;
// Call DPREPI to do matrix preprocessing operations. -------------------
    DPREPI(neq, y, s, &ARRAY1D(rwork, lyh), &ARRAY1D(rwork, lsavf), &ARRAY1D(rwork, lewt),
    &ARRAY1D(rwork, lacor), ia, ja, ic, jc, &ARRAY1D(rwork, lwm), &ARRAY1D(rwork, lwm), ipflag,
    res, jac, adda, user_data);
    lenwk = std::max(lreq, lwmin);
    if (ipflag < 0) return;
// If DPREPI was successful, move YH to end of required space for WM. ---
    lyhn = lwm + lenwk;
    if (lyhn > lyh) return;
    lyhd = lyh - lyhn;
    if (lyhd == 0) goto LABEL_20;
    imax = lyhn - 1 + lenyhm;
    for (i = lyhn; i <= imax; ++i) {
        ARRAY1D(rwork, i) = ARRAY1D(rwork, i + lyhd);
    }
    lyh = lyhn;
// Reset pointers for SAVR, EWT, and ACOR. ------------------------------
LABEL_20:
    lsavf = lyh + lenyh;
    lentn = lsavf + n;
    lacor = lentn + n;
    if (istatc == 3) goto LABEL_40;
// If ISTATE = 1, move EWT (left) to its new position. ------------------
    if (lentn > lewt) return;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(rwork, i + lentn - 1) = ARRAY1D(rwork, i + lewt - 1);
    }
LABEL_40:
    lewt = lentn;
    return;
}


/**
 * @fn DPREPI
 * 
C-----------------------------------------------------------------------
C This routine performs preprocessing related to the sparse linear
C systems that must be solved.
C The operations that are performed here are:
C  * compute sparseness structure of the iteration matrix
C      P = A - con*J  according to MOSS,
C  * compute grouping of column indices (MITER = 2),
C  * compute a new ordering of rows and columns of the matrix,
C  * reorder JA corresponding to the new ordering,
C  * perform a symbolic LU factorization of the matrix, and
C  * set pointers for segments of the IWK/WK array.
C In addition to variables described previously, DPREPI uses the
C following for communication:
C YH     = the history array.  Only the first column, containing the
C          current Y vector, is used.  Used only if MOSS .ne. 0.
C S      = array of length NEQ, identical to YDOTI in the driver, used
C          only if MOSS .ne. 0.
C SAVR   = a work array of length NEQ, used only if MOSS .ne. 0.
C EWT    = array of length NEQ containing (inverted) error weights.
C          Used only if MOSS = 2 or 4 or if ISTATE = MOSS = 1.
C RTEM   = a work array of length NEQ, identical to ACOR in the driver,
C          used only if MOSS = 2 or 4.
C WK     = a real work array of length LENWK, identical to WM in
C          the driver.
C IWK    = integer work array, assumed to occupy the same space as WK.
C LENWK  = the length of the work arrays WK and IWK.
C ISTATC = a copy of the driver input argument ISTATE (= 1 on the
C          first call, = 3 on a continuation call).
C IYS    = flag value from ODRV or CDRV.
C IPPER  = output error flag , with the following values and meanings:
C        =   0  no error.
C        =  -1  insufficient storage for internal structure pointers.
C        =  -2  insufficient storage for JGROUP.
C        =  -3  insufficient storage for ODRV.
C        =  -4  other error flag from ODRV (should never occur).
C        =  -5  insufficient storage for CDRV.
C        =  -6  other error flag from CDRV.
C        =  -7  if the RES routine returned error flag IRES = IER = 2.
C        =  -8  if the RES routine returned error flag IRES = IER = 3.
C-----------------------------------------------------------------------
 */
void Odepack::DPREPI(int neq, odepack_cpp_real *y, odepack_cpp_real *s, odepack_cpp_real *yh, odepack_cpp_real *savr, odepack_cpp_real* ewt, odepack_cpp_real *rtem,
    int *ia, int *ja, int *ic, int *jc, odepack_cpp_real *wk, void *iwk_in, int &ipper,
    ODEPACK_RESIDUAL res, ODEPACK_JACOBIAN6 jac, ODEPACK_ADDA3 adda, void *user_data)
{
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
// DLSS01
    int &iplost = dlss_.iplost, &iesp = dlss_.iesp, &istatc = dlss_.istatc, &iys = dlss_.iys, &iba = dlss_.iba, &ibian = dlss_.ibian, &ibjan = dlss_.ibjan, &ibjgp = dlss_.ibjgp,
        &ipian = dlss_.ipian, &ipjan = dlss_.ipjan, &ipjgp = dlss_.ipjgp, &ipigp = dlss_.ipigp, &ipr = dlss_.ipr, &ipc = dlss_.ipc, &ipic = dlss_.ipic, &ipisp = dlss_.ipisp, &iprsp = dlss_.iprsp, &ipa = dlss_.ipa,
        &lenyh = dlss_.lenyh, &lenyhm = dlss_.lenyhm, &lenwk = dlss_.lenwk, &lreq = dlss_.lreq, &lrat = dlss_.lrat, &lrest = dlss_.lrest, &lwmin = dlss_.lwmin, &moss = dlss_.moss, &msbj = dlss_.msbj,
        &nslj = dlss_.nslj, &ngp = dlss_.ngp, &nlu = dlss_.nlu, &nnz = dlss_.nnz, &nsp = dlss_.nsp, &nzl = dlss_.nzl, &nzu = dlss_.nzu;
//
    int i, ibr, ier, ipil, ipiu, iptt1, iptt2, j, k, knew, kamax,
        kamin, kcmax, kcmin, ldif, lenigp, lenwk1, liwk, ljfo, maxg,
        np1, nzsut;
    odepack_cpp_real erwt, fac, yj;
//
    int *iwk = static_cast<int*>(iwk_in);
//
    ibian = lrat * 2;
    ipian = ibian + 1;
    np1 = n + 1;
    ipjan = ipian + np1;
    ibjan = ipjan - 1;
    lenwk1 = lenwk - n;
    liwk = lenwk * lrat;
    if (moss == 0) liwk -= n;
    if (moss == 1 || moss == 2) liwk = lenwk1 * lrat;
    if (ipjan + n - 1 > liwk) goto LABEL_310;
    if (moss == 0) goto LABEL_30;
//
    if (istatc == 3) goto LABEL_20;
// ISTATE = 1 and MOSS .ne. 0.  Perturb Y for structure determination.
// Initialize S with random nonzero elements for structure determination.
    for (i = 1; i <= n; ++i) {
        erwt = 1.0 / ARRAY1D(ewt, i);
        fac = 1.0 + 1.0 / static_cast<odepack_cpp_real>(i + 1);
        ARRAY1D(y, i) += fac * std::copysign(erwt, ARRAY1D(y, i));
        ARRAY1D(s, i) = 1.0 + fac * erwt;
    }
    if (moss == 1) {
        goto LABEL_70;
    } else if (moss == 2) {
        goto LABEL_100;
    } else if (moss == 3) {
        goto LABEL_150;
    } else if (moss == 4) {
        goto LABEL_200;
    }
//
LABEL_20:
// ISTATE = 3 and MOSS .ne. 0. Load Y from YH(*,1) and S from YH(*,2). --
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = ARRAY1D(yh, i);
        ARRAY1D(s, i) = ARRAY1D(yh, n + i);
    }
    if (moss == 1) {
        goto LABEL_70;
    } else if (moss == 2) {
        goto LABEL_100;
    } else if (moss == 3) {
        goto LABEL_150;
    } else if (moss == 4) {
        goto LABEL_200;
    }
//
// MOSS = 0. Process user's IA,JA and IC,JC. ----------------------------
LABEL_30:
    knew = ipjan;
    kamin = ARRAY1D(ia, 1);
    kcmin = ARRAY1D(ic, 1);
    ARRAY1D(iwk, ipian) = 1;
    for (j = 1; j <= n; ++j) {
        for (i = 1; i <= n; ++i) {
            ARRAY1D(iwk, liwk + i) = 0;
        }
        kamax = ARRAY1D(ia, j + 1) - 1;
        if (kamin > kamax) goto LABEL_45;
        for (k = kamin; k <= kamax; ++k) {
            i = ARRAY1D(ja, k);
            ARRAY1D(iwk, liwk + i) = 1;
            if (knew > liwk) goto LABEL_310;
            ARRAY1D(iwk, knew) = i;
            knew++;
        }
LABEL_45:
        kamin = kamax + 1;
        kcmax = ARRAY1D(ic, j + 1) - 1;
        if (kcmin > kcmax) goto LABEL_55;
        for (k = kcmin; k <= kcmax; ++k) {
            i = ARRAY1D(jc, k);
            if (ARRAY1D(iwk, liwk + i) != 0) continue;
            if (knew > liwk) goto LABEL_310;
            ARRAY1D(iwk, knew) = i;
            knew++;
        }
LABEL_55:
        ARRAY1D(iwk, ipian + j) = knew + 1 - ipjan;
        kcmin = kcmax + 1;
    }
    goto LABEL_240;
//
// MOSS = 1. Compute structure from user-supplied Jacobian routine JAC. -
LABEL_70:
// A dummy call to RES allows user to create temporaries for use in JAC.
    ier = 1;
    (*res)(neq, tn, y, s, savr, ier, user_data);
    if (ier > 1) goto LABEL_370;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(savr, i) = 0.0;
        ARRAY1D(wk, lenwk1 + i) = 0.0;
    }
    k = ipjan;
    ARRAY1D(iwk, ipian) = 1;
    for (j = 1; j <= n; ++j) {
        (*adda)(neq, tn, y, j, &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan), &ARRAY1D(wk, lenwk1 + 1), user_data);
        (*jac)(neq, tn, y, s, j, &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan), savr, user_data);
        for (i = 1; i <= n; ++i) {
            ljfo = lenwk1 + i;
            if (ARRAY1D(wk, ljfo) == 0.0) goto LABEL_80;
            ARRAY1D(wk, ljfo) = 0.0;
            ARRAY1D(savr, i) = 0.0;
            goto LABEL_85;
LABEL_80:
            if (ARRAY1D(savr, i) == 0.0) continue;
            ARRAY1D(savr, i) = 0.0;
LABEL_85:
            if (k > liwk) goto LABEL_310;
            ARRAY1D(iwk, k) = i;
            k++;
        }
        ARRAY1D(iwk, ipian + j) = k + 1 - ipjan;
    }
    goto LABEL_240;
//
// MOSS = 2. Compute structure from results of N + 1 calls to RES. ------
LABEL_100:
    for (i = 1; i <= n; ++i) {
        ARRAY1D(wk, lenwk1 + i) = 0.0;
    }
    k = ipjan;
    ARRAY1D(iwk, ipian) = 1;
    ier = -1;
    if (miter == 1) ier = 1;
    (*res)(neq, tn, y, s, savr, ier, user_data);
    if (ier > 1) goto LABEL_370;
    for (j = 1; j <= n; ++j) {
        (*adda)(neq, tn, y, j, &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan), &ARRAY1D(wk, lenwk1 + 1), user_data);
        yj = ARRAY1D(y, j);
        erwt = 1.0 /  ARRAY1D(ewt, j);
        ARRAY1D(y, j) = yj + std::copysign(erwt, yj);
        (*res)(neq, tn, y, s, rtem, ier, user_data);
        if (ier > 1) return;
        ARRAY1D(y, j) = yj;
        for (i = 1; i <= n; ++i) {
            ljfo = lenwk1 + i;
            if (ARRAY1D(wk, ljfo) == 0.0) goto LABEL_110;
            ARRAY1D(wk, ljfo) = 0.0;
            goto LABEL_115;
LABEL_110:
            if (ARRAY1D(rtem, i) == ARRAY1D(savr, i)) continue;
LABEL_115:
            if (k > liwk) goto LABEL_310;
            ARRAY1D(iwk, k) = i;
            k++;
        }
        ARRAY1D(iwk, ipian + j) = k + 1 - ipjan;
    }
    goto LABEL_240;
//
// MOSS = 3. Compute structure from the user's IA/JA and JAC routine. ---
LABEL_150:
// A dummy call to RES allows user to create temporaries for use in JAC.
    ier = 1;
    (*res)(neq, tn, y, s, savr, ier, user_data);
    if (ier > 1) goto LABEL_370;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(savr, i) = 0.0;
    }
    knew = ipjan;
    kamin = ARRAY1D(ia, 1);
    ARRAY1D(iwk, ipian) = 1;
    for (j = 1; j <= n; ++j) {
        (*jac)(neq, tn, y, s, j, &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan), savr, user_data);
        kamax = ARRAY1D(ia, j + 1) - 1;
        if (kamin > kamax) goto LABEL_170;
        for (k = kamin; k <= kamax; ++k) {
            i = ARRAY1D(ja, k);
            ARRAY1D(savr, i) = 0.0;
            if (knew > liwk) goto LABEL_310;
            ARRAY1D(iwk, knew) = i;
            knew++;
        }
LABEL_170:
        kamin = kamax + 1;
        for (i = 1; i <= n; ++i) {
            if (ARRAY1D(savr, i) == 0.0) continue;
            ARRAY1D(savr, i) = 0.0;
            if (knew > liwk) goto LABEL_310;
            ARRAY1D(iwk, knew) = i;
            knew++;
        }
        ARRAY1D(iwk, ipian + j) = knew + 1 - ipjan;
    }
    goto LABEL_240;
//
// MOSS = 4. Compute structure from user's IA/JA and N + 1 RES calls. ---
LABEL_200:
    knew = ipjan;
    kamin = ARRAY1D(ia, 1);
    ARRAY1D(iwk, ipian) = 1;
    ier = -1;
    if (miter == 1) ier = 1;
    (*res)(neq, tn, y, s, savr, ier, user_data);
    if (ier > 1) goto LABEL_370;
    for (j = 1; j <= n; ++j) {
        yj = ARRAY1D(y, j);
        erwt = 1.0 / ARRAY1D(ewt, j);
        ARRAY1D(y, j) = yj + std::copysign(erwt, yj);
        (*res)(neq, tn, y, s, rtem, ier, user_data);
        if (ier > 1) return;
        ARRAY1D(y, j) = yj;
        kamax = ARRAY1D(ia, j + 1) - 1;
        if (kamin > kamax) goto LABEL_225;
        for (k = kamin; k <= kamax; ++k) {
            i = ARRAY1D(ja, k);
            ARRAY1D(rtem, i) = ARRAY1D(savr, i);
            if (knew > liwk) goto LABEL_310;
            ARRAY1D(iwk, knew) = i;
            knew++;
        }
LABEL_225:
        kamin = kamax + 1;
        for (i = 1; i <= n; ++i) {
            if (ARRAY1D(rtem, i) == ARRAY1D(savr, i)) continue;
            if (knew > liwk) goto LABEL_310;
            ARRAY1D(iwk, knew) = i;
            knew++;
        }
        ARRAY1D(iwk, ipian + j) = knew + 1 - ipjan;
    }
//
LABEL_240:
    if (moss == 0 || istatc == 3) goto LABEL_250;
// If ISTATE = 0 or 1 and MOSS .ne. 0, restore Y from YH. ---------------
    for (i = 1; i <= n; ++i) {
        ARRAY1D(y, i) = ARRAY1D(yh, i);
    }
LABEL_250:
    nnz = ARRAY1D(iwk, ipian + n) - 1;
    ipper = 0;
    ngp = 0;
    lenigp = 0;
    ipigp = ipjan + nnz;
    if (miter != 2) goto LABEL_260;
//
// Compute grouping of column indices (MITER = 2). ----------------------
//
    maxg = np1;
    ipjgp = ipjan + nnz;
    ibjgp = ipjgp - 1;
    ipigp = ipjgp + n;
    iptt1 = ipigp + np1;
    iptt2 = iptt1 + n;
    lreq = iptt2 + n - 1;
    if (lreq > liwk) goto LABEL_320;
    JGROUP(n, &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan), maxg, ngp, &ARRAY1D(iwk, ipigp),
        &ARRAY1D(iwk, ipjgp), &ARRAY1D(iwk, iptt1), &ARRAY1D(iwk, iptt2), ier);
    if (ier != 0) goto LABEL_320;
    lenigp = ngp + 1;
//
// Compute new ordering of rows/columns of Jacobian. --------------------
LABEL_260:
    ipr = ipigp + lenigp;
    ipc = ipr;
    ipic = ipc + n;
    ipisp = ipic + n;
    iprsp = (ipisp - 2) / lrat + 2;
    iesp = lenwk + 1 - iprsp;
    if (iesp < 0) goto LABEL_330;
    ibr = ipr - 1;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(iwk, ibr + i) = i;
    }
    nsp = liwk + 1 - ipisp;
    ODRV(n, &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan), wk, &ARRAY1D(iwk, ipr), &ARRAY1D(iwk, ipic), nsp,
        &ARRAY1D(iwk, ipisp), 1, iys);
    if (iys == 11 * n + 1) goto LABEL_340;
    if (iys != 0) goto LABEL_330;
//
// Reorder JAN and do symbolic LU factorization of matrix. --------------
    ipa = lenwk + 1 - nnz;
    nsp = ipa - iprsp;
    lreq = std::max(12 * n / lrat, 6 * n / lrat + 2 * n + nnz) + 3;
    lreq = lreq + iprsp - 1 + nnz;
    if (lreq > lenwk) goto LABEL_350;
    iba = ipa - 1;
    for (i = 1; i <= nnz; ++i) {
        ARRAY1D(wk, iba + i) = 0.0;
    }
    ipisp = lrat * (iprsp - 1) + 1;
    CDRV(n, &ARRAY1D(iwk, ipr), &ARRAY1D(iwk, ipc), &ARRAY1D(iwk, ipic), &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan),
        &ARRAY1D(wk, ipa), &ARRAY1D(wk, ipa), &ARRAY1D(wk, ipa), nsp, &ARRAY1D(iwk, ipisp), &ARRAY1D(wk, iprsp), iesp, 5, iys);
    lreq = lenwk - iesp;
    if (iys == 10 * n + 1) goto LABEL_350;
    if (iys != 0) goto LABEL_360;
    ipil = ipisp;
    ipiu = ipil + 2 * n + 1;
    nzu = ARRAY1D(iwk, ipil + n) - ARRAY1D(iwk, ipil);
    nzl = ARRAY1D(iwk, ipiu + n) - ARRAY1D(iwk, ipiu);
    if (lrat > 1) goto LABEL_290;
    ADJLR(n, &ARRAY1D(iwk, ipisp), ldif);
    lreq += ldif;
LABEL_290:
    if (lrat == 2 && nnz == n) lreq += 1;
    nsp = nsp + lreq - lenwk;
    ipa = lreq + 1 - nnz;
    iba = ipa - 1;
    ipper = 0;
    return;
//
LABEL_310:
    ipper = -1;
    lreq = 2 + (2 * n + 1) / lrat;
    lreq = std::max(lenwk + 1, lreq);
    return;
//
LABEL_320:
    ipper = -2;
    lreq = (lreq - 1) / lrat + 1;
    return;
//
LABEL_330:
    ipper = -3;
    CNTNZU(n, &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan), nzsut);
    lreq = lenwk - iesp + (3 * n + 4 * nzsut - 1) / lrat + 1;
    return;
//
LABEL_340:
    ipper = -4;
    return;
//
LABEL_350:
    ipper = -5;
    return;
//
LABEL_360:
    ipper = -6;
    lreq = lenwk;
    return;
//
LABEL_370:
    ipper = -ier - 5;
    lreq = 2 + (2 * n + 1) / lrat;
    return;
}


/**
 * @fn DAINVGS
 * 
C-----------------------------------------------------------------------
C This subroutine computes the initial value of the vector YDOT
C satisfying
C     A * YDOT = g(t,y)
C when A is nonsingular.  It is called by DLSODIS for initialization
C only, when ISTATE = 0.  The matrix A is subjected to LU
C decomposition in CDRV.  Then the system A*YDOT = g(t,y) is solved
C in CDRV.
C In addition to variables described previously, communication
C with DAINVGS uses the following:
C Y     = array of initial values.
C WK    = real work space for matrices.  On output it contains A and
C         its LU decomposition.  The LU decomposition is not entirely
C         sparse unless the structure of the matrix A is identical to
C         the structure of the Jacobian matrix dr/dy.
C         Storage of matrix elements starts at WK(3).
C         WK(1) = SQRT(UROUND), not used here.
C IWK   = integer work space for matrix-related data, assumed to
C         be equivalenced to WK.  In addition, WK(IPRSP) and WK(IPISP)
C         are assumed to have identical locations.
C TEM   = vector of work space of length N (ACOR in DSTODI).
C YDOT  = output vector containing the initial dy/dt. YDOT(i) contains
C         dy(i)/dt when the matrix A is non-singular.
C IER   = output error flag with the following values and meanings:
C       = 0  if DAINVGS was successful.
C       = 1  if the A-matrix was found to be singular.
C       = 2  if RES returned an error flag IRES = IER = 2.
C       = 3  if RES returned an error flag IRES = IER = 3.
C       = 4  if insufficient storage for CDRV (should not occur here).
C       = 5  if other error found in CDRV (should not occur here).
C-----------------------------------------------------------------------
 */
void Odepack::DAINVGS(int neq, odepack_cpp_real t, odepack_cpp_real *y, odepack_cpp_real *wk, void *iwk_in, odepack_cpp_real *tem, odepack_cpp_real *ydot,int &ier,
    ODEPACK_RESIDUAL res, ODEPACK_ADDA3 adda, void *user_data)
{
    int i, imul, j, k, kmin, kmax;
    odepack_cpp_real rlss;
// DLSS01
    int &iplost = dlss_.iplost, &iesp = dlss_.iesp, &istatc = dlss_.istatc, &iys = dlss_.iys, &iba = dlss_.iba, &ibian = dlss_.ibian, &ibjan = dlss_.ibjan, &ibjgp = dlss_.ibjgp,
        &ipian = dlss_.ipian, &ipjan = dlss_.ipjan, &ipjgp = dlss_.ipjgp, &ipigp = dlss_.ipigp, &ipr = dlss_.ipr, &ipc = dlss_.ipc, &ipic = dlss_.ipic, &ipisp = dlss_.ipisp, &iprsp = dlss_.iprsp, &ipa = dlss_.ipa,
        &lenyh = dlss_.lenyh, &lenyhm = dlss_.lenyhm, &lenwk = dlss_.lenwk, &lreq = dlss_.lreq, &lrat = dlss_.lrat, &lrest = dlss_.lrest, &lwmin = dlss_.lwmin, &moss = dlss_.moss, &msbj = dlss_.msbj,
        &nslj = dlss_.nslj, &ngp = dlss_.ngp, &nlu = dlss_.nlu, &nnz = dlss_.nnz, &nsp = dlss_.nsp, &nzl = dlss_.nzl, &nzu = dlss_.nzu;
//
    int *iwk = static_cast<int*>(iwk_in);
//
    for (i = 1; i <= nnz; ++i) {
        ARRAY1D(wk, iba + i) = 0.0;
    }
//
    ier = 1;
    (*res)(neq, t, y, &ARRAY1D(wk, ipa), ydot, ier, user_data);
    if (ier > 1) return;
//
    kmin = ARRAY1D(iwk, ipian);
    for (j = 1; j <= neq; ++j) {
        kmax = ARRAY1D(iwk, ipian + j) - 1;
        for (k = kmin; k <= kmax; ++k) {
            i = ARRAY1D(iwk, ibjan + k);
            ARRAY1D(tem, i) = 0.0;
        }
        (*adda)(neq, t, y, j, &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan), tem, user_data);
        for (k = kmin; k <= kmax; ++k) {
            i = ARRAY1D(iwk, ibjan + k);
            ARRAY1D(wk, iba + k) = ARRAY1D(tem, i);
        }
        kmin = kmax + 1;
    }
    nlu++;
    ier = 0;
    for (i = 1; i <= neq; ++i) {
        ARRAY1D(tem, i) = 0.0;
    }
//
// Numerical factorization of matrix A. ---------------------------------
    CDRV(neq, &ARRAY1D(iwk, ipr), &ARRAY1D(iwk, ipc), &ARRAY1D(iwk, ipic), &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan),
        &ARRAY1D(wk, ipa), tem, tem, nsp, &ARRAY1D(iwk, ipisp), &ARRAY1D(wk, iprsp), iesp, 2, iys);
    if (iys == 0) goto LABEL_50;
    imul = (iys - 1) / neq;
    ier = 5;
    if (imul == 8) ier = 1;
    if (imul == 10) ier = 4;
    return;
//
// Solution of the linear system. ---------------------------------------
LABEL_50:    
    CDRV(neq, &ARRAY1D(iwk, ipr), &ARRAY1D(iwk, ipc), &ARRAY1D(iwk, ipic), &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan),
            &ARRAY1D(wk, ipa), ydot, ydot, nsp, &ARRAY1D(iwk, ipisp), &ARRAY1D(wk, iprsp), iesp, 4, iys);
    if (iys != 0) ier = 5;
    return;
}


/**
 * @fn DPRJIS
 * 
C-----------------------------------------------------------------------
C DPRJIS is called to compute and process the matrix
C P = A - H*EL(1)*J, where J is an approximation to the Jacobian dr/dy,
C where r = g(t,y) - A(t,y)*s.  J is computed by columns, either by
C the user-supplied routine JAC if MITER = 1, or by finite differencing
C if MITER = 2.  J is stored in WK, rescaled, and ADDA is called to
C generate P.  The matrix P is subjected to LU decomposition in CDRV.
C P and its LU decomposition are stored separately in WK.
C
C In addition to variables described previously, communication
C with DPRJIS uses the following:
C Y     = array containing predicted values on entry.
C RTEM  = work array of length N (ACOR in DSTODI).
C SAVR  = array containing r evaluated at predicted y. On output it
C         contains the residual evaluated at current values of t and y.
C S     = array containing predicted values of dy/dt (SAVF in DSTODI).
C WK    = real work space for matrices.  On output it contains P and
C         its sparse LU decomposition.  Storage of matrix elements
C         starts at WK(3).
C         WK also contains the following matrix-related data.
C         WK(1) = SQRT(UROUND), used in numerical Jacobian increments.
C IWK   = integer work space for matrix-related data, assumed to be
C         equivalenced to WK.  In addition,  WK(IPRSP) and IWK(IPISP)
C         are assumed to have identical locations.
C EL0   = EL(1) (input).
C IERPJ = output error flag (in COMMON).
C         =  0 if no error.
C         =  1 if zero pivot found in CDRV.
C         = IRES (= 2 or 3) if RES returned IRES = 2 or 3.
C         = -1 if insufficient storage for CDRV (should not occur).
C         = -2 if other error found in CDRV (should not occur here).
C JCUR  = output flag = 1 to indicate that the Jacobian matrix
C         (or approximation) is now current.
C This routine also uses other variables in Common.
C-----------------------------------------------------------------------
 */
void Odepack::DPRJIS(int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *ewt, odepack_cpp_real *rtem, odepack_cpp_real *savr,
    odepack_cpp_real *s, odepack_cpp_real *wk, int *iwk, ODEPACK_RESIDUAL res, ODEPACK_JACOBIAN6 jac, ODEPACK_ADDA3 adda, void *user_data)
{
#ifndef YH
#define YH(i, j) ARRAY2D(yh, nyh, i, j)
#endif
// DLS001
    odepack_cpp_real &ccmax = dls1_.ccmax, &el0 = dls1_.el0, &h = dls1_.h, &hmin = dls1_.hmin, &hmxi = dls1_.hmxi, &hu = dls1_.hu, &rc = dls1_.rc, &tn = dls1_.tn, &uround = dls1_.uround;
    int &icf = dls1_.icf, &ierpj = dls1_.ierpj, &iersl = dls1_.iersl, &jcur = dls1_.jcur, &jstart = dls1_.jstart, &kflag = dls1_.kflag, &l = dls1_.l, 
        &lyh = dls1_.lyh, &lewt = dls1_.lewt, &lacor = dls1_.lacor, &lsavf = dls1_.lsavf, &lwm = dls1_.lwm, &liwm = dls1_.liwm, &meth = dls1_.meth, &miter = dls1_.miter, 
        &maxord = dls1_.maxord, &maxcor = dls1_.maxcor, &msbp = dls1_.msbp, &mxncf = dls1_.mxncf, &n = dls1_.n, &nq = dls1_.nq, &nst = dls1_.nst, &nfe = dls1_.nfe, &nje = dls1_.nje, &nqu = dls1_.nqu;
// DLSS01
    int &iplost = dlss_.iplost, &iesp = dlss_.iesp, &istatc = dlss_.istatc, &iys = dlss_.iys, &iba = dlss_.iba, &ibian = dlss_.ibian, &ibjan = dlss_.ibjan, &ibjgp = dlss_.ibjgp,
        &ipian = dlss_.ipian, &ipjan = dlss_.ipjan, &ipjgp = dlss_.ipjgp, &ipigp = dlss_.ipigp, &ipr = dlss_.ipr, &ipc = dlss_.ipc, &ipic = dlss_.ipic, &ipisp = dlss_.ipisp, &iprsp = dlss_.iprsp, &ipa = dlss_.ipa,
        &lenyh = dlss_.lenyh, &lenyhm = dlss_.lenyhm, &lenwk = dlss_.lenwk, &lreq = dlss_.lreq, &lrat = dlss_.lrat, &lrest = dlss_.lrest, &lwmin = dlss_.lwmin, &moss = dlss_.moss, &msbj = dlss_.msbj,
        &nslj = dlss_.nslj, &ngp = dlss_.ngp, &nlu = dlss_.nlu, &nnz = dlss_.nnz, &nsp = dlss_.nsp, &nzl = dlss_.nzl, &nzu = dlss_.nzu;
//
    int i, imul, ires, j, jj, jmax, jmin, k, kmax, kmin, ng;
    odepack_cpp_real con, fac, hl0, r, srur;
//
    hl0 = h * el0;
    con = -hl0;
    jcur = 1;
    nje++;
    if (miter == 1) {
        goto LABEL_100;
    } else if (miter == 2) {
        goto LABEL_200;
    }
//
// If MITER = 1, call RES, then call JAC and ADDA for each column. ------
LABEL_100:
    ires = 1;
    (*res)(neq, tn, y, s, savr, ires, user_data);
    nfe++;
    if (ires > 1) goto LABEL_600;
    kmin = ARRAY1D(iwk, ipian);
    for (j = 1; j <= n; ++j) {
        kmax = ARRAY1D(iwk, ipian + j) - 1;
        for (i = 1; i <= n; ++i) {
            ARRAY1D(rtem, i) = 0.0;
        }
        (*jac)(neq, tn, y, s, j, &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan), rtem, user_data);
        for (i = 1; i <= n; ++i) {
            ARRAY1D(rtem, i) *= con;
        }
        (*adda)(neq, tn, y, j, &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan), rtem, user_data);
        for (k = kmin; k <= kmax; ++k) {
            i = ARRAY1D(iwk, ibjan + k);
            ARRAY1D(wk, iba + k) = ARRAY1D(rtem, i);
        }
        kmin = kmax + 1;
    }
    goto LABEL_290;
//
// If MITER = 2, make NGP + 1 calls to RES to approximate J and P. ------
LABEL_200:
    ires = -1;
    (*res)(neq, tn, y, s, savr, ires, user_data);
    nfe++;
    if (ires > 1) goto LABEL_600;
    srur = ARRAY1D(wk, 1);
    jmin = ARRAY1D(iwk, ipigp);
    for (ng = 1; ng <= ngp; ++ng) {
        jmax = ARRAY1D(iwk, ipigp + ng) - 1;
        for (j = jmin; j <= jmax; ++j) {
            jj = ARRAY1D(iwk, ibjgp + j);
            r = std::max(srur * std::abs(ARRAY1D(y, jj)), 0.01 / ARRAY1D(ewt, jj));
            ARRAY1D(y, jj) += r;
        }
        (*res)(neq, tn, y, s, rtem, ires, user_data);
        nfe++;
        if (ires > 1) goto LABEL_600;
        for (j = jmin; j <= jmax; ++j) {
            jj = ARRAY1D(iwk, ibjgp + j);
            ARRAY1D(y, jj) = YH(jj, 1);
            r = std::max(srur * std::abs(ARRAY1D(y, jj)), 0.01 / ARRAY1D(ewt, jj));
            fac = -hl0 / r;
            kmin = ARRAY1D(iwk, ibian + jj);
            kmax = ARRAY1D(iwk, ibian + jj + 1) - 1;
            for (k = kmin; k <= kmax; ++k) {
                i = ARRAY1D(iwk, ibjan + k);
                ARRAY1D(rtem, i) = (ARRAY1D(rtem, i) - ARRAY1D(savr, i)) * fac;
            }
            (*adda)(neq, tn, y, jj, &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan), rtem, user_data);
            for (k = kmin; k <= kmax; ++k) {
                i = ARRAY1D(iwk, ibjan + k);
                ARRAY1D(wk, iba + k) = ARRAY1D(rtem, i);
            }
        }
        jmin = jmax + 1;
    }
    ires = 1;
    (*res)(neq, tn, y, s, savr, ires, user_data);
    nfe++;
    if (ires > 1) goto LABEL_600;
//
// Do numerical factorization of P matrix. ------------------------------
LABEL_290:
    nlu++;
    ierpj = 0;
    for (i = 1; i <= n; ++i) {
        ARRAY1D(rtem, i) = 0.0;
    }
    CDRV(n, &ARRAY1D(iwk, ipr), &ARRAY1D(iwk, ipc), &ARRAY1D(iwk, ipic), &ARRAY1D(iwk, ipian), &ARRAY1D(iwk, ipjan),
        &ARRAY1D(wk, ipa), rtem, rtem, nsp, &ARRAY1D(iwk, ipisp), &ARRAY1D(wk, iprsp), iesp, 2, iys);
    if (iys == 0) return;
    imul = (iys - 1) / n;
    ierpj = -2;
    if (imul == 8) ierpj = 1;
    if (imul == 10) ierpj = -1;
    return;
// Error return for IRES = 2 or IRES = 3 return from RES. ---------------
LABEL_600:
    ierpj = ires;
    return;
//
#ifdef YH
#undef YH
#endif
}

} // end namespace odepack_cpp
