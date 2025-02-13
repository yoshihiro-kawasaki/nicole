/**
 * @note
 * @date 2025/01/21 kawasaki
 * @date 2025/02/07 kawasaki
*/

#include "odepack.hpp"

namespace odepack_cpp
{

/**
 * @fn DUMACH
 * @brief Compute the unit roundoff of the machine.
 * @details
 *  The unit roundoff is defined as the smallest positive machine
 *  number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
 *  in a machine-independent manner.
 * @return the unit roundoff of the machine.
 * @date 2025/02/07 kawasaki
 */
double Odepack::DUMACH()
{
    double u, comp;
    u = 1.0;
    while (1) {
        u *= 0.5;
        comp = 1.0 + u;
        if (comp == 1.0) {
            break;
        }
    }
    return 2.0*u;
}

/**
 * @fn DUMSUM
 * @brief Routine to force normal storing of A + B, for DUMACH.
 * @date 2025/02/07 kawasaki
 */
void Odepack::DUMSUM(const double a, const double b, double &c)
{
    c = a + b;
    return;
}

/**
 * @fn DCFODE
 * @brief Set ODE integrator coefficients.
 * @date 2025/02/07 kawasaki
 */
void Odepack::DCFODE(const int meth, double *elco, double *tesco)
{
#ifndef ELCO
#define ELCO(i, j) MATF(elco, 13, i, j)
#endif
#ifndef TESCO
#define TESCO(i, j) MATF(tesco, 3, i, j)
#endif
//
    int i, ib, nq, nqm1, nqp1;
    double agamq, fnq, fnqm1, pc[12], pint, ragq, rqfac, rq1fac, tsign, xpin;
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
    TESCO(1,  1) = 0.0;
    TESCO(2,  1) = 2.0;
    TESCO(1,  2) = 1.0;
    TESCO(3, 12) = 0.0;
    ARRAYF(pc, 1) = 1.0;
    rqfac = 1.0;
    for (nq = 2; nq <= 12; ++nq) {
//-----------------------------------------------------------------------
// The PC array will contain the coefficients of the polynomial
// p(x) = (x+1)*(x+2)*...*(x+nq-1).
// Initially, p(x) = 1.
//-----------------------------------------------------------------------
        rq1fac = rqfac;
        rqfac  = rqfac / static_cast<double>(nq);
        nqm1   = nq - 1;
        fnqm1  = static_cast<double>(nqm1);
        nqp1   = nq + 1;
// Form coefficients of p(x)*(x+nq-1). ----------------------------------
        ARRAYF(pc, nq) = 0.0;
        for (ib = 1; ib <= nqm1; ++ib) {
            i = nqp1 - ib;
            ARRAYF(pc, i) = ARRAYF(pc, i-1) + fnqm1 * ARRAYF(pc, i);
        }
        ARRAYF(pc, 1) = fnqm1 * ARRAYF(pc, 1);
// Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        pint = ARRAYF(pc, 1);
        xpin = ARRAYF(pc, 1) / 2.0;
        tsign = 1.0;
        for (i = 2; i <= nq; ++i) {
            tsign = -tsign;
            pint = pint + tsign * ARRAYF(pc, i) / static_cast<double>(i);
            xpin = xpin + tsign * ARRAYF(pc, i) / static_cast<double>(i + 1);
        }
// Store coefficients in elco and tesco. --------------------------------
        ELCO(1, nq) = pint * rq1fac;
        ELCO(2, nq) = 1.0;
        for (i = 2; i <= nq; ++i) {
            ELCO(i + 1, nq) = rq1fac * ARRAYF(pc, i) / static_cast<double>(i);
        }
        agamq = rqfac * xpin;
        ragq = 1.0 / agamq;
        TESCO(2, nq) = ragq;
        if (nq < 12) TESCO(1, nqp1) = ragq * rqfac / static_cast<double>(nqp1);
        TESCO(3, nqm1) = ragq;
    }
    return;
//
LABEL_200:
    ARRAYF(pc, 1) = 1.0;
    rq1fac = 1.0;
    for (nq = 1; nq <= 5; ++nq) {
//-----------------------------------------------------------------------
// The PC array will contain the coefficients of the polynomial
//     p(x) = (x+1)*(x+2)*...*(x+nq).
// Initially, p(x) = 1.
//-----------------------------------------------------------------------
        fnq = static_cast<double>(nq);
        nqp1 = nq + 1;
// Form coefficients of p(x)*(x+nq). ------------------------------------
        ARRAYF(pc, nqp1) = 0.0;
        for (ib = 1; ib <= nq; ++ib) {
            i = nq + 2 - ib;
            ARRAYF(pc, i) = ARRAYF(pc, i - 1) + fnq * ARRAYF(pc, i);
        }
        ARRAYF(pc, 1) = fnq * ARRAYF(pc, 1);
// Store coefficients in elco and tesco. --------------------------------
        for (i = 1; i <= nqp1; ++i) {
            ELCO(i, nq) = ARRAYF(pc, i) / ARRAYF(pc, 2);
        }
        ELCO(2, nq) = 1.0;
        TESCO(1, nq) = rq1fac;
        TESCO(2, nq) = (static_cast<double>(nqp1)) / ELCO(1, nq);
        TESCO(3, nq) = (static_cast<double>(nq + 2)) / ELCO(1, nq);
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
 * @fn DINTDY
 * @brief Interpolate solution derivatives.
 * @date 2025/02/07 kawasaki
 */
void Odepack::DINTDY(const double t, const int k,  const double *yh, const int nyh,  double *dky, int &iflag)
{
#ifndef YH
#define YH(i, j) MATF(yh, nyh, i, j)
#endif
//
    int i, ic, j, jj, jb, jb2, jj1, jp1;
    double c, r, s, tp;
    std::string msg;
//
//***FIRST EXECUTABLE STATEMENT  DINTDY
    iflag = 0;
    if (k < 0 || k > dls1_.nq) goto LABEL_80;
    tp = dls1_.tn - dls1_.hu - 100.0 * dls1_.uround * (std::abs(dls1_.tn) + std::abs(dls1_.hu)) * SIGN(dls1_.hu);
    if ((t - tp) * (t - dls1_.tn) > 0.0) goto LABEL_90;
//
    s = (t - dls1_.tn) / dls1_.h;
    ic = 1;
    if (k == 0) goto LABEL_15;
    jj1 = dls1_.l - k;
    for (jj = jj1; jj <= dls1_.nq; ++jj) {
        ic = ic * jj;
    }
LABEL_15:
    c = static_cast<double>(ic);
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(dky, i) = c * YH(i, dls1_.l);
    }
    if (k == dls1_.nq) goto LABEL_55;
    jb2 = dls1_.nq - k;
    for (jb = 1; jb <= jb2; ++jb) {
        j = dls1_.nq - jb;
        jp1 = j + 1;
        ic = 1;
        if (k == 0) goto LABEL_35;
        jj1 = jp1 - k;
        for (jj = jj1; jj <= j; ++jj) {
            ic = ic * jj;
        }
LABEL_35:
        c = static_cast<double>(ic);
        for (i = 1; i <= dls1_.n; ++i) {
            ARRAYF(dky, i) = c * YH(i, jp1) + s * ARRAYF(dky, i);
        }
    }
    if (k == 0) return;
LABEL_55:
    r = std::pow(dls1_.h, -static_cast<double>(k));
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(dky, i) = r * ARRAYF(dky, i);
    }
    return;
//
LABEL_80:
    msg = "DLSODE_CPP::DINTDY-  K (=I1) illegal      ";
    XERRWD(msg, 30, 51, 0, 1, k, 0, 0, 0.0, 0.0);
    iflag = -1;
    return;
LABEL_90:
    msg = "DLSODE_CPP::DINTDY-  T (=R1) illegal      ";
    XERRWD(msg, 30, 51, 0, 1, k, 0, 0, 0.0, 0.0);
    iflag = -2;
    return;
//
#ifdef YH
#undef YH
#endif
}

/**
 * @fn DPREPJ
 * @brief Compute and process Newton iteration matrix.
 * @date 2025/02/07 kawasaki
 */
void Odepack::DPREPJ(const int neq,  double *y, double *yh, const int nyh, const double *ewt, 
    double *ftem, double *savf, double *wm, int *iwm, ODEPACK_FUNCTION f, 
    ODEPACK_JACOBIAN1 jac, void *user_data)
{
#ifndef YH
#define YH(i, j) MATF(yh, nyh, i, j)
#endif
//
    int i, i1, i2, ier, ii, j, j1, jj, lenp, mba, mband, meb1, meband, ml, ml3, mu, np1;
    double con, di, fac, hl0, r, r0, srur, yi, yj, yjj;
//***FIRST EXECUTABLE STATEMENT  DPREPJ
    dls1_.nje++;
    dls1_.ierpj = 0;
    dls1_.jcur  = 1;
    hl0 = dls1_.h * dls1_.el0;
    if (dls1_.miter == 1) {
        goto LABEL_100;
    } else if (dls1_.miter == 2) {
        goto LABEL_200;
    } else if (dls1_.miter == 3) {
        goto LABEL_300;
    } else if (dls1_.miter == 4) {
        goto LABEL_400;
    } else if (dls1_.miter == 5) {
        goto LABEL_500;
    }
// If MITER = 1, call JAC and multiply by scalar. -----------------------
LABEL_100:
    lenp = dls1_.n * dls1_.n;
    for (i = 1; i <= lenp; ++i) {
        ARRAYF(wm, i + 2) = 0.0;
    }
    (*jac)(neq, dls1_.tn, y, 0, 0, &ARRAYF(wm, 3), dls1_.n, user_data);
    con = -hl0;
    for (i = 1; i <= lenp; ++i) {
        ARRAYF(wm, i + 2) *= con;
    }
    goto LABEL_240;
// If MITER = 2, make N calls to F to approximate J. --------------------
LABEL_200:
    fac = DVNORM(dls1_.n, savf, ewt);
    r0 = 1000.0 * std::abs(dls1_.h) * dls1_.uround * static_cast<double>(dls1_.n) * fac;
    if (r0 == 0.0) r0 = 1.0;
    srur = ARRAYF(wm, 1);
    j1 = 2;
    for (j = 1; j <= dls1_.n; ++j) {
        yj = ARRAYF(y, j);
        r = std::max(srur * std::abs(yj), r0 / ARRAYF(ewt, j));
        ARRAYF(y, j) += r;
        fac = -hl0 / r;
        (*f)(neq, dls1_.tn, y, ftem, user_data);
        for (i = 1; i <= dls1_.n; ++i) {
            ARRAYF(wm, i + j1) = (ARRAYF(ftem, i) - ARRAYF(savf, i)) * fac;
        }
        ARRAYF(y, j) = yj;
        j1 += dls1_.n;
    }
    dls1_.nfe += dls1_.n;
// Add identity matrix. -------------------------------------------------
LABEL_240:
    j = 3;
    np1 = dls1_.n + 1;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(wm, j) += 1.0;
        j += np1;
    }
// Do LU decomposition on P. --------------------------------------------
    DGEFA(&ARRAYF(wm, 3), dls1_.n, dls1_.n, &ARRAYF(iwm, 21), ier);
    if (ier != 0) dls1_.ierpj = 1;
    return;
// If MITER = 3, construct a diagonal approximation to J and P. ---------
LABEL_300:
    ARRAYF(wm, 2) = hl0;
    r = dls1_.el0 * 0.1;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = ARRAYF(y, i) + r * (dls1_.h * ARRAYF(savf, i) - YH(i, 2));
    }
    (*f)(neq, dls1_.tn, y, &ARRAYF(wm, 3), user_data);
    dls1_.nfe++;
    for (i = 1; i <= dls1_.n; ++i) {
        r0 = dls1_.h * ARRAYF(savf, i) - YH(i, 2);
        di = 0.1 * r0 - dls1_.h * (ARRAYF(wm, i + 2) - ARRAYF(savf, i));
        ARRAYF(wm, i + 2) = 1.0;
        if (std::abs(r0) < dls1_.uround / ARRAYF(ewt, i)) continue;
        if (std::abs(di) == 0.0) goto LABEL_330;
        ARRAYF(wm, i + 2) = 0.1 * r0 / di;
    }
    return;
LABEL_330:
    dls1_.ierpj = 1;
    return;
// If MITER = 4, call JAC and multiply by scalar. -----------------------
LABEL_400:
    ml = ARRAYF(iwm, 1);
    mu = ARRAYF(iwm, 2);
    ml3 = ml + 3;
    mband = ml + mu + 1;
    meband = mband + ml;
    lenp = meband * dls1_.n;
    for (i = 1; i <= lenp; ++i) {
        ARRAYF(wm, i + 2) = 0.0;
    }
    (*jac)(neq, dls1_.tn, y, ml, mu, &ARRAYF(wm, ml3), meband, user_data);
    con = -hl0;
    for (i = 1; i <= lenp; ++i) {
        ARRAYF(wm, i + 2) *= con;
    }
    goto LABEL_570;
// If MITER = 5, make MBAND calls to F to approximate J. ----------------
LABEL_500:
    ml = ARRAYF(iwm, 1);
    mu = ARRAYF(iwm, 2);
    mband = ml + mu + 1;
    mba = std::min(mband, dls1_.n);
    meband = mband + ml;
    meb1 = meband - 1;
    srur = ARRAYF(wm, 1);
    fac = DVNORM(dls1_.n, savf, ewt);
    r0 = 1000.0 * std::abs(dls1_.h) * dls1_.uround * static_cast<double>(dls1_.n) * fac;
    if (r0 == 0.0) r0 = 1.0;
    for (j = 1; j <= mba; ++j) {
        for (i = j; i <= dls1_.n; i += mband) {
            yi = ARRAYF(y, i);
            r = std::max(srur * std::abs(yi), r0 / ARRAYF(ewt, i));
            ARRAYF(y, i) += r;
        }
        (*f)(neq, dls1_.tn, y, ftem, user_data);
        for (jj = j; jj <= dls1_.n; jj += mband) {
            ARRAYF(y, jj) = YH(jj, 1);
            yjj = ARRAYF(y, jj);
            r = std::max(srur * std::abs(yjj), r0 / ARRAYF(ewt, jj));
            fac = -hl0 / r;
            i1 = std::max(jj - mu, 1);
            i2 = std::min(jj + ml, dls1_.n);
            ii = jj * meb1 - ml + 2;
            for (i = i1; i <= i2; ++i) {
                ARRAYF(wm, ii + i) = (ARRAYF(ftem, i) - ARRAYF(savf, i)) * fac;
            }
        }
    }
    dls1_.nfe += mba;
// Add identity matrix. -------------------------------------------------
LABEL_570:
    ii = mband + 2;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(wm, ii) += 1.0;
        ii += meband;
    }
// Do LU decomposition of P. --------------------------------------------
    DGBFA(&ARRAYF(wm, 3), meband, dls1_.n, ml, mu, &ARRAYF(iwm, 21), ier);
    if (ier != 0) dls1_.ierpj = 1;
    return;
//
#ifdef YH
#undef YH
#endif
}

/**
 * @fn DSOLSY
 * @brief ODEPACK linear system solver.
 * @date 2025/02/07 kawasaki
 */
void Odepack::DSOLSY(double *wm, int *iwm, double *x, double *tem)
{
    int i, meband, ml, mu;
    double di, hl0, phl0, r;
//***FIRST EXECUTABLE STATEMENT  DSOLSY
    dls1_.iersl = 0;
    if (dls1_.miter == 1 || dls1_.miter == 2) {
        goto LABEL_100;
    } else if (dls1_.miter == 3) {
        goto LABEL_300;
    } else if (dls1_.miter == 4 || dls1_.miter == 5) {
        goto LABEL_400;
    }
LABEL_100:
    DGESL(&ARRAYF(wm, 3), dls1_.n, dls1_.n, &ARRAYF(iwm, 21), x, 0);
    return;
//
LABEL_300:
    phl0 = ARRAYF(wm, 2);
    hl0 = dls1_.h * dls1_.el0;
    ARRAYF(wm, 2) = hl0;
    if (hl0 == phl0) goto LABEL_330;
    r = hl0 / phl0;
    for (i = 1; i <= dls1_.n; ++i) {
        di = 1.0 - r * (1.0 - 1.0 / ARRAYF(wm, i + 2));
        if (std::abs(di) == 0.0) goto LABEL_390;
        ARRAYF(wm, i + 2) = 1.0 / di;
    }
LABEL_330:
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(x, i) *= ARRAYF(wm, i + 2);
    }
    return;
LABEL_390:
    dls1_.iersl = 1;
    return;
//
LABEL_400:
    ml = ARRAYF(iwm, 1);
    mu = ARRAYF(iwm, 2);
    meband = 2 * ml + mu + 1;
    DGBSL(&ARRAYF(wm, 3), meband, dls1_.n, ml, mu, &ARRAYF(iwm, 21), x, 0);
    return;
}

/**
 * @fn DSRCOM
 * @brief Save/restore ODEPACK COMMON blocks.
 * @date 2025/02/07 kawasaki
 */
void Odepack::DSRCOM(double *rsav, int *isav, const int job)
{
    int i;
    const int lenrls = 218;
    const int lenils = 37; 
//***FIRST EXECUTABLE STATEMENT  DSRCOM
    if (job == 2) goto LABEL_100;
//
    for (i = 1; i <= lenrls; ++i) {
        ARRAYF(rsav, i) = ARRAYF(dls1_.rls, i);
    }
    for (i = 1; i <= lenils; ++i) {
        ARRAYF(isav, i) = ARRAYF(dls1_.ils, i);
    }
    return;
//
LABEL_100:
    for (i = 1; i <= lenrls; ++i) {
        ARRAYF(dls1_.rls, i) = ARRAYF(rsav, i);
    }
    for (i = 1; i <= lenils; ++i) {
        ARRAYF(dls1_.ils, i) = ARRAYF(isav, i);
    }
    return;
}

/**
 * @fn DSTODE
 * @brief Performs one step of an ODEPACK integration.
 * @date 2025/02/07 kawasaki
 */
template <typename ODEPACK_JACOBIAN>
void Odepack::DSTODE(const int neq, double *y, double *yh, int nyh, double*yh1, 
    double *ewt, double *savf, double *acor, double *wm, void *iwm_in, 
    ODEPACK_FUNCTION f, ODEPACK_JACOBIAN jac, FUNC_PJAC<ODEPACK_JACOBIAN> pjac, FUNC_SLVS slvs, 
    void *user_data)
{
#ifndef YH
#define YH(i, j) MATF(yh, nyh, i, j)
#endif
#ifndef ELCO
#define ELCO(i, j) MATF(dls1_.elco, 13, i, j)
#endif
#ifndef TESCO
#define TESCO(i, j) MATF(dls1_.tesco, 3, i, j)
#endif
//
    int i, i1, j, jb, iredo, iret, m, ncf, newq;
    double dcon, ddn, del, delp, dsm, dup, exdn, exup, r, rh, rhdn, rhup, told, exsm, rhsm;
//
    int *iwm = static_cast<int*>(iwm_in);
//***FIRST EXECUTABLE STATEMENT  DSTODE
    dls1_.kflag = 0;
    told = dls1_.tn;
    ncf = 0;
    dls1_.ierpj = 0;
    dls1_.iersl = 0;
    dls1_.jcur = 0;
    dls1_.icf = 0;
    delp = 0.0;
    if (dls1_.jstart > 0)   goto LABEL_200;
    if (dls1_.jstart == -1) goto LABEL_100;
    if (dls1_.jstart == -2) goto LABEL_160;
//-----------------------------------------------------------------------
// On the first call, the order is set to 1, and other variables are
// initialized.  RMAX is the maximum ratio by which H can be increased
// in a single step.  It is initially 1.E4 to compensate for the small
// initial H, but then is normally equal to 10.  If a failure
// occurs (in corrector convergence or error test), RMAX is set to 2
// for the next increase.
//-----------------------------------------------------------------------
    dls1_.lmax = dls1_.maxord + 1;
    dls1_.nq = 1;
    dls1_.l = 2;
    dls1_.ialth = 2;
    dls1_.rmax = 10000.0;
    dls1_.rc = 0.0;
    dls1_.el0 = 1.0;
    dls1_.crate = 0.7;
    dls1_.hold = dls1_.h;
    dls1_.meo = dls1_.meth;
    dls1_.nslp = 0;
    dls1_.ipup = dls1_.miter;
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
    dls1_.ipup = dls1_.miter;
    dls1_.lmax = dls1_.maxord + 1;
    if (dls1_.ialth == 1) dls1_.ialth = 2;
    if (dls1_.meth == dls1_.meo) goto LABEL_110;
    DCFODE(dls1_.meth, dls1_.elco, dls1_.tesco);
    dls1_.meo = dls1_.meth;
    if (dls1_.nq > dls1_.maxord) goto LABEL_120;
    dls1_.ialth = dls1_.l;
    iret = 1;
    goto LABEL_150;
LABEL_110:
    if (dls1_.nq <= dls1_.maxord) goto LABEL_160;
LABEL_120:
    dls1_.nq = dls1_.maxord;
    dls1_.l = dls1_.lmax;
    for (i = 1; i <= dls1_.l; ++i) {
        ARRAYF(dls1_.el, i) =  ELCO(i, dls1_.nq);
    }
    dls1_.nqnyh = dls1_.nq * nyh;
    dls1_.rc = dls1_.rc * ARRAYF(dls1_.el, 1) / dls1_.el0;
    dls1_.el0 = ARRAYF(dls1_.el, 1);
    dls1_.coint = 0.5 / static_cast<double>(dls1_.nq + 2);
    ddn = DVNORM(dls1_.n, savf, ewt) / TESCO(1, dls1_.l);
    exdn = 1.0 / static_cast<double>(dls1_.l);
    rhdn = 1.0 / (1.3 * std::pow(ddn, exdn) + 0.0000013);
    rh = std::min(rhdn, 1.0);
    iredo = 3;
    if (dls1_.h == dls1_.hold) goto LABEL_170;
    rh = std::min(rh, std::abs(dls1_.h / dls1_.hold));
    dls1_.h = dls1_.hold;
    goto LABEL_175;
//-----------------------------------------------------------------------
// DCFODE is called to get all the integration coefficients for the
// current METH.  Then the EL vector and related constants are reset
// whenever the order NQ is changed, or at the start of the problem.
//-----------------------------------------------------------------------
LABEL_140:
    DCFODE(dls1_.meth, dls1_.elco, dls1_.tesco);
LABEL_150:
    for (i = 1; i <= dls1_.l; ++i) {
        ARRAYF(dls1_.el, i) = ELCO(i, dls1_.nq);
    }
    dls1_.nqnyh = dls1_.nq * nyh;
    dls1_.rc = dls1_.rc * ARRAYF(dls1_.el, 1) / dls1_.el0;
    dls1_.el0 = ARRAYF(dls1_.el, 1);
    dls1_.coint = 0.5 / static_cast<double>(dls1_.nq + 2);
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
    if (dls1_.h == dls1_.hold) goto LABEL_200;
    rh = dls1_.h / dls1_.hold;
    dls1_.h = dls1_.hold;
    iredo = 3;
    goto LABEL_175;
LABEL_170:
    rh = std::max(rh, dls1_.hmin / std::abs(dls1_.h));
LABEL_175:
    rh = std::min(rh, dls1_.rmax);
    rh = rh / std::max(1.0, std::abs(dls1_.h) * dls1_.hmxi * rh);
    r = 1.0;
    for (j = 2; j <= dls1_.l; ++j) {
        r = r * rh;
        for (i = 1; i <= dls1_.n; ++i) {
            YH(i, j) *= r;
        }
    }
    dls1_.h *= rh;
    dls1_.rc *= rh;
    dls1_.ialth = dls1_.l;
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
    if (std::abs(dls1_.rc - 1.0) > dls1_.ccmax) dls1_.ipup = dls1_.miter; 
    if (dls1_.nst >= (dls1_.nslp + dls1_.msbp)) dls1_.ipup = dls1_.miter;
    dls1_.tn += dls1_.h;
    i1 = dls1_.nqnyh + 1;
    for (jb = 1; jb <= dls1_.nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= dls1_.nqnyh; ++i) {
            ARRAYF(yh1, i) += ARRAYF(yh1, i + dls1_.nyh);
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
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = YH(i, 1);
    }
    (*f)(neq, dls1_.tn, y, savf, user_data);
    dls1_.nfe++;
    if (dls1_.ipup <= 0) goto LABEL_250;
//-----------------------------------------------------------------------
// If indicated, the matrix P = I - h*el(1)*J is reevaluated and
// preprocessed before starting the corrector iteration.  IPUP is set
// to 0 as an indicator that this has been done.
//-----------------------------------------------------------------------
    (this->*pjac)(neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac, user_data);
    dls1_.ipup = 0;
    dls1_.rc = 1.0;
    dls1_.nslp = dls1_.nst;
    dls1_.crate = 0.7;
    if (dls1_.ierpj != 0) goto LABEL_430;
LABEL_250:
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(acor, i) = 0.0;
    }
LABEL_270:
    if (dls1_.miter != 0) goto LABEL_350;
//-----------------------------------------------------------------------
// In the case of functional iteration, update Y directly from
// the result of the last function evaluation.
//-----------------------------------------------------------------------
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(savf, i) = dls1_.h * ARRAYF(savf, i) - YH(i, 2);
        ARRAYF(y, i) = ARRAYF(savf, i) - ARRAYF(acor, i);
    }
    del = DVNORM(dls1_.n, y, ewt);
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = YH(i, 1) + ARRAYF(dls1_.el, 1) * ARRAYF(savf, i);
        ARRAYF(acor, i) = ARRAYF(savf, i);
    }
    goto LABEL_400;
//-----------------------------------------------------------------------
// In the case of the chord method, compute the corrector error,
// and solve the linear system with that as right-hand side and
// P as coefficient matrix.
//-----------------------------------------------------------------------
LABEL_350:
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = dls1_.h * ARRAYF(savf, i) - (YH(i, 2) + ARRAYF(acor, i));
    }
    (this->*slvs)(wm, iwm, y, savf);
    if (dls1_.iersl < 0) goto LABEL_430;
    if (dls1_.iersl > 0) goto LABEL_410;
    del = DVNORM(dls1_.n, y, ewt);
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(acor, i) += ARRAYF(y, i);
        ARRAYF(y, i) = YH(i, 1) + ARRAYF(dls1_.el, 1) * ARRAYF(acor, i);
    }
//-----------------------------------------------------------------------
// Test for convergence.  If M.gt.0, an estimate of the convergence
// rate constant is stored in CRATE, and this is used in the test.
//-----------------------------------------------------------------------
LABEL_400:
    if (m != 0) dls1_.crate = std::max(0.2 * dls1_.crate, del / delp);
    dcon = del * std::min(1.0, 1.5 * dls1_.crate) / (TESCO(2, dls1_.nq) * dls1_.coint);
    if (dcon <= 1.0) goto LABEL_450;
    m++;
    if (m == dls1_.maxcor) goto LABEL_410;
    if (m >= 2 && del > 2.0 * delp) goto LABEL_410;
    delp = del;
    (*f)(neq, dls1_.tn, y, savf, user_data);
    dls1_.nfe++;
    goto LABEL_270;
//-----------------------------------------------------------------------
// The corrector iteration failed to converge.
// If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
// the next try.  Otherwise the YH array is retracted to its values
// before prediction, and H is reduced, if possible.  If H cannot be
// reduced or MXNCF failures have occurred, exit with KFLAG = -2.
//-----------------------------------------------------------------------
LABEL_410:
    if (dls1_.miter == 0 || dls1_.jcur == 1) goto LABEL_430;
    dls1_.icf = 1;
    dls1_.ipup = dls1_.miter;
    goto LABEL_220;
LABEL_430:
    dls1_.icf = 2;
    ncf++;
    dls1_.rmax = 2.0;
    dls1_.tn = told;
    i1 = dls1_.nqnyh + 1;
    for (jb = 1; jb <= dls1_.nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= dls1_.nqnyh; ++i) {
            ARRAYF(yh1, i) -= ARRAYF(yh1, i + nyh);
        }
    }
    if (dls1_.ierpj < 0 || dls1_.iersl < 0) goto LABEL_680;
    if (std::abs(dls1_.h) <= dls1_.hmin * 1.00001) goto LABEL_670;
    if (ncf == dls1_.mxncf) goto LABEL_670;
    rh = 0.25;
    dls1_.ipup = dls1_.miter;
    iredo = 1;
    goto LABEL_170;
//-----------------------------------------------------------------------
// The corrector has converged.  JCUR is set to 0
// to signal that the Jacobian involved may need updating later.
// The local error test is made and control passes to statement 500
// if it fails.
//-----------------------------------------------------------------------
LABEL_450:
    dls1_.jcur = 0;
    if (m == 0) dsm = del / TESCO(2, dls1_.nq);
    if (m > 0) dsm = DVNORM(dls1_.n, acor, ewt) / TESCO(2, dls1_.nq);
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
    dls1_.kflag = 0;
    iredo = 0;
    dls1_.nst++;
    dls1_.hu = dls1_.h;
    dls1_.nqu = dls1_.nq;
    for (j = 1; j <= dls1_.l; ++j) {
        for (i = 1; i <= dls1_.n; ++i) {
            YH(i, j) = YH(i, j) + ARRAYF(dls1_.el, j) * ARRAYF(acor, i);
        }
    }
    dls1_.ialth--;
    if (dls1_.ialth == 0) goto LABEL_520;
    if (dls1_.ialth > 1) goto LABEL_700;
    if (dls1_.l == dls1_.lmax) goto LABEL_700;
    for (i = 1; i <= dls1_.n; ++i) {
        YH(i, dls1_.lmax) = ARRAYF(acor, i);
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
    dls1_.kflag--;
    dls1_.tn = told;
    i1 = dls1_.nqnyh + 1;
    for (jb = 1; jb <= dls1_.nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= dls1_.nqnyh; ++i) {
            ARRAYF(yh1, i) -= ARRAYF(yh1, i + nyh);
        }
    }
    dls1_.rmax = 2.0;
    if (std::abs(dls1_.h) <= dls1_.hmin * 1.00001) goto LABEL_660;
    if (dls1_.kflag <= -3) goto LABEL_640;
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
    if (dls1_.l == dls1_.lmax) goto LABEL_540;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(savf, i) = ARRAYF(acor, i) - YH(i, dls1_.lmax);
    }
    dup = DVNORM(dls1_.n, savf, ewt) / TESCO(3, dls1_.nq);
    exup = 1.0 / static_cast<double>(dls1_.l + 1);
    rhup = 1.0 / (1.4 * std::pow(dup, exup) + 0.0000014);
LABEL_540:
    exsm = 1.0 / static_cast<double>(dls1_.l);
    rhsm = 1.0 / (1.2 * std::pow(dsm, exsm) + 0.0000012);
    rhdn = 0.0;
    if (dls1_.nq == 1) goto LABEL_560;
    ddn = DVNORM(dls1_.n, &YH(1, dls1_.l), ewt) / TESCO(1, dls1_.nq);
    exdn = 1.0 / static_cast<double>(dls1_.nq);
    rhdn = 1.0 / (1.3 * std::pow(ddn, exdn) + 0.0000013);
LABEL_560:
    if (rhsm >= rhup) goto LABEL_570;
    if (rhup > rhdn) goto LABEL_590;
    goto LABEL_580;
LABEL_570:
    if (rhsm < rhdn) goto LABEL_580;
    newq = dls1_.nq;
    rh = rhsm;
    goto LABEL_620;
LABEL_580:
    newq = dls1_.nq - 1;
    rh = rhdn;
    if (dls1_.kflag < 0 && rh > 1.0) rh = 1.0;
    goto LABEL_620;
LABEL_590:
    newq = dls1_.l;
    rh = rhup;
    if (rh < 1.1) goto LABEL_610;
    r = ARRAYF(dls1_.el, dls1_.l) / static_cast<double>(dls1_.l);
    for (i = 1; i <= dls1_.n; ++i) {
        YH(i, newq + 1) = ARRAYF(acor, i) * r;
    }
    goto LABEL_630;
LABEL_610:
    dls1_.ialth = 3;
    goto LABEL_700;
LABEL_620:
    if (dls1_.kflag == 0 && rh < 1.1) goto LABEL_610;
    if (dls1_.kflag <= -2) rh = std::min(rh, 0.2);
//-----------------------------------------------------------------------
// If there is a change of order, reset NQ, l, and the coefficients.
// In any case H is reset according to RH and the YH array is rescaled.
// Then exit from 690 if the step was OK, or redo the step otherwise.
//-----------------------------------------------------------------------
    if (newq == dls1_.nq) goto LABEL_170;
LABEL_630:
    dls1_.nq = newq;
    dls1_.l = dls1_.nq + 1;
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
    if (dls1_.kflag == -10) goto LABEL_660;
    rh = 0.1;
    rh = std::max(dls1_.hmin / std::abs(dls1_.h), rh);
    dls1_.h *= rh;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = YH(i, 1);
    }
    (*f)(neq, dls1_.tn, y, savf, user_data);
    dls1_.nfe++;
    for (i = 1; i <= dls1_.n; ++i) {
        YH(i, 2) = dls1_.h * ARRAYF(savf, i);
    }
    dls1_.ipup = dls1_.miter;
    dls1_.ialth = 5;
    if (dls1_.nq == 1) goto LABEL_200;
    dls1_.nq = 1;
    dls1_.l = 2;
    iret = 3;
    goto LABEL_150;
//-----------------------------------------------------------------------
// All returns are made through this section.  H is saved in HOLD
// to allow the caller to change H on the next step.
//-----------------------------------------------------------------------
LABEL_660:
    dls1_.kflag = -1;
    goto LABEL_720;
LABEL_670:
    dls1_.kflag = -2;
    goto LABEL_720;
LABEL_680:
    dls1_.kflag = -3;
    goto LABEL_720;
LABEL_690:
    dls1_.rmax = 10.0;
LABEL_700:
    r = 1.0 / TESCO(2, dls1_.nqu);
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(acor, i) *= r;
    }
LABEL_720:
    dls1_.hold   = dls1_.h;
    dls1_.jstart = 1;
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
template void Odepack::DSTODE<ODEPACK_JACOBIAN1>(const int neq, double *y, double *yh, int nyh, double*yh1, 
    double *ewt, double *savf, double *acor, double *wm, void *iwm_in, 
    ODEPACK_FUNCTION f, ODEPACK_JACOBIAN1 jac, FUNC_PJAC<ODEPACK_JACOBIAN1> pjac, FUNC_SLVS slvs, 
    void *user_data);

template void Odepack::DSTODE<ODEPACK_JACOBIAN2>(const int neq, double *y, double *yh, int nyh, double*yh1, 
    double *ewt, double *savf, double *acor, double *wm, void *iwm_in, 
    ODEPACK_FUNCTION f, ODEPACK_JACOBIAN2 jac, FUNC_PJAC<ODEPACK_JACOBIAN2> pjac, FUNC_SLVS slvs, 
    void *user_data);

/**
 * @fn DEWSET
 * @brief Set error weight vector.
 * @date 2025/02/07 kawasaki
 */
void Odepack::DEWSET(const int n, const int itol, const double *rtol, const double *atol, const double *ycur, double *ewt)
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
        ARRAYF(ewt, i) = ARRAYF(rtol, 1) * std::abs(ARRAYF(ycur, i)) + ARRAYF(atol, 1);
    }
    return;  
LABEL_20:
    for (i = 1; i <= n; ++i) {   
        ARRAYF(ewt, i) = ARRAYF(rtol, 1) * std::abs(ARRAYF(ycur, i)) + ARRAYF(atol, i);
    }
    return;
LABEL_30:
    for (i = 1; i <= n; ++i) {   
        ARRAYF(ewt, i) = ARRAYF(rtol, i) * std::abs(ARRAYF(ycur, i)) + ARRAYF(atol, 1);
    }
    return;
LABEL_40:
    for (i = 1; i <= n; ++i) {   
        ARRAYF(ewt, i) = ARRAYF(rtol, i) * std::abs(ARRAYF(ycur, i)) + ARRAYF(atol, i);
    }
    return;
}

/**
 * @fn DVNORM
 * @brief Weighted root-mean-square vector norm.
 * @date 2025/02/07 kawasaki
 */
double Odepack::DVNORM(const int n, const double *v, const double *w)
{
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += v[i]*v[i] * w[i]*w[i];
    }
    return std::sqrt(sum / static_cast<double>(n));
}

/**
 * @fn DIPREP
 * @brief This routine serves as an interface between the driver and Subroutine DPREP
 * @date 2025/02/07 kawasaki
 */
void Odepack::DIPREP(const int neq, double *y, double *rwork, int *ia, int *ja,
    int &ipflag, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN2 jac, void *user_data)
{
    int i, imax, lewtn, lyhd, lyhn;
//
    ipflag = 0;
// Call DPREP to do matrix preprocessing operations. --------------------
    DPREP(neq, y, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lsavf), &ARRAYF(rwork, dls1_.lewt),
        &ARRAYF(rwork, dls1_.lacor), ia, ja, &ARRAYF(rwork, dls1_.lwm), &ARRAYF(rwork, dls1_.lwm), 
        ipflag, f, jac, user_data);
    dlss_.lenwk = std::max(dlss_.lreq, dlss_.lwmin);
    if (ipflag < 0) return;
// If DPREP was successful, move YH to end of required space for WM. ----
    lyhn = dls1_.lwm + dlss_.lenwk;
    if (lyhn > dls1_.lyh) return;
    lyhd = dls1_.lyh - lyhn;
    if (lyhd == 0) goto LABEL_20;
    imax = lyhn - 1 + dlss_.lenyhm;
    for (i = lyhn; i <= imax; ++i) {
        ARRAYF(rwork, i) = ARRAYF(rwork, i + lyhd);
    }
    dls1_.lyh = lyhn;
// Reset pointers for SAVF, EWT, and ACOR. ------------------------------
LABEL_20:
    dls1_.lsavf = dls1_.lyh + dlss_.lenyh;
    lewtn = dls1_.lsavf + dls1_.n;
    dls1_.lacor = lewtn + dls1_.n;
    if (dlss_.istatc == 3) goto LABEL_40;
// If ISTATE = 1, move EWT (left) to its new position. ------------------
    if (lewtn > dls1_.lewt) return;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(rwork, i + lewtn - 1) = ARRAYF(rwork, i + dls1_.lewt - 1);
    }
LABEL_40:
    dls1_.lewt = lewtn;
    return;
}

/**
 * @fn DPREP
 * @brief This routine performs preprocessing related to the sparse linear systems
 * @date 2025/02/07 kawasaki
 */
void Odepack::DPREP(const int neq, double *y, double *yh, double *savf, double *ewt,
        double *ftem, int *ia, int *ja, double *wk, void *iwk_in, int &ipper, 
        ODEPACK_FUNCTION f, ODEPACK_JACOBIAN2 jac, void *user_data)
{
    int i, ibr, ier, ipil, ipiu, iptt1, iptt2, j, jfound, k, 
        knew, kmax, kmin, ldif, lenigp, liwk, maxg, np1, nzsut;
    double dq, dyj, erwt, fac, yj;
//
    int *iwk = static_cast<int*>(iwk_in);
//
    dlss_.ibian = dlss_.lrat * 2;
    dlss_.ipian = dlss_.ibian + 1;
    np1         = dls1_.n + 1;
    dlss_.ipjan = dlss_.ipian + np1;
    dlss_.ibjan = dlss_.ipjan - 1;
    liwk        = dlss_.lenwk * dlss_.lrat;
    if ((dlss_.ipjan + dls1_.n - 1) > liwk) goto LABEL_210;
    if (dlss_.moss == 0) goto LABEL_30;
//
    if (dlss_.istatc == 3) goto LABEL_20;
// ISTATE = 1 and MOSS .ne. 0.  Perturb Y for structure determination. --
    for (i = 1; i <= dls1_.n; ++i) {
        erwt = 1.0 / ARRAYF(ewt, i);
        fac = 1.0 + 1.0 / (static_cast<double>(i) + 1.0);
        ARRAYF(y, i) += fac * std::abs(erwt) * SIGN(ARRAYF(y, i));
    }
    if (dlss_.moss == 1) {
        goto LABEL_70;
    } else if (dlss_.moss == 2) {
        goto LABEL_100;
    }
//
LABEL_20:
// ISTATE = 3 and MOSS .ne. 0.  Load Y from YH(*,1). --------------------
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = ARRAYF(yh, i);
    }
    if (dlss_.moss == 1) {
        goto LABEL_70;
    } else if (dlss_.moss == 2) {
        goto LABEL_100;
    }
// MOSS = 0.  Process user's IA,JA.  Add diagonal entries if necessary. -
LABEL_30:
    knew = dlss_.ipjan;
    kmin = ARRAYF(ia, 1);
    ARRAYF(iwk, dlss_.ipian) = 1;
    for (j = 1; j <= dls1_.n; ++j) {
        jfound = 0;
        kmax = ARRAYF(ia, j + 1) - 1;
        if (kmin > kmax) goto LABEL_45;
        for (k = kmin; k <= kmax; ++k) {
            i = ARRAYF(ja, k);
            if (i == j) jfound = 1;
            if (knew > liwk) goto LABEL_210;
            ARRAYF(iwk, knew) = i;
            knew++;
        }
        if (jfound == 1) goto LABEL_50;
LABEL_45:
        if (knew > liwk) goto LABEL_210;
        ARRAYF(iwk, knew) = j;
        knew++;
LABEL_50:
        ARRAYF(iwk, dlss_.ipian + j) = knew + 1 - dlss_.ipjan;
        kmin = kmax + 1;
    }
    goto LABEL_140;
//
// MOSS = 1.  Compute structure from user-supplied Jacobian routine JAC.
LABEL_70:
// A dummy call to F allows user to create temporaries for use in JAC. --
    (*f)(neq, dls1_.tn, y, savf, user_data);
    k = dlss_.ipjan;
    ARRAYF(iwk, dlss_.ipian) = 1;
    for (j = 1; j <= dls1_.n; ++j) {
        if (k > liwk) goto LABEL_210;
        ARRAYF(iwk, k) = j;
        k++;
        for (i = 1; i <= dls1_.n; ++i) {
            ARRAYF(savf, i) = 0.0;
        }
        (*jac)(neq, dls1_.tn, y, j, &ARRAYF(iwk, dlss_.ipian), &ARRAYF(iwk, dlss_.ipjan), savf, user_data);
        for (i = 1; i <= dls1_.n; ++i) {
            if (std::abs(ARRAYF(savf, i)) <= dlss_.seth) continue;
            if (i == j) continue;
            if (k > liwk) goto LABEL_210;
            ARRAYF(iwk, k) = i;
            k++;
        }
        ARRAYF(iwk, dlss_.ipian + j) = k + 1 - dlss_.ipjan;
    }
    goto LABEL_140;
//
// MOSS = 2.  Compute structure from results of N + 1 calls to F. -------
LABEL_100:
    k = dlss_.ipjan;
    ARRAYF(iwk, dlss_.ipian) = 1;
    (*f)(neq, dls1_.tn, y, savf, user_data);
    for (j = 1; j <= dls1_.n; ++j) {
        if (k > liwk) goto LABEL_210;
        ARRAYF(iwk, k) = j;
        k++;
        yj = ARRAYF(y, j);
        erwt = 1.0 / ARRAYF(ewt, j);
        dyj = std::abs(erwt) * SIGN(yj);
        ARRAYF(y, j) = yj + dyj;
        (*f)(neq, dls1_.tn, y, ftem, user_data);
        ARRAYF(y, j) = yj;
        for (i = 1; i <= dls1_.n; ++i) {
            dq = (ARRAYF(ftem, i) - ARRAYF(savf, i)) / dyj;
            if (std::abs(dq) <= dlss_.seth) continue;
            if (i == j) continue;
            if (k > liwk) goto LABEL_210;
            ARRAYF(iwk, k) = i;
            k++;
        }
        ARRAYF(iwk, dlss_.ipian + j) = k + 1 - dlss_.ipjan;
    }
//
LABEL_140:
    if (dlss_.moss == 0 || dlss_.istatc != 1) goto LABEL_150;
// If ISTATE = 1 and MOSS .ne. 0, restore Y from YH. --------------------
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = ARRAYF(yh, i);
    }
LABEL_150:
    dlss_.nnz = ARRAYF(iwk, dlss_.ipian + dls1_.n) - 1;
    lenigp = 0;
    dlss_.ipigp = dlss_.ipjan + dlss_.nnz;
    if (dls1_.miter != 2) goto LABEL_160;
//
// Compute grouping of column indices (MITER = 2). ----------------------
    maxg        = np1;
    dlss_.ipjgp = dlss_.ipjan + dlss_.nnz;
    dlss_.ibjgp = dlss_.ipjgp - 1;
    dlss_.ipigp = dlss_.ipjgp + dls1_.n;
    iptt1       = dlss_.ipigp + np1;
    iptt2       = iptt1 + dls1_.n;
    dlss_.lreq  = iptt2 + dls1_.n - 1;
    if (dlss_.lreq > liwk) goto LABEL_220;
    JGROUP(dls1_.n, &ARRAYF(iwk, dlss_.ipian), &ARRAYF(iwk, dlss_.ipjan), maxg, dlss_.ngp, &ARRAYF(iwk, dlss_.ipigp), 
        &ARRAYF(iwk, dlss_.ipjgp), &ARRAYF(iwk, iptt1), &ARRAYF(iwk, iptt2), ier);
    if (ier != 0) goto LABEL_220;
    lenigp = dlss_.ngp + 1;
//
// Compute new ordering of rows/columns of Jacobian. --------------------
LABEL_160:
    dlss_.ipr = dlss_.ipigp + lenigp;
    dlss_.ipc = dlss_.ipr;
    dlss_.ipic = dlss_.ipc  + dls1_.n;
    dlss_.ipisp = dlss_.ipic + dls1_.n;
    dlss_.iprsp = (dlss_.ipisp - 2) / dlss_.lrat + 2;
    dlss_.iesp  = dlss_.lenwk + 1 - dlss_.iprsp;
    if (dlss_.iesp < 0) goto LABEL_230;
    ibr = dlss_.ipr - 1;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(iwk, ibr + i) = i;
    }
    dlss_.nsp = liwk + 1 - dlss_.ipisp;
    ODRV(dls1_.n, &ARRAYF(iwk, dlss_.ipian), &ARRAYF(iwk, dlss_.ipjan), wk, &ARRAYF(iwk, dlss_.ipr), &ARRAYF(iwk, dlss_.ipic), 
        dlss_.nsp, &ARRAYF(iwk, dlss_.ipisp), 1, dlss_.iys);
    if (dlss_.iys == 11 * dls1_.n + 1) goto LABEL_240;
    if (dlss_.iys != 0) goto LABEL_230;
//
// Reorder JAN and do symbolic LU factorization of matrix. --------------
    dlss_.ipa = dlss_.lenwk + 1 - dlss_.nnz;
    dlss_.nsp = dlss_.ipa - dlss_.iprsp;
    dlss_.lreq = std::max(12 * dls1_.n / dlss_.lrat, 6 * dls1_.n / dlss_.lrat + 2 * dls1_.n + dlss_.nnz) + 3;
    dlss_.lreq = dlss_.lreq + dlss_.iprsp - 1 + dlss_.nnz;
    if (dlss_.lreq > dlss_.lenwk) goto LABEL_250;
    dlss_.iba = dlss_.ipa - 1;
    for (i = 1; i <= dlss_.nnz; ++i) {
        ARRAYF(wk, dlss_.iba + i) = 0.0;
    }
    dlss_.ipisp = dlss_.lrat * (dlss_.iprsp - 1) + 1;
    CDRV(dls1_.n, &ARRAYF(iwk, dlss_.ipr), &ARRAYF(iwk, dlss_.ipc), &ARRAYF(iwk, dlss_.ipic), &ARRAYF(iwk, dlss_.ipian), &ARRAYF(iwk, dlss_.ipjan),
        &ARRAYF(wk, dlss_.ipa), &ARRAYF(wk, dlss_.ipa), &ARRAYF(wk, dlss_.ipa), dlss_.nsp, &ARRAYF(iwk, dlss_.ipisp), &ARRAYF(wk, dlss_.iprsp), dlss_.iesp, 5, dlss_.iys);
    dlss_.lreq = dlss_.lenwk - dlss_.iesp;
    if (dlss_.iys == 10 * dls1_.n + 1) goto LABEL_250;
    if (dlss_.iys != 0) goto LABEL_260;
    ipil      = dlss_.ipisp;
    ipiu      = ipil + 2 * dls1_.n + 1;
    dlss_.nzu = ARRAYF(iwk, ipil + dls1_.n) - ARRAYF(iwk, ipil);
    dlss_.nzl = ARRAYF(iwk, ipiu + dls1_.n) - ARRAYF(iwk, ipiu);
    if (dlss_.lrat > 1) goto LABEL_190;
    ADJLR(dls1_.n, &ARRAYF(iwk, dlss_.ipisp), ldif);
    dlss_.lreq += ldif;
LABEL_190:
    if (dlss_.lrat == 2 && dlss_.nnz == dls1_.n) dlss_.lreq++;
    dlss_.nsp = dlss_.nsp + dlss_.lreq - dlss_.lenwk;
    dlss_.ipa = dlss_.lreq + 1 - dlss_.nnz;
    dlss_.iba = dlss_.ipa  - 1;
    ipper = 0;
    return;
//
LABEL_210:
    ipper = -1;
    dlss_.lreq = 2 + (2 * dls1_.n + 1) / dlss_.lrat;
    dlss_.lreq = std::max(dlss_.lenwk + 1, dlss_.lreq);
    return;
//
LABEL_220:
    ipper = -2;
    dlss_.lreq = (dlss_.lreq - 1) / dlss_.lrat + 1;
    return;
//
LABEL_230:
    ipper = -3;
    CNTNZU(dls1_.n, &ARRAYF(iwk, dlss_.ipian), &ARRAYF(iwk, dlss_.ipjan), nzsut);
    dlss_.lreq = dlss_.lenwk - dlss_.iesp + (3 * dls1_.n + 4 * nzsut - 1) / dlss_.lrat + 1;
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
    dlss_.lreq = dlss_.lenwk;
    return;
}

/**
 * @fn JGROUP
 * @brief
 *  This subroutine constructs groupings of the column indices of
 *  the Jacobian matrix, used in the numerical evaluation of the
 *  Jacobian by finite differences.
 * @date 2025/02/07 kawasaki
 */
void Odepack::JGROUP(const int n, int *ia, int *ja, const int maxg, int &ngrp, int *igp, int *jgp, int *incl, int *jdone, int &ier)
{
    int i, j, k, kmin, kmax, ncol, ng;
    bool is_goto_50 = false;
//
    ier = 0;
    for (j = 1; j <= n; ++j) {
        ARRAYF(jdone, j) = 0;
    }
    ncol = 1;
    for (ng = 1; ng <= maxg; ++ng) {
        ARRAYF(igp, ng) = ncol;
        for (i = 1; i <= n; ++i) {
            ARRAYF(incl, i) = 0;
        }
        is_goto_50 = false;
        for (j = 1; j <= n; ++j) {
// Reject column J if it is already in a group.--------------------------
            if (ARRAYF(jdone, j) == 1) continue;
            kmin = ARRAYF(ia, j);
            kmax = ARRAYF(ia, j + 1) - 1;
            for (k = kmin; k <= kmax; ++k) {
// Reject column J if it overlaps any column already in this group.------
                i = ARRAYF(ja, k);
                if (ARRAYF(incl, i) == 1) {
                    is_goto_50 = true;
                    break;
                };
            }
            if (is_goto_50) continue;
// Accept column J into group NG.----------------------------------------
            ARRAYF(jgp, ncol) = j;
            ncol++;
            ARRAYF(jdone, j) = 1;
            for (k = kmin; k <= kmax; ++k) {
                i = ARRAYF(ja, k);
                ARRAYF(incl, i) = 1;
            }
        }
// Stop if this group is empty (grouping is complete).-------------------
        if (ncol == ARRAYF(igp, ng)) goto LABEL_70;
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
 */
void Odepack::ADJLR(const int n, const int *isp, int &ldif)
{
    int ip, jlmax, jumax, lnfc, lsfc, nzlu;
//
    ip = 2 * n + 1;
// Get JLMAX = IJL(N) and JUMAX = IJU(N) (sizes of JL and JU). ----------
    jlmax = ARRAYF(isp, ip);
    jumax = ARRAYF(isp, ip + ip);
// NZLU = (size of L) + (size of U) = (IL(N+1)-IL(1)) + (IU(N+1)-IU(1)).
    nzlu = ARRAYF(isp, n + 1) - ARRAYF(isp, 1) + ARRAYF(isp, ip + n + 1) - ARRAYF(isp, ip + 1);
    lsfc = 12 * n + 3 + 2 * std::max(jlmax, jumax);
    lnfc = 9 * n + 2 + jlmax + jumax + nzlu;
    ldif = std::max(0, lsfc - lnfc);
    return;
}

/**
 * @fn CNTNZU
 */
void Odepack::CNTNZU(const int n, const int *ia, const int *ja, int &nzsut)
{
    int ii, jj, j, jmin, jmax, k, kmin, kmax, num;
    bool is_goto_40 = false;
//
    num = 0;
    for (ii = 1; ii <= n; ++ii) {
        jmin = ARRAYF(ia, ii);
        jmax = ARRAYF(ia, ii + 1) - 1;
        if (jmin > jmax) continue;
        for (j = jmin; j <= jmax; ++j) {
            if (ARRAYF(ja, j) < ii) {
                goto LABEL_10;
            } else if (ARRAYF(ja, j) == ii) {
                continue;
            } else {
                goto LABEL_30;
            }
LABEL_10:
            jj   = ARRAYF(ja,  j);
            kmin = ARRAYF(ia, jj);
            kmax = ARRAYF(ia, jj + 1) - 1;
            if (kmin > kmax) goto LABEL_30;
            is_goto_40 = false;
            for (k = kmin; k <= kmax; ++k) {
                if (ARRAYF(ja, k) == ii) {
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
 * @brief DPRJS is called to compute and process the matrix
 * @date 2025/02/08 kawasaki
 */
void Odepack::DPRJS(const int neq, double *y, double *yh, const int nyh, const double *ewt, 
    double *ftem, double *savf, double *wk, int *iwk, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN2 jac, 
    void *user_data)
{
#ifndef YH
#define YH(i, j) MATF(yh, nyh, i, j)
#endif
//
    int i, imul, j, jj, jok, jmax, jmin, k, kmax, kmin, ng;
    double con, di, fac, hl0, pij, r, r0, rcon, rcont, srur;
//
    hl0 = dls1_.h * dls1_.el0;
    con = - hl0;
    if (dls1_.miter == 3) goto LABEL_300;
// See whether J should be reevaluated (JOK = 0) or not (JOK = 1). ------
    jok = 1;
    if (dls1_.nst == 0 || dls1_.nst >= (dlss_.nslj + dlss_.msbj)) jok = 0;
    if (dls1_.icf == 1 && std::abs(dls1_.rc - 1.0) < dlss_.ccmxj) jok = 0;
    if (dls1_.icf == 2) jok = 0;
    if (jok == 1) goto LABEL_250;
//
// MITER = 1 or 2, and the Jacobian is to be reevaluated. ---------------
LABEL_20:
    dls1_.jcur = 1;
    dls1_.nje++;
    dlss_.nslj = dls1_.nst;
    dlss_.iplost = 0;
    dlss_.conmin = std::abs(con);
    if (dls1_.miter == 1) {
        goto LABEL_100;
    } else if (dls1_.miter == 2) {
        goto LABEL_200;
    }
//
// If MITER = 1, call JAC, multiply by scalar, and add identity. --------
LABEL_100:
    kmin = ARRAYF(iwk, dlss_.ipian);
    for (j = 1; j <= dls1_.n; ++j) {
        kmax = ARRAYF(iwk, dlss_.ipian + j) - 1;
        for (i = 1; i <= dls1_.n; ++i) {
            ARRAYF(ftem, i) = 0.0;
        }
        (*jac)(neq, dls1_.tn, y, j, &ARRAYF(iwk, dlss_.ipian), &ARRAYF(iwk, dlss_.ipjan), ftem, user_data);
        for (k = kmin; k <= kmax; ++k) {
            i = ARRAYF(iwk, dlss_.ibjan + k);
            ARRAYF(wk, dlss_.iba + k) = ARRAYF(ftem, i) * con;
            if (i == j) ARRAYF(wk, dlss_.iba + k) += 1.0;
        }
        kmin = kmax + 1;
    }
    goto LABEL_290;
//
// If MITER = 2, make NGP calls to F to approximate J and P. ------------
LABEL_200:
    fac = DVNORM(dls1_.n, savf, ewt);
    r0 = 1000.0 * std::abs(dls1_.h) * dls1_.uround * static_cast<double>(dls1_.n) * fac;
    if (r0 == 0.0) r0 = 1.0;
    srur = ARRAYF(wk, 1);
    jmin = ARRAYF(iwk, dlss_.ipigp);
    for (ng = 1; ng <= dlss_.ngp; ++ng) {
        jmax = ARRAYF(iwk, dlss_.ipigp + ng) - 1;
        for (j = jmin; j <= jmax; ++j) {
            jj = ARRAYF(iwk, dlss_.ibjgp + j);
            r  = std::max(srur * std::abs(ARRAYF(y, jj)), r0 / ARRAYF(ewt, jj));
            ARRAYF(y, jj) += r;
        }
        (*f)(neq, dls1_.tn, y, ftem, user_data);
        for (j = jmin; j <= jmax; ++j) {
            jj = ARRAYF(iwk, dlss_.ibjgp + j);
            ARRAYF(y, jj) = YH(jj, 1);
            r = std::max(srur * std::abs(ARRAYF(y, jj)), r0 / ARRAYF(ewt, jj));
            fac  = - hl0 / r;
            kmin = ARRAYF(iwk, dlss_.ibian + jj);
            kmax = ARRAYF(iwk, dlss_.ibian + jj + 1) - 1;
            for (k = kmin; k <= kmax; ++k) {
                i = ARRAYF(iwk, dlss_.ibjan + k);
                ARRAYF(wk, dlss_.iba + k) = (ARRAYF(ftem, i) - ARRAYF(savf, i)) * fac;
                if (i == jj) ARRAYF(wk, dlss_.iba + k) += 1.0;
            }
        }
        jmin = jmax + 1;
    }
    dls1_.nfe += dlss_.ngp;
    goto LABEL_290;
//
// If JOK = 1, reconstruct new P from old P. ----------------------------
LABEL_250:
    dls1_.jcur = 0;
    rcon = con / dlss_.con0;
    rcont = std::abs(con) / dlss_.conmin;
    if (rcont > dlss_.rbig && dlss_.iplost == 1) goto LABEL_20;
    kmin = ARRAYF(iwk, dlss_.ipian);
    for (j = 1; j <= dls1_.n; ++j) {
        kmax = ARRAYF(iwk, dlss_.ipian + j) - 1;
        for (k = kmin; k <= kmax; ++k) {
            i = ARRAYF(iwk, dlss_.ibjan + k);
            pij = ARRAYF(wk, dlss_.iba + k);
            if (i != j) goto LABEL_260;
            pij = pij - 1.0;
            if (std::abs(pij) >= dlss_.psmall) goto LABEL_260;
            dlss_.iplost = 1;
            dlss_.conmin = std::min(std::abs(dlss_.con0), dlss_.conmin);
LABEL_260:
            pij *= rcon;
            if (i == j) pij += 1.0;
            ARRAYF(wk, dlss_.iba + k) = pij;
        }
        kmin = kmax + 1;
    }
//
// Do numerical factorization of P matrix. ------------------------------
LABEL_290:
    dlss_.nlu++;
    dlss_.con0 = con;
    dls1_.ierpj = 0;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(ftem, i) = 0.0;
    }
    CDRV(dls1_.n, &ARRAYF(iwk, dlss_.ipr), &ARRAYF(iwk, dlss_.ipc), &ARRAYF(iwk, dlss_.ipic), &ARRAYF(iwk, dlss_.ipian), &ARRAYF(iwk, dlss_.ipjan),
        &ARRAYF(wk, dlss_.ipa), ftem, ftem, dlss_.nsp, &ARRAYF(iwk, dlss_.ipisp), &ARRAYF(wk, dlss_.iprsp), dlss_.iesp, 2, dlss_.iys);
    if (dlss_.iys == 0) return;
    imul = (dlss_.iys - 1) / dls1_.n;
    dls1_.ierpj = -2;
    if (imul == 8) dls1_.ierpj = 1;
    if (imul == 10) dls1_.ierpj = -1;
    return;
//
// If MITER = 3, construct a diagonal approximation to J and P. ---------
LABEL_300:
    dls1_.jcur = 1;
    dls1_.nje++;
    ARRAYF(wk, 2) = hl0;
    dls1_.ierpj = 0;
    r = dls1_.el0 * 0.1;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = ARRAYF(y, i) + r * (dls1_.h * ARRAYF(savf, i) - YH(i, 2));
    }
    (*f)(neq, dls1_.tn, y, &ARRAYF(wk, 3), user_data);
    dls1_.nfe++;
    for (i = 1; i <= dls1_.n; ++i) {
        r0 = dls1_.h * ARRAYF(savf, i) - YH(i, 2);
        di = 0.1 * r0 - dls1_.h * (ARRAYF(wk, i + 2) - ARRAYF(savf, i));
        ARRAYF(wk, i + 2) = 1.0;
        if (std::abs(r0) < dls1_.uround / ARRAYF(ewt, i)) continue;
        if (std::abs(di) == 0.0) goto LABEL_330;
        ARRAYF(wk, i + 2) = 0.1 * r0 / di;
    }
    return;
LABEL_330:
    dls1_.ierpj = 2;
    return;
//
#ifdef YH
#undef YH
#endif
}

/**
 * @fn DSOLSS
 * @brief This routine manages the solution of the linear system arising from a chord iteration.
 * @date 2025/02/08 kawasaki
 */
void Odepack::DSOLSS(double *wk, int *iwk, double *x, double *tem)
{
    int i;
    double di, hl0, phl0, r;
//
    dls1_.iersl = 0;
    if (dls1_.miter == 1 || dls1_.miter == 2) {
        goto LABEL_100;
    } else if (dls1_.miter == 3) {
        goto LABEL_300;
    }
LABEL_100:
    CDRV(dls1_.n, &ARRAYF(iwk, dlss_.ipr), &ARRAYF(iwk, dlss_.ipc), &ARRAYF(iwk, dlss_.ipic), &ARRAYF(iwk, dlss_.ipian), &ARRAYF(iwk, dlss_.ipjan),
        &ARRAYF(wk, dlss_.ipa), x, x, dlss_.nsp, &ARRAYF(iwk, dlss_.ipisp), &ARRAYF(wk, dlss_.iprsp), dlss_.iesp, 4, dls1_.iersl);
    if (dls1_.iersl != 0) dls1_.iersl = -1;
    return;
//
LABEL_300:
    phl0 = ARRAYF(wk, 2);
    hl0 = dls1_.h * dls1_.el0;
    ARRAYF(wk, 2) = hl0;
    if (hl0 == phl0) goto LABEL_330;
    r = hl0 / phl0;
    for (i = 1; i <= dls1_.n; ++i) {
        di = 1.0 - r * (1.0 - 1.0 / ARRAYF(wk, i + 2));
        if (std::abs(di) == 0.0) goto LABEL_390;
        ARRAYF(wk, i + 2) = 1.0 / di;
    }
LABEL_330:
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(x, i) *= ARRAYF(wk, i + 2);
    }
    return;
LABEL_390:
    dls1_.iersl = 1;
    return;
}

/**
 * @fn DSRCMS
 * @brief 
 * @date kawasaki 2025/02/08
 */
void Odepack::DSRCMS(double *rsav, int *isav, const int job)
{
    int i;
    const int lenrls = 218;
    const int lenils = 37;
    const int lenrss = 6;
    const int leniss = 34;
//
    if (job == 2) goto LABEL_100;
    for (i = 1; i <= lenrls; ++i) {
        ARRAYF(rsav, i) = ARRAYF(dls1_.rls, i);
    }
    for (i = 1; i <= lenrss; ++i) {
        ARRAYF(rsav, lenrls + i) = ARRAYF(dls1_.rlss, i);
    }
//
    for (i = 1; i <= lenils; ++i) {
        ARRAYF(isav, i) = ARRAYF(dls1_.ils, i);
    }
    for (i = 1; i <= leniss; ++i) {
        ARRAYF(isav, lenils + i) = ARRAYF(dls1_.ilss, i);
    }
//
    return;
//
LABEL_100:
    for (i = 1; i <= lenrls; ++i) {
        ARRAYF(dls1_.rls, i) = ARRAYF(rsav, i);
    }
    for (i = 1; i <= lenrss; ++i) {
        ARRAYF(dls1_.rlss, i) = ARRAYF(rsav, lenrls + i);
    }
//
    for (i = 1; i <= lenils; ++i) {
        ARRAYF(dls1_.ils, i) = ARRAYF(isav, i);
    }
    for (i = 1; i <= leniss; ++i) {
        ARRAYF(dls1_.ilss, i) = ARRAYF(isav, lenils + i);
    }
//
    return;
}

/**
 * @fn ODRV
 * @brief driver for sparse matrix reordering routines
 * @date 2025/02/08 kawasaki
 */
void Odepack::ODRV(const int n, int *ia, int *ja, double *a, int *p, int *ip, int nsp, int *isp, int path, int &flag)
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
    max  = (nsp - n) / 2;
    v    = 1;
    l    = v + max;
    head = l + max;
    next = head + n;
    if (max < n) goto LABEL_110;
//
    MD(n, ia, ja, max, &ARRAYF(isp, v), &ARRAYF(isp, l), &ARRAYF(isp, head), p, ip, &ARRAYF(isp, v), flag);
    if (flag != 0) goto LABEL_100;
//
//----allocate storage and symmetrically reorder matrix
LABEL_1:
    if ((path - 2) * (path - 3) * (path - 4) * (path - 5) != 0) goto LABEL_2;
    tmp = (nsp + 1) - n;
    q   = tmp - (ARRAYF(ia, n + 1) - 1);
    if (q < 1) goto LABEL_110;
//
    dflag = (path == 4 || path == 5);
    SRO(n, ip, ia, ja, a, &ARRAYF(isp, tmp), &ARRAYF(isp, q), dflag);
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
 * @fn MD
 * @brief minimum degree algorithm (based on element model)
 * @date 2025/02/08 kawasaki
 */
void Odepack::MD(const int n, int *ia, int *ja, int max, int *v, int *l, int *head, int *last, int *next, int *mark, int &flag)
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
    if (ARRAYF(head, dmin) > 0) goto LABEL_3;
    dmin++;
    goto LABEL_2;
//
//------remove vertex vk of minimum degree from degree list
LABEL_3:
    vk = ARRAYF(head, dmin);
    ARRAYF(head, dmin) = ARRAYF(next, vk);
    if (ARRAYF(head, dmin) > 0) ARRAYF(last, ARRAYF(head, dmin)) = -dmin;
//
//------number vertex vk, adjust tag, and tag vk
    k++;
    ARRAYF(next,  vk) = -k;
    ARRAYF(last, *ek) = dmin - 1;
    tag += ARRAYF(last, *ek);
    ARRAYF(mark,  vk) = tag;
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
        ARRAYF(next, k) = - ARRAYF(next, k);
        ARRAYF(last, ARRAYF(next, k)) = k;
    }
//
    return;
}

/**
 * @fn MDI
 * @brief mdi -- initialization
 * @date 2025/02/08 kawasaki
 */
void Odepack::MDI(const int n, int *ia, int *ja, int max, int *v, int *l, int *head, int *last, int *next, int *mark, int tag, int &flag)
{
    int sfs, j, jmin, jmax, vi, dvi, vj, lvk, k, kmax, nextvi;
    bool is_goto_5 = false;
//
//----initialize degrees, element lists, and degree lists
    for (vi = 1; vi <= n; ++vi) {
        ARRAYF(mark, vi) = 1;
        ARRAYF(l, vi)    = 0;
        ARRAYF(head, vi) = 0;
    }
    sfs = n + 1;
//
//----create nonzero structure
//----for each nonzero entry a(vi,vj)
    for (vi = 1; vi <= n; ++vi) {
        jmin = ARRAYF(ia, vi);
        jmax = ARRAYF(ia, vi + 1) - 1;
        if (jmin > jmax) continue;
        for (j = jmin; j <= jmax; ++j) {
            vj = ARRAYF(ja, j);
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
            kmax = ARRAYF(mark, vi) - 1;
            if (kmax == 0) goto LABEL_4;
            is_goto_5 = false;
            for (k = 1; k <= kmax; ++k) {
                lvk = ARRAYF(l, lvk);
                if (ARRAYF(v, lvk) == vj) {
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
            ARRAYF(mark, vi) += 1;
            ARRAYF(v, sfs) = vj;
            ARRAYF(l, sfs) = ARRAYF(l, vi);
            ARRAYF(l, vi) = sfs;
            sfs++;
//
//------enter vi in element list for vj
            ARRAYF(mark, vj) += 1;
            ARRAYF(v, sfs) = vi;
            ARRAYF(l, sfs) = ARRAYF(l, vj);
            ARRAYF(l, vj) = sfs;
            sfs++;
        }
    }
//
//----create degree lists and initialize mark vector
    for (vi = 1; vi <= n; ++vi) {
        dvi = ARRAYF(mark, vi);
        ARRAYF(next, vi) = ARRAYF(head, dvi);
        ARRAYF(head, dvi) = vi;
        ARRAYF(last, vi) = -dvi;
        nextvi = ARRAYF(next,  vi);
        if (nextvi > 0) ARRAYF(last, nextvi) = vi;
        ARRAYF(mark, vi) = tag;
    }
//
    return;
// ** error-  insufficient storage
LABEL_101:
    flag = 9 * n + vi;
    return;
}

/**
 * @fn MDM
 * @brief mdm -- form element from uneliminated neighbors of vk
 * @date 2025/02/08 kawasaki
 */
void Odepack::MDM(const int vk, int &tail, int *v, int *l, int *last, int *next, int *mark)
{
    int tag, s, ls, vs, b, lb, vb, blp, blpmax;
    int *es = &vs;
//
//----initialize tag and list of uneliminated neighbors
    tag = ARRAYF(mark, vk);
    tail = vk;
//
//----for each vertex/element vs/es in element list of vk
    ls = ARRAYF(l, vk);
LABEL_1:
    s = ls;
    if (s == 0) goto LABEL_5;
    ls = ARRAYF(l, s);
    vs = ARRAYF(v, s);
    if (ARRAYF(next, vs) < 0) goto LABEL_2;
//
//------if vs is uneliminated vertex, then tag and append to list of
//------uneliminated neighbors
    ARRAYF(mark, vs) = tag;
    ARRAYF(l, tail) = s;
    tail = s;
    goto LABEL_4;
//
//------if es is active element, then ...
//--------for each vertex vb in boundary list of element es
LABEL_2:
    lb = ARRAYF(l, *es);
    blpmax = ARRAYF(last, *es);
    for (blp = 1; blp <= blpmax; ++blp) {
        b = lb;
        lb = ARRAYF(l, b);
        vb = ARRAYF(v, b);
//
//----------if vb is untagged vertex, then tag and append to list of
//----------uneliminated neighbors
        if (ARRAYF(mark, vb) >= tag) continue;
        ARRAYF(mark, vb) = tag;
        ARRAYF(l, tail) = b;
        tail = b;
    }
//
//--------mark es inactive
    ARRAYF(mark, *es) = tag;
//
LABEL_4:
    goto LABEL_1;
//
//----terminate list of uneliminated neighbors
LABEL_5:
    ARRAYF(l, tail) = 0;
//
    return;
}

/**
 * @fn MDP
 * @brief mdp -- purge inactive elements and do mass elimination
 * @date 2025/02/08 kawasaki
 */
void Odepack::MDP(int &k, int ek, int &tail, int *v, int *l, int *head, int *last, int *next, int *mark)
{
    int tag, free, li, vi, lvi, evi, s, ls, es, ilp, ilpmax, i;
//
//----initialize tag
    tag = ARRAYF(mark, ek);
//
//----for each vertex vi in ek
    li = ek;
    ilpmax = ARRAYF(last, ek);
    if (ilpmax <= 0) goto LABEL_12;
    for (ilp = 1; ilp <= ilpmax; ++ilp) {
        i = li;
        li = ARRAYF(l, i);
        vi = ARRAYF(v, li);
//
//------remove vi from degree list
        if (ARRAYF(last, vi) == 0) goto LABEL_3;
        if (ARRAYF(last, vi) > 0) goto LABEL_1;
        ARRAYF(head, -ARRAYF(last, vi)) = ARRAYF(next, vi);
        goto LABEL_2;
LABEL_1:
        ARRAYF(next, ARRAYF(last, vi)) = ARRAYF(next, vi);
LABEL_2:
        if (ARRAYF(next, vi) > 0) ARRAYF(last, ARRAYF(next, vi)) = ARRAYF(last, vi);
//
//------remove inactive items from element list of vi
LABEL_3:
        ls = vi;
LABEL_4:
        s = ls;
        ls = ARRAYF(l, s);
        if (ls == 0) goto LABEL_6;
        es = ARRAYF(v, ls);
        if (ARRAYF(mark, es) < tag) goto LABEL_5;
        free = ls;
        ARRAYF(l, s) = ARRAYF(l, ls);
        ls = s;
LABEL_5:
        goto LABEL_4;
//
//------if vi is interior vertex, then remove from list and eliminate
LABEL_6:
        lvi = ARRAYF(l, vi);
        if (lvi != 0) goto LABEL_7;
        ARRAYF(l, i) = ARRAYF(l, li);
        li = i;
//
        k++;
        ARRAYF(next, vi) = -k;
        ARRAYF(last, ek) -= 1;
        continue;
//
//------else ...
//--------classify vertex vi
LABEL_7:
        if (ARRAYF(l, lvi) != 0) goto LABEL_9;
        evi = ARRAYF(v, lvi);
        if (ARRAYF(next, evi) >= 0) goto LABEL_9;
        if (ARRAYF(mark, evi) < 0) goto LABEL_8;
//
//----------if vi is prototype vertex, then mark as such, initialize
//----------overlap count for corresponding element, and move vi to end
//----------of boundary list
        ARRAYF(last, vi) = evi;
        ARRAYF(mark, evi) = -1;
        ARRAYF(l, tail) = li;
        tail = li;
        ARRAYF(l, i) = ARRAYF(l, li);
        li = i;
        goto LABEL_10;
//
//----------else if vi is duplicate vertex, then mark as such and adjust
//----------overlap count for corresponding element
LABEL_8:
        ARRAYF(last, vi) = 0;
        ARRAYF(mark, evi) -= 1;
        goto LABEL_10;
//
//----------else mark vi to compute degree
LABEL_9:
        ARRAYF(last, vi) = -ek;
//
//--------insert ek in element list of vi
LABEL_10:
        ARRAYF(v, free) = ek;
        ARRAYF(l, free) = ARRAYF(l, vi);
        ARRAYF(l, vi) = free;
    }
//
//----terminate boundary list
LABEL_12:
    ARRAYF(l, tail) = 0;
//
    return;
}

/**
 * @fn MDU
 * @brief mdu -- update degrees of uneliminated vertices in ek
 * @date 2025/02/08 kawasaki
 */
void Odepack::MDU(int ek, int &dmin, int *v, int *l, int *head, int *last, int *next, int *mark)
{
    int tag, vi, evi, dvi, s, vs, b, vb, ilp, ilpmax, blp, blpmax, i;
    int *es = &vs;
//
//----initialize tag
    tag = ARRAYF(mark, ek) - ARRAYF(last, ek);
//
//----for each vertex vi in ek
    i = ek;
    ilpmax = ARRAYF(last, ek);
    if (ilpmax <= 0) goto LABEL_11;
    for (ilp = 1; ilp <= ilpmax; ++ilp) {
        i = ARRAYF(l, i);
        vi = ARRAYF(v, i);
        if (ARRAYF(last, vi) < 0) {
            goto LABEL_1;
        } else if (ARRAYF(last, vi) == 0) {
            continue;
        } else {
            goto LABEL_8;
        }
//
//------if vi neither prototype nor duplicate vertex, then merge elements
//------to compute degree
LABEL_1:
        tag++;
        dvi = ARRAYF(last, ek);
//
//--------for each vertex/element vs/es in element list of vi
        s = ARRAYF(l, vi);
LABEL_2:
        s = ARRAYF(l, s);
        if (s == 0) goto LABEL_9;
        vs = ARRAYF(v, s);
        if (ARRAYF(next, vs) < 0) goto LABEL_3;
//
//----------if vs is uneliminated vertex, then tag and adjust degree
        ARRAYF(mark, vs) = tag;
        dvi++;
        goto LABEL_5;
//
//----------if es is active element, then expand
//------------check for outmatched vertex
LABEL_3:
        if (ARRAYF(mark, *es) < 0) goto LABEL_6;
//
//------------for each vertex vb in es
        b = *es;
        blpmax = ARRAYF(last, *es);
        for (blp = 1; blp <= blpmax; ++blp) {
            b = ARRAYF(l, b);
            vb = ARRAYF(v, b);
//
//--------------if vb is untagged, then tag and adjust degree
            if (ARRAYF(mark, vb) >= tag) continue;
            ARRAYF(mark, vb) = tag;
            dvi++;
        }
//
LABEL_5:
        goto LABEL_2;
//
//------else if vi is outmatched vertex, then adjust overlaps but do not
//------compute degree
LABEL_6:
        ARRAYF(last, vi)  = 0;
        ARRAYF(mark, *es) -= 1;
LABEL_7:
        s = ARRAYF(l, s);
        if (s == 0) continue;
        *es = ARRAYF(v, s);
        if (ARRAYF(mark, *es) < 0) ARRAYF(mark, *es) -= 1;
        goto LABEL_7;
//
//------else if vi is prototype vertex, then calculate degree by
//------inclusion/exclusion and reset overlap count
LABEL_8:
        evi = ARRAYF(last, vi);
        dvi = ARRAYF(last, ek) + ARRAYF(last, evi) + ARRAYF(mark, evi);
        ARRAYF(mark, evi) = 0;
//
//------insert vi in appropriate degree list
LABEL_9:
        ARRAYF(next, vi) = ARRAYF(head, dvi);
        ARRAYF(head, dvi) = vi;
        ARRAYF(last, vi) = -dvi;
        if (ARRAYF(next, vi) > 0) ARRAYF(last, ARRAYF(next, vi)) = vi;
        if (dvi < dmin) dmin = dvi;
    }
//
LABEL_11:
    return;
}

/**
 * @fn SRO
 * @brief sro -- symmetric reordering of sparse symmetric matrix
 * @date 2025/02/08 kawasaki
 */
void Odepack::SRO(const int n, int *ip, int *ia, int *ja, double *a, int *q, int *r, bool dflag)
{
    int i, j, jmin, jmax, jdummy, k, ilast, jak;
    double ak;
//
//
//--phase 1 -- find row in which to store each nonzero
//----initialize count of nonzeroes to be stored in each row
    for (i = 1; i <= n; ++i) {
        ARRAYF(q, i) = 0;
    }
//
//----for each nonzero element a(j)
    for (i = 1; i <= n; ++i) {
        jmin = ARRAYF(ia, i);
        jmax = ARRAYF(ia, i + 1) - 1;
        if (jmin > jmax) continue;
        for (j = jmin; j <= jmax; ++j) {
//
//--------find row (=r(j)) and column (=ja(j)) in which to store a(j)
            k = ARRAYF(ja, j);
            if (ARRAYF(ip, k)  < ARRAYF(ip, i)) ARRAYF(ja, j) = i;
            if (ARRAYF(ip, k) >= ARRAYF(ip, i)) k = i;
            ARRAYF(r, j) = k;
//
//--------.... and increment count of nonzeroes (=q(r(j)) in that row
            ARRAYF(q, k) += 1;
        }
    }
//
//
//--phase 2 -- find new ia and permutation to apply to (ja,a)
//----determine pointers to delimit rows in permuted (ja,a)
    for (i = 1; i <= n; ++i) {
        ARRAYF(ia, i + 1) = ARRAYF(ia,     i) + ARRAYF(q, i);
        ARRAYF( q,     i) = ARRAYF(ia, i + 1);
    }
//
//----determine where each (ja(j),a(j)) is stored in permuted (ja,a)
//----for each nonzero element (in reverse order)
    ilast = 0;
    jmin  = ARRAYF(ia,     1);
    jmax  = ARRAYF(ia, n + 1) - 1;
    j     = jmax;
    for (jdummy = jmin; jdummy <= jmax; ++jdummy) {
        i = ARRAYF(r, j);
        if (!dflag || ARRAYF(ja, j) != i || i == ilast) goto LABEL_5;
//
//------if dflag, then put diagonal nonzero at beginning of row
        ARRAYF(r, j) = ARRAYF(ia, i);
        ilast        = i;
        continue;
//
//------put (off-diagonal) nonzero in last unused location in row
LABEL_5:
        ARRAYF(q, i) -= 1;
        ARRAYF(r, j)  = ARRAYF(q, i);
//
        j--;
    }
//
//
//--phase 3 -- permute (ja,a) to upper triangular form (wrt new ordering)
    for (j = jmin; j <= jmax; ++j) {
LABEL_7:
        if (ARRAYF(r, j) == j) continue;
        k             = ARRAYF(r, j);
        ARRAYF(r, j)  = ARRAYF(r, k);
        ARRAYF(r, k)  = k;
        jak           = ARRAYF(ja, k);
        ARRAYF(ja, k) = ARRAYF(ja, j);
        ARRAYF(ja, j) = jak;
        ak            = ARRAYF(a, k);
        ARRAYF(a, k)  = ARRAYF(a, j);
        ARRAYF(a, j)  = ak;
        goto LABEL_7;
    }
//
    return;
}

void Odepack::CDRV(const int n, int *r, int *c, int *ic, int *ia, int *ja, 
    double *a, double *b, double *z, int nsp, int *isp, double *rsp, 
    int &esp, const int path, int &flag)
{
    int i, il, ijl, iu, iju, irl, jrl, jl, max, jlmax, ira, jra, irac, iru, jru, ju, jutmp, jumax, j, l, lmax;
    int d, u, q, row, tmp, ar, umax;
//
    const int lratio = 2;
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
        if (ARRAYF(c, i) != i) goto LABEL_2;
    }
    goto LABEL_3;
LABEL_2:
    ar = nsp + 1 - n;
    NROC(n, ic, ia, ja, a, &ARRAYF(isp, il), &ARRAYF(rsp, ar), &ARRAYF(isp, iu), flag);
    if (flag != 0) goto LABEL_100;
//
LABEL_3:
    NSFC(n, r, ic, ia, ja, 
        jlmax, &ARRAYF(isp, il), &ARRAYF(isp, jl), &ARRAYF(isp, ijl),
        jumax, &ARRAYF(isp, iu), &ARRAYF(isp, jutmp), &ARRAYF(isp, iju),
        &ARRAYF(isp, q), &ARRAYF(isp, ira), &ARRAYF(isp, jra), &ARRAYF(isp, irac),
        &ARRAYF(isp, irl), &ARRAYF(isp, jrl), &ARRAYF(isp, iru), &ARRAYF(isp, jru), flag);
    if (flag != 0) goto LABEL_100;
//  ******  move ju next to jl  *****************************************
    jlmax = ARRAYF(isp, ijl + n - 1);
    ju    = jl + jlmax;
    jumax = ARRAYF(isp, iju + n - 1);
    if (jumax <= 0) goto LABEL_5;
    for (j = 1; j <= jumax; ++j) {
        ARRAYF(isp, ju + j - 1) = ARRAYF(isp, jutmp + j - 1);
    }
//
//  ******  call remaining subroutines  *********************************
LABEL_5:
    jlmax = ARRAYF(isp, ijl + n - 1);
    ju    = jl + jlmax;
    jumax = ARRAYF(isp, iju + n - 1);
    l     = (ju + jumax - 2 + lratio) / lratio + 1;
    lmax  = ARRAYF(isp, il + n) - 1;
    d     = l + lmax;
    u     = d + n;
    row   = nsp + 1 - n;
    tmp   = row - n;
    umax  = tmp - u;
    esp   = umax - (ARRAYF(isp, iu + n) - 1);
//
    if ((path - 1) * (path - 2) != 0) goto LABEL_6;
    if (umax < 0) goto LABEL_110;
    NNFC(n, r, c, ic, ia, ja, a, z, b, 
        lmax, &ARRAYF(isp, il), &ARRAYF(isp, jl), &ARRAYF(isp, ijl), &ARRAYF(rsp, l), &ARRAYF(rsp, d),
        umax, &ARRAYF(isp, iu), &ARRAYF(isp, ju), &ARRAYF(isp, iju), &ARRAYF(rsp, u),
        &ARRAYF(rsp, row), &ARRAYF(rsp, tmp), &ARRAYF(isp, irl), &ARRAYF(isp, jrl), flag);
    if (flag != 0) goto LABEL_100;
//
LABEL_6:
    if ((path - 3) != 0) goto LABEL_7;
    NNSC(n, r, c, &ARRAYF(isp, il), &ARRAYF(isp, jl), &ARRAYF(isp, ijl), &ARRAYF(rsp, l),
         &ARRAYF(rsp, d), &ARRAYF(isp, iu), &ARRAYF(isp, ju), &ARRAYF(isp, iju), &ARRAYF(rsp, u),
         z, b, &ARRAYF(rsp, tmp));
//
LABEL_7:
    if ((path - 4) != 0) goto LABEL_8;
    NNTC(n, r, c, &ARRAYF(isp, il), &ARRAYF(isp, jl), &ARRAYF(isp, ijl), &ARRAYF(rsp, l),
        &ARRAYF(rsp, d), &ARRAYF(isp, iu), &ARRAYF(isp, ju), &ARRAYF(isp, iju), &ARRAYF(rsp, u),
        z, b, &ARRAYF(rsp, tmp));
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

void Odepack::NROC(const int n, int *ic, int *ia, int *ja, double *a, int *jar, 
    double *ar, int *p, int &flag)
{
    int i, j, jmin, jmax, newj, k;
//
//  ******  for each nonempty row  *******************************
    for (k = 1; k <= n; ++k) {
        jmin = ARRAYF(ia, k);
        jmax = ARRAYF(ia, k + 1) - 1;
        if (jmin > jmax) continue;
        ARRAYF(p, n + 1) = n + 1;
//  ******  insert each element in the list  *********************
        for (j = jmin; j <= jmax; ++j) {
            newj = ARRAYF(ic, ARRAYF(ja, j));
            i    = n + 1;
LABEL_1:
            if (ARRAYF(p, i) >= newj) goto LABEL_2;
            i = ARRAYF(p, i);
            goto LABEL_1;
LABEL_2:
            if (ARRAYF(p, i) == newj) goto LABEL_102;
            ARRAYF(  p, newj) = ARRAYF( p, i);
            ARRAYF(  p,    i) = newj;
            ARRAYF(jar, newj) = ARRAYF(ja, j);
            ARRAYF( ar, newj) = ARRAYF( a, j);
        }
//  ******  replace old row in ja and a  *************************
        i = n + 1;
        for (j = jmin; j <= jmax; ++j) {
            i             = ARRAYF(  p, i);
            ARRAYF(ja, j) = ARRAYF(jar, i);
            ARRAYF( a, j) = ARRAYF( ar, i);
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

void Odepack::NSFC(const int n, int *r, const int *ic, const int *ia, const int *ja, 
    const int jlmax, int *il, int *jl, int *ijl, const int jumax, int *iu, int *ju, 
    int *iju, int *q, int *ira, int *jra, int *irac, int *irl, int *jrl, int *iru, 
    int *jru, int &flag)
{
    int i, i1, np1, j, jlmin, jlptr, jumin, juptr, k, rk, iak, jaiak, luk, vj, qm, m, lastid, lasti,
        jmin, jmax, llong, jtmp, irll, cend, irul, rend, irai, jairai;
//
//  ******  initialize pointers  ****************************************
    np1           = n + 1;
    jlmin         = 1;
    jlptr         = 0;
    ARRAYF(il, 1) = 1;
    jumin         = 1;
    juptr         = 0;
    ARRAYF(iu, 1) = 1;
    for (k = 1; k <= n; ++k) {
        ARRAYF(irac, k) = 0;
        ARRAYF( jra, k) = 0;
        ARRAYF( jrl, k) = 0;
        ARRAYF( jru, k) = 0;
    }
//  ******  initialize column pointers for a  ***************************
    for (k = 1; k <= n; ++k) {
        rk  = ARRAYF( r,  k);
        iak = ARRAYF(ia, rk);
        if (iak >= ARRAYF(ia, rk + 1)) goto LABEL_101;
        jaiak = ARRAYF(ic, ARRAYF(ja, iak));
        if (jaiak > k) goto LABEL_105;
        ARRAYF( jra,     k) = ARRAYF(irac, jaiak);
        ARRAYF(irac, jaiak) = k;
        ARRAYF( ira,     k) = iak;
    }
//
//  ******  for each column of l and row of u  **************************
    for (k = 1; k <= n; ++k) {
//
//  ******  initialize q for computing kth column of l  *****************
        ARRAYF(q, np1) = np1;
        luk            = -1;
//  ******  by filling in kth column of a  ******************************
        vj = ARRAYF(irac, k);
        if (vj == 0) goto LABEL_5;
LABEL_3:
        qm = np1;
LABEL_4:
        m  = qm;
        qm = ARRAYF(q, m);
        if (qm  < vj) goto LABEL_4;
        if (qm == vj) goto LABEL_102;
        luk++;
        ARRAYF(q,  m) = vj;
        ARRAYF(q, vj) = qm;
        vj            = ARRAYF(jra, vj);
        if (vj != 0) goto LABEL_3;
//  ******  link through jru  *******************************************
LABEL_5:
        lastid         = 0;
        lasti          = 0;
        ARRAYF(ijl, k) = jlptr;
        i              = k;
LABEL_6:
        i     = ARRAYF(jru, i);
        if (i == 0) goto LABEL_10;
        qm    = np1;
        jmin  = ARRAYF(irl, i);
        jmax  = ARRAYF(ijl, i) + ARRAYF(il, i + 1) - ARRAYF(il, i) - 1;
        llong = jmax - jmin; // llong = long
        if (llong < 0) goto LABEL_6;
        jtmp  = ARRAYF(jl, jmin);
        if (jtmp != k) llong++;
        if (jtmp == k) ARRAYF(r, i) = - ARRAYF(r, i);
        if (lastid >= llong) goto LABEL_7;
        lasti  = i;
        lastid = llong;
//  ******  and merge the corresponding columns into the kth column  ****
LABEL_7:
        for (j = jmin; j <= jmax; ++j) {
            vj = ARRAYF(jl, j);
LABEL_8:
            m  = qm;
            qm = ARRAYF(q, m);
            if (qm  < vj) goto LABEL_8;
            if (qm == vj) continue;
            luk++;
            ARRAYF(q,  m) = vj;
            ARRAYF(q, vj) = qm;
            qm            = vj;
        }
        goto LABEL_6;
//  ******  lasti is the longest column merged into the kth  ************
//  ******  see if it equals the entire kth column  *********************
LABEL_10:
        qm = ARRAYF(q, np1);
        if (qm  != k) goto LABEL_105;
        if (luk == 0) goto LABEL_17;
        if (lastid != luk) goto LABEL_11;
//  ******  if so, jl can be compressed  ********************************
        irll           = ARRAYF(irl, lasti);
        ARRAYF(ijl, k) = irll + 1;
        if (ARRAYF(jl, irll) != k) ARRAYF(ijl, k) -= 1;
        goto LABEL_17;
//  ******  if not, see if kth column can overlap the previous one  *****
LABEL_11:
        if (jlmin > jlptr) goto LABEL_15;
        qm = ARRAYF(q, qm);
        for (j = jlmin; j <= jlptr; ++j) {
            if (ARRAYF(jl, j) < qm) {
                continue;
            } else if (ARRAYF(jl, j) == qm) {
                goto LABEL_13;
            } else {
                goto LABEL_15;
            }
        }
        goto LABEL_15;
LABEL_13:
        ARRAYF(ijl, k) = j;
        for (i = j; i <= jlptr; ++i) {
            if (ARRAYF(jl, i) != qm) goto LABEL_15;
            qm = ARRAYF(q, qm);
            if (qm > n) goto LABEL_17;
        }
        jlptr = j - 1;
//  ******  move column indices from q to jl, update vectors  ***********
LABEL_15:
        jlmin = jlptr + 1;
        ARRAYF(ijl, k) = jlmin;
        if (luk == 0) goto LABEL_17;
        jlptr += luk;
        if (jlptr > jlmax) goto LABEL_103;
        qm = ARRAYF(q, np1);
        for (j = jlmin; j <= jlptr; ++j) {
            qm            = ARRAYF(q, qm);
            ARRAYF(jl, j) = qm;
        }
LABEL_17:
        ARRAYF(irl,     k) = ARRAYF(ijl, k);
        ARRAYF( il, k + 1) = ARRAYF( il, k) + luk;
//
//  ******  initialize q for computing kth row of u  ********************
        ARRAYF(q, np1) = np1;
        luk  = -1;
//  ******  by filling in kth row of reordered a  ***********************
        rk   = ARRAYF(  r, k);
        jmin = ARRAYF(ira, k);
        jmax = ARRAYF( ia, rk + 1) - 1;
        if (jmin > jmax) goto LABEL_20;
        for (j = jmin; j <= jmax; ++j) {
            vj = ARRAYF(ic, ARRAYF(ja, j));
            qm = np1;
LABEL_18:
            m  = qm;
            qm = ARRAYF(q, m);
            if (qm  < vj) goto LABEL_18;
            if (qm == vj) goto LABEL_102;
            luk++;
            ARRAYF(q,  m) = vj;
            ARRAYF(q, vj) = qm;
        }
//  ******  link through jrl,  ******************************************
LABEL_20:
        lastid         = 0;
        lasti          = 0;
        ARRAYF(iju, k) = juptr;
        i              = k;
        i1             = ARRAYF(jrl, k);
LABEL_21:
        i     = i1;
        if (i == 0) goto LABEL_26;
        i1    = ARRAYF(jrl, i);
        qm    = np1;
        jmin  = ARRAYF(iru, i);
        jmax  = ARRAYF(iju, i) + ARRAYF(iu, i + 1) - ARRAYF(iu, i) - 1;
        llong = jmax - jmin;
        if (llong < 0) goto LABEL_21;
        jtmp = ARRAYF(ju, jmin);
        if (jtmp == k) goto LABEL_22;
//  ******  update irl and jrl, *****************************************
        llong          += 1;
        cend            = ARRAYF(ijl, i) + ARRAYF(il, i + 1) - ARRAYF(il, i);
        ARRAYF(irl, i) += 1;
        if (ARRAYF(irl, i) >= cend) goto LABEL_22;
        j              = ARRAYF(jl, ARRAYF(irl, i));
        ARRAYF(jrl, i) = ARRAYF(jrl, j);
        ARRAYF(jrl, j) = i;
LABEL_22:
        if (lastid >= llong) goto LABEL_23;
        lasti  = i;
        lastid = llong;
//  ******  and merge the corresponding rows into the kth row  **********
LABEL_23:
        for (j = jmin; j <= jmax; ++j) {
            vj = ARRAYF(ju, j);
LABEL_24:
            m  = qm;
            qm = ARRAYF(q, m);
            if (qm  < vj) goto LABEL_24;
            if (qm == vj) continue;;
            luk++;
            ARRAYF(q,  m) = vj;
            ARRAYF(q, vj) = qm;
            qm            = vj;
        }
        goto LABEL_21;
//  ******  update jrl(k) and irl(k)  ***********************************
LABEL_26:
        if (ARRAYF(il, k + 1) <= ARRAYF(il, k)) goto LABEL_27;
        j              = ARRAYF(jl, ARRAYF(irl, k));
        ARRAYF(jrl, k) = ARRAYF(jrl, j);
        ARRAYF(jrl, j) = k;
//  ******  lasti is the longest row merged into the kth  ***************
//  ******  see if it equals the entire kth row  ************************
LABEL_27:
        qm = ARRAYF(q, np1);
        if (qm  != k) goto LABEL_105;
        if (luk == 0) goto LABEL_34;
        if (lastid != luk) goto LABEL_28;
//  ******  if so, ju can be compressed  ********************************
        irul           = ARRAYF(iru, lasti);
        ARRAYF(iju, k) = irul + 1;
        if (ARRAYF(ju, irul) != k) ARRAYF(iju, k) -= 1;
        goto LABEL_34;
//  ******  if not, see if kth row can overlap the previous one  ********
LABEL_28:
        if (jumin > juptr) goto LABEL_32;
        qm = ARRAYF(q, qm);
        for (j = jumin; j <= juptr; ++j) {
            if (ARRAYF(ju, j) < qm) {
                continue;
            } else if (ARRAYF(ju, j) == qm) {
                goto LABEL_30;
            } else {
                goto LABEL_32;
            }
        }
        goto LABEL_32;
LABEL_30:
        ARRAYF(iju, k) = j;
        for (i = j; i <= juptr; ++i) {
            if (ARRAYF(ju, i) != qm) goto LABEL_32;
            qm = ARRAYF(q, qm);
            if (qm > n) goto LABEL_34;
        }
        juptr = j - 1;
//  ******  move row indices from q to ju, update vectors  **************
LABEL_32:
        jumin          = juptr + 1;
        ARRAYF(iju, k) = jumin;
        if (luk == 0) goto LABEL_34;
        juptr         += luk;
        if (juptr > jumax) goto LABEL_106;
        qm = ARRAYF(q, np1);
        for (j = jumin; j <= juptr; ++j) {
            qm            = ARRAYF(q, qm);
            ARRAYF(ju, j) = qm;
        }
LABEL_34:
        ARRAYF(iru,     k) = ARRAYF(iju, k);
        ARRAYF( iu, k + 1) = ARRAYF( iu, k) + luk;
//
//  ******  update iru, jru  ********************************************
        i = k;
LABEL_35:
        i1 = ARRAYF(jru, i);
        if (ARRAYF(r, i) < 0) goto LABEL_36;
        rend = ARRAYF(iju, i) + ARRAYF(iu, i + 1) - ARRAYF(iu, i);
        if (ARRAYF(iru, i) >= rend) goto LABEL_37;
        j              = ARRAYF(ju, ARRAYF(iru, i));
        ARRAYF(jru, i) = ARRAYF(jru, j);
        ARRAYF(jru, j) = i;
        goto LABEL_37;
LABEL_36:
        ARRAYF(r, i) = - ARRAYF(r, i);
LABEL_37:
        i = i1;
        if (i == 0) goto LABEL_38;
        ARRAYF(iru, i) += 1;
        goto LABEL_35;
//
//  ******  update ira, jra, irac  **************************************
LABEL_38:
        i = ARRAYF(irac, k);
        if (i == 0) continue;
LABEL_39:
        i1 = ARRAYF(jra, i);
        ARRAYF(ira, i) += 1;
        if (ARRAYF(ira, i) >= ARRAYF(ia, ARRAYF(r, i) + 1)) goto LABEL_40;
        irai   = ARRAYF(ira, i);
        jairai = ARRAYF( ic, ARRAYF(ja, irai));
        if (jairai > i) goto LABEL_40;
        ARRAYF( jra,      i) = ARRAYF(irac, jairai);
        ARRAYF(irac, jairai) = i;
LABEL_40:
        i = i1;
        if (i != 0) goto LABEL_39;
    }
//
    ARRAYF(ijl, n) = jlptr;
    ARRAYF(iju, n) = juptr;
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

void Odepack::NNFC(const int n, const int *r, const int *c, const int *ic, const int *ia, 
    const int *ja, const double *a, double *z, const double *b, const int lmax, 
    const int *il, const int *jl, const int *ijl, double *l, double *d, const int umax,
    const int *iu, const int *ju, const int *iju, double *u, double *row, double *tmp, 
    int *irl, int *jrl, int &flag)
{
    int i, i1, i2, k, j, jmin, jmax, rk, mu, ijlb;
    double lki, sum, dk;
//
//  ******  initialize pointers and test storage  ***********************
    if (ARRAYF(il, n + 1) - 1 > lmax) goto LABEL_104;
    if (ARRAYF(iu, n + 1) - 1 > umax) goto LABEL_107;
    for (k = 1; k <= n; ++k) {
        ARRAYF(irl, k) = ARRAYF(il, k);
        ARRAYF(jrl, k) = 0;
    }
//
//  ******  for each row  ***********************************************
    for (k = 1; k <= n; ++k) {
//  ******  reverse jrl and zero row where kth row of l will fill in  ***
        ARRAYF(row, k) = 0;
        i1             = 0;
        if (ARRAYF(jrl, k) == 0) goto LABEL_3;
        i              = ARRAYF(jrl, k);
LABEL_2:
        i2             = ARRAYF(jrl, i);
        ARRAYF(jrl, i) = i1;
        i1             = i;
        ARRAYF(row, i) = 0;
        i              = i2;
        if (i != 0) goto LABEL_2;
//  ******  set row to zero where u will fill in  ***********************
LABEL_3:
        jmin = ARRAYF(iju, k);
        jmax = jmin + ARRAYF(iu, k + 1) - ARRAYF(iu, k) - 1;
        if (jmin > jmax) goto LABEL_5;
        for (j = jmin; j <= jmax; ++j) {
            ARRAYF(row, ARRAYF(ju, j)) = 0;
        }
//  ******  place kth row of a in row  **********************************
LABEL_5:
        rk   = ARRAYF( r, k);
        jmin = ARRAYF(ia, rk);
        jmax = ARRAYF(ia, rk + 1) - 1;
        for (j = jmin; j <= jmax; ++j) {
            ARRAYF(row, ARRAYF(ic, ARRAYF(ja, j))) = ARRAYF(a, j);
        }
//  ******  initialize sum, and link through jrl  ***********************
        sum = ARRAYF(b, rk);
        i   = i1;
        if (i == 0) goto LABEL_10;
//  ******  assign the kth row of l and adjust row, sum  ****************
LABEL_7:
        lki = - ARRAYF(row, i);
//  ******  if l is not required, then comment out the following line  **
        ARRAYF(l, ARRAYF(irl, i)) = -lki;
        sum  = sum + lki * ARRAYF(tmp, i);
        jmin = ARRAYF(iu, i);
        jmax = ARRAYF(iu, i + 1) - 1;
        if (jmin > jmax) goto LABEL_9;
        mu   = ARRAYF(iju, i) - jmin;
        for (j = jmin; j <= jmax; ++j) {
            ARRAYF(row, ARRAYF(ju, mu + j)) += lki * ARRAYF(u, j);
        }
LABEL_9:
        i = ARRAYF(jrl, i);
        if (i != 0) goto LABEL_7;
//
//  ******  assign kth row of u and diagonal d, set tmp(k)  *************
LABEL_10:
        if (ARRAYF(row, k) == 0.0) goto LABEL_108;
        dk             = 1.0 / ARRAYF(row, k);
        ARRAYF(  d, k) = dk;
        ARRAYF(tmp, k) = sum * dk;
        if (k == n) continue;
        jmin = ARRAYF(iu, k);
        jmax = ARRAYF(iu, k + 1) - 1;
        if (jmin > jmax) goto LABEL_12;
        mu   = ARRAYF(iju, k) - jmin;
        for (j = jmin; j <= jmax; ++j) {
            ARRAYF(u, j) = ARRAYF(row, ARRAYF(ju, mu + j)) * dk;
        }
LABEL_12:
//
//  ******  update irl and jrl, keeping jrl in decreasing order  ********
        i = i1;
        if (i == 0) goto LABEL_18;
LABEL_14:
        ARRAYF(irl, i) += 1;
        i1   = ARRAYF(jrl, i);
        if (ARRAYF(irl, i) >= ARRAYF(il, i + 1)) goto LABEL_17;
        ijlb = ARRAYF(irl, i) - ARRAYF(il, i) + ARRAYF(ijl, i);
        j    = ARRAYF(jl, ijlb);
LABEL_15:
        if (i > ARRAYF(jrl, j)) goto LABEL_16;
        j = ARRAYF(jrl, j);
        goto LABEL_15;
LABEL_16:
        ARRAYF(jrl, i) = ARRAYF(jrl, j);
        ARRAYF(jrl, j) = i;
LABEL_17:
        i = i1;
        if (i != 0) goto LABEL_14;
LABEL_18:
        if (ARRAYF(irl, k) >= ARRAYF(il, k + 1)) continue;
        j              = ARRAYF(jl, ARRAYF(ijl, k));
        ARRAYF(jrl, k) = ARRAYF(jrl, j);
        ARRAYF(jrl, j) = k;
    }
//
//  ******  solve  ux = tmp  by back substitution  **********************
    k = n;
    for (i = 1; i <= n; ++i) {
        sum  = ARRAYF(tmp, k);
        jmin = ARRAYF( iu, k);
        jmax = ARRAYF( iu, k + 1) - 1;
        if (jmin > jmax) goto LABEL_21;
        mu = ARRAYF(iju, k) - jmin;
        for (j = jmin; j <= jmax; ++j) {
            sum = sum - ARRAYF(u, j) * ARRAYF(tmp, ARRAYF(ju, mu + j));
        }
LABEL_21:
        ARRAYF(tmp, k)          = sum;
        ARRAYF(z, ARRAYF(c, k)) = sum;
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

void Odepack::NNSC(const int n, const int *r, const int *c, const int *il, const int *jl, 
    const int *ijl, const double *l, const double *d, const int *iu, const int *ju, 
    const int *iju, const double *u, double *z, const double *b, double *tmp)
{
    int i, k, j, jmin, jmax, ml, mu;
    double tmpk, sum;
//
//  ******  set tmp to reordered b  *************************************
    for (k = 1; k <= n; ++k) {
        ARRAYF(tmp, k) = ARRAYF(b, ARRAYF(r, k));
    }
//  ******  solve  ly = b  by forward substitution  *********************
    for (k = 1; k <= n; ++k) {
        jmin           = ARRAYF(il, k);
        jmax           = ARRAYF(il, k + 1) - 1;
        tmpk           = - ARRAYF(d, k) * ARRAYF(tmp, k);
        ARRAYF(tmp, k) = - tmpk;
        if (jmin > jmax) continue;
        ml = ARRAYF(ijl, k) - jmin;
        for (j = jmin; j <= jmax; ++j) {
            ARRAYF(tmp, ARRAYF(jl, ml + j)) += tmpk * ARRAYF(l, j);
        }
    }
//  ******  solve  ux = y  by back substitution  ************************
    k = n;
    for (i = 1; i <= n; ++i) {
        sum  = - ARRAYF(tmp, k);
        jmin = ARRAYF(iu, k);
        jmax = ARRAYF(iu, k + 1) - 1;
        if (jmin > jmax) goto LABEL_5;
        mu   = ARRAYF(iju, k) - jmin;
        for (j = jmin; j <= jmax; ++j) {
            sum += ARRAYF(u, j) * ARRAYF(tmp, ARRAYF(ju, mu + j));
        }
LABEL_5:
        ARRAYF(tmp, k)          = - sum;
        ARRAYF(z, ARRAYF(c, k)) = - sum;
        k--;
    }
    return;
}

void Odepack::NNTC(const int n, const int *r, const int *c, const int *il, const int *jl, 
    const int *ijl, const double *l, const double *d, const int *iu, const int *ju, 
    const int *iju, const double *u, double *z, const double *b, double *tmp)
{
    int i, j, jmin, jmax, k, mu, ml;
    double tmpk, sum;
//
//  ******  set tmp to reordered b  *************************************
    for (k = 1; k <= n; ++k) {
        ARRAYF(tmp, k) = ARRAYF(b, ARRAYF(c, k));
    }
//  ******  solve  ut y = b  by forward substitution  *******************
    for (k = 1; k <= n; ++k) {
        jmin = ARRAYF(iu, k);
        jmax = ARRAYF(iu, k + 1) - 1;
        tmpk = - ARRAYF(tmp, k);
        if (jmin > jmax) continue;
        mu   = ARRAYF(iju, k) - jmin;
        for (j = jmin; j <= jmax; ++j) {
            ARRAYF(tmp, ARRAYF(ju, mu + j)) += tmpk * ARRAYF(u, j);
        }
    }
//  ******  solve  lt x = y  by back substitution  **********************
    k = n;
    for (i = 1; i <= n; ++i) {
        sum  = - ARRAYF(tmp, k);
        jmin = ARRAYF(il, k);
        jmax = ARRAYF(il, k + 1) - 1;
        if (jmin > jmax) goto LABEL_5;
        ml   = ARRAYF(ijl, k) - jmin;
        for (j = jmin; j <= jmax; ++j) {
            sum += ARRAYF(l, j) * ARRAYF(tmp, ARRAYF(jl, ml + j));
        }
LABEL_5:
        ARRAYF(tmp, k)          = - sum * ARRAYF(d, k);
        ARRAYF(z, ARRAYF(r, k)) = ARRAYF(tmp, k);
        k--;
    }
    return;
}

void Odepack::DSTODA(const int neq, double *y, double *yh, const int nyh, double *yh1, 
        double *ewt, double *savf, double *acor, double *wm, int *iwm, 
        ODEPACK_FUNCTION f, ODEPACK_JACOBIAN1 jac, FUNC_PJAC<ODEPACK_JACOBIAN1> pjac, FUNC_SLVS slvs, 
        void *user_data)
{
#ifndef YH
#define YH(i, j) MATF(yh, nyh, i, j)
#endif
#ifndef ELCO
#define ELCO(i, j) MATF(dls1_.elco, 13, i, j)
#endif
#ifndef TESCO
#define TESCO(i, j) MATF(dls1_.tesco, 3, i, j)
#endif
//
    int i, j, i1, jb, iredo, iret, m, ncf, newq, lm1, lm1p1, lm2, lm2p1, nqm1, nqm2;
    double dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup, r, rh, rhdn, rhsm, rhup, told, dmnorm, rate, rm,
       rh1, rh2, rh1it, exm2, dm2, exm1, dm1, alpha;
    double pdh, pnorm;
    const double sm1[12] = {
        0.5, 0.575, 0.55, 0.45, 0.35, 0.25, 0.20, 0.15, 0.10, 0.075, 0.050, 0.025
    };
//
    dls1_.kflag = 0;
    told   = dls1_.tn;
    ncf    = 0;
    dls1_.ierpj = 0;
    dls1_.iersl = 0;
    dls1_.jcur  = 0;
    dls1_.icf   = 0;
    delp   = 0.0;
    if (dls1_.jstart > 0)   goto LABEL_200;
    if (dls1_.jstart == -1) goto LABEL_100;
    if (dls1_.jstart == -2) goto LABEL_160;
//-----------------------------------------------------------------------
// On the first call, the order is set to 1, and other variables are
// initialized.  RMAX is the maximum ratio by which H can be increased
// in a single step.  It is initially 1.E4 to compensate for the small
// initial H, but then is normally equal to 10.  If a failure
// occurs (in corrector convergence or error test), RMAX is set at 2
// for the next increase.
// DCFODE is called to get the needed coefficients for both methods.
//-----------------------------------------------------------------------
    dls1_.lmax  = dls1_.maxord + 1;
    dls1_.nq    = 1;
    dls1_.l     = 2;
    dls1_.ialth = 2;
    dls1_.rmax  = 10000.0;
    dls1_.rc    = 0.0;
    dls1_.el0   = 1.0;
    dls1_.crate = 0.7;
    dls1_.hold  = dls1_.h;
    dls1_.nslp  = 0;
    dls1_.ipup  = dls1_.miter;
    iret   = 3;
// Initialize switching parameters.  METH = 1 is assumed initially. -----
    dlsa_.icount = 20;
    dlsa_.irflag = 0;
    dlsa_.pdest  = 0.0;
    dlsa_.pdlast = 0.0;
    dlsa_.ratio  = 5.0;
    DCFODE(2, dls1_.elco, dls1_.tesco);
    for (i = 1; i <= 5; ++i) {
        ARRAYF(dlsa_.cm2, i) = TESCO(2, i) * ELCO(i+1, i);
    }
    DCFODE(1, dls1_.elco, dls1_.tesco);
    for (i = 1; i <= 12; ++i) {
        ARRAYF(dlsa_.cm1, i) = TESCO(2, i) * ELCO(i+1, i);
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
    dls1_.ipup = dls1_.miter;
    dls1_.lmax = dls1_.maxord + 1;
    if (dls1_.ialth == 1) dls1_.ialth = 2;
    if (dls1_.meth == dlsa_.mused) goto LABEL_160;
    DCFODE(dls1_.meth, dls1_.elco, dls1_.tesco);
    dls1_.ialth = dls1_.l;
    iret = 1;
//-----------------------------------------------------------------------
// The el vector and related constants are reset
// whenever the order NQ is changed, or at the start of the problem.
//-----------------------------------------------------------------------
LABEL_150:
    for (i = 1; i <= dls1_.l; ++i) {
        ARRAYF(dls1_.el, i) = ELCO(i, dls1_.nq);
    }
    dls1_.nqnyh = dls1_.nq * dls1_.nyh;
    dls1_.rc    = dls1_.rc * ARRAYF(dls1_.el, 1) / dls1_.el0;
    dls1_.el0   = ARRAYF(dls1_.el, 1);
    dls1_.coint = 0.5 / static_cast<double>(dls1_.nq + 2);
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
    if (dls1_.h == dls1_.hold) goto LABEL_200;
    rh = dls1_.h / dls1_.hold;
    dls1_.h = dls1_.hold;
    iredo = 3;
    goto LABEL_175;
LABEL_170:
    rh = std::max(rh, dls1_.hmin / std::abs(dls1_.h));
LABEL_175:
    rh = std::min(rh, dls1_.rmax);
    rh = rh / std::max(1.0, std::abs(dls1_.h) * dls1_.hmxi * rh);
//-----------------------------------------------------------------------
// If METH = 1, also restrict the new step size by the stability region.
// If this reduces H, set IRFLAG to 1 so that if there are roundoff
// problems later, we can assume that is the cause of the trouble.
//-----------------------------------------------------------------------
    if (dls1_.meth == 2) goto LABEL_178;
    dlsa_.irflag = 0;
    pdh = std::max(std::abs(dls1_.h) * dlsa_.pdlast, 0.000001);
    if (rh * pdh * 1.00001 < ARRAYF(sm1, dls1_.nq)) goto LABEL_178;
    rh = ARRAYF(sm1, dls1_.nq) / pdh;
    dlsa_.irflag = 1;
LABEL_178:
    r = 1.0;
    for (j = 2; j <= dls1_.l; ++j) {
        r *= rh;
        for (i = 1; i <= dls1_.n; ++i) {
            YH(i, j) *= r;
        }
    }
    dls1_.h    *= rh;
    dls1_.rc   *= rh;
    dls1_.ialth = dls1_.l;
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
    if (std::abs(dls1_.rc - 1.0) > dls1_.ccmax) dls1_.ipup = dls1_.miter; 
    if (dls1_.nst >= dls1_.nslp + dls1_.msbp) dls1_.ipup = dls1_.miter;
    dls1_.tn += dls1_.h;
    i1 = dls1_.nqnyh + 1;
    for (jb = 1; jb <= dls1_.nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= dls1_.nqnyh; ++i) {
            ARRAYF(yh1, i) += ARRAYF(yh1, i + nyh);
        }
    }
    pnorm = DMNORM(dls1_.n, yh1, ewt);
//-----------------------------------------------------------------------
// Up to MAXCOR corrector iterations are taken.  A convergence test is
// made on the RMS-norm of each correction, weighted by the error
// weight vector EWT.  The sum of the corrections is accumulated in the
// vector ACOR(i).  The YH array is not altered in the corrector loop.
//-----------------------------------------------------------------------
 LABEL_220:
    m    = 0;
    rate = 0.0;
    del  = 0.0;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = YH(i, 1);
    }
    f(neq, dls1_.tn, y, savf, user_data);
    dls1_.nfe++;
    if (dls1_.ipup <= 0) goto LABEL_250;
//-----------------------------------------------------------------------
// If indicated, the matrix P = I - H*EL(1)*J is reevaluated and
// preprocessed before starting the corrector iteration.  IPUP is set
// to 0 as an indicator that this has been done.
//-----------------------------------------------------------------------
    (this->*pjac)(neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac, user_data);
    dls1_.ipup  = 0;
    dls1_.rc    = 1.0;
    dls1_.nslp  = dls1_.nst;
    dls1_.crate = 0.7;
    if (dls1_.ierpj != 0) goto LABEL_430;
LABEL_250:
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(acor, i) = 0.0;
    }
LABEL_270:
    if (dls1_.miter != 0) goto LABEL_350;
//-----------------------------------------------------------------------
// In the case of functional iteration, update Y directly from
// the result of the last function evaluation.
//-----------------------------------------------------------------------
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(savf, i) = dls1_.h * ARRAYF(savf, i) - YH(i, 2);
        ARRAYF(   y, i) = ARRAYF(savf, i) - ARRAYF(acor, i);
    }
    del = DMNORM(dls1_.n, y, ewt);
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(   y, i) = YH(i, 1) + ARRAYF(dls1_.el, 1) * ARRAYF(savf, i);
        ARRAYF(acor, i) = ARRAYF(savf, i);
    }
    goto LABEL_400;
//-----------------------------------------------------------------------
// In the case of the chord method, compute the corrector error,
// and solve the linear system with that as right-hand side and
// P as coefficient matrix.
//-----------------------------------------------------------------------
LABEL_350:
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = dls1_.h * ARRAYF(savf, i) - (YH(i, 2) + ARRAYF(acor, i));
    }
    (this->*slvs)(wm, iwm, y, savf);
    if (dls1_.iersl < 0) goto LABEL_430;
    if (dls1_.iersl > 0) goto LABEL_410;
    del = DMNORM(dls1_.n, y, ewt);
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(acor, i) += ARRAYF(y, i);
        ARRAYF(   y, i) = YH(i, 1) + ARRAYF(dls1_.el, 1) * ARRAYF(acor, i);
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
    if (del <= 100.0 * pnorm * dls1_.uround) goto LABEL_450;
    if (m == 0 && dls1_.meth == 1) goto LABEL_405;
    if (m == 0) goto LABEL_402;
    rm = 1024.0;
    if (del <= 1024.0 * delp) rm = del / delp;
    rate = std::max(rate, rm);
    dls1_.crate = std::max(0.2 * dls1_.crate, rm);
LABEL_402:
    dcon = del * std::min(1.0, 1.5 * dls1_.crate) / (TESCO(2, dls1_.nq) * dls1_.coint);
    if (dcon > 1.0) goto LABEL_405;
    dlsa_.pdest = std::max(dlsa_.pdest, rate / std::abs(dls1_.h * ARRAYF(dls1_.el, 1)));
    if (dlsa_.pdest != 0.0) dlsa_.pdlast = dlsa_.pdest;
    goto LABEL_450;
LABEL_405:
    m++;
    if (m == dls1_.maxcor) goto LABEL_410;
    if (m >= 2 && del > 2.0 * delp) goto LABEL_410;
    delp = del;
    f(neq, dls1_.tn, y, savf, user_data);
    dls1_.nfe++;
    goto LABEL_270;
//-----------------------------------------------------------------------
// The corrector iteration failed to converge.
// If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
// the next try.  Otherwise the YH array is retracted to its values
// before prediction, and H is reduced, if possible.  If H cannot be
// reduced or MXNCF failures have occurred, exit with KFLAG = -2.
//-----------------------------------------------------------------------
LABEL_410:
    if (dls1_.miter == 0 || dls1_.jcur == 1) goto LABEL_430;
    dls1_.icf  = 1;
    dls1_.ipup = dls1_.miter;
    goto LABEL_220;
LABEL_430:
    dls1_.icf = 2;
    ncf++;
    dls1_.rmax = 2.0;
    dls1_.tn = told;
    i1 = dls1_.nqnyh + 1;
    for (jb = 1; jb <= dls1_.nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= dls1_.nqnyh; ++i) {
            ARRAYF(yh1, i) -= ARRAYF(yh1, i + nyh);
        }
    }
    if (dls1_.ierpj < 0 || dls1_.iersl < 0) goto LABEL_680;
    if (std::abs(dls1_.h) <= dls1_.hmin * 1.00001) goto LABEL_670;
    if (ncf == dls1_.mxncf) goto LABEL_670;
    rh    = 0.25;
    dls1_.ipup = dls1_.miter;
    iredo = 1;
    goto LABEL_170;
//-----------------------------------------------------------------------
// The corrector has converged.  JCUR is set to 0
// to signal that the Jacobian involved may need updating later.
// The local error test is made and control passes to statement 500
// if it fails.
//-----------------------------------------------------------------------
LABEL_450:
    dls1_.jcur = 0;
    if (m == 0) dsm = del / TESCO(2, dls1_.nq);
    if (m > 0)  dsm = DMNORM(dls1_.n, acor, ewt) / TESCO(2, dls1_.nq);
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
    dls1_.kflag = 0;
    iredo  = 0;
    dls1_.nst++;
    dls1_.hu    = dls1_.h;
    dls1_.nqu   = dls1_.nq;
    dlsa_.mused = dls1_.meth;
    for (j = 1; j <= dls1_.l; ++j) {
        for (i = 1; i <= dls1_.n; ++i) {
            YH(i, j) += ARRAYF(dls1_.el, j) * ARRAYF(acor, i);
        }
    }
    dlsa_.icount--;
    if (dlsa_.icount >= 0) goto LABEL_488;
    if (dls1_.meth == 2) goto LABEL_480;
//-----------------------------------------------------------------------
// We are currently using an Adams method.  Consider switching to BDF.
// If the current order is greater than 5, assume the problem is
// not stiff, and skip this section.
// If the Lipschitz constant and error estimate are not polluted
// by roundoff, go to 470 and perform the usual test.
// Otherwise, switch to the BDF methods if the last step was
// restricted to insure stability (irflag = 1), and stay with Adams
// method if not.  When switching to BDF with polluted error estimates,
// in the absence of other information, double the step size.
//
// When the estimates are OK, we make the usual test by computing
// the step size we could have (ideally) used on this step,
// with the current (Adams) method, and also that for the BDF.
// If NQ .gt. MXORDS, we consider changing to order MXORDS on switching.
// Compare the two step sizes to decide whether to switch.
// The step size advantage must be at least RATIO = 5 to switch.
//-----------------------------------------------------------------------
    if (dls1_.nq > 5) goto LABEL_488;
    if (dsm > 100.0 * pnorm * dls1_.uround && dlsa_.pdest != 0.0) goto LABEL_470;
    if (dlsa_.irflag == 0) goto LABEL_488;
    rh2 = 2.0;
    nqm2 = std::min(dls1_.nq, dlsa_.mxords);
    goto LABEL_478;
LABEL_470:
    exsm  = 1.0 / static_cast<double>(dls1_.l);
    rh1   = 1.0 / (1.2 * std::pow(dsm, exsm) + 0.0000012);
    rh1it = 2.0 * rh1;
    pdh   = dlsa_.pdlast * std::abs(dls1_.h);
    if (pdh * rh1 > 0.00001) rh1it = ARRAYF(sm1, dls1_.nq) / pdh;
    rh1   = std::min(rh1, rh1it);
    if (dls1_.nq <= dlsa_.mxords) goto LABEL_474;
    nqm2  = dlsa_.mxords;
    lm2   = dlsa_.mxords + 1;
    exm2  = 1.0 / static_cast<double>(lm2);
    lm2p1 = lm2 + 1;
    dm2   = DMNORM(dls1_.n, &YH(1, lm2p1), ewt) / ARRAYF(dlsa_.cm2, dlsa_.mxords);
    rh2   = 1.0 / (1.2 * std::pow(dm2, exm2) + 0.0000012);
    goto LABEL_476;
LABEL_474:
    dm2  = dsm * (ARRAYF(dlsa_.cm1, dls1_.nq) / ARRAYF(dlsa_.cm2, dls1_.nq));
    rh2  = 1.0 / (1.2 * std::pow(dm2, exsm) + 0.0000012);
    nqm2 = dls1_.nq;
LABEL_476:
    if (rh2 < dlsa_.ratio * rh1) goto LABEL_488;
// THE SWITCH TEST PASSED.  RESET RELEVANT QUANTITIES FOR BDF. ----------
LABEL_478:
    rh      = rh2;
    dlsa_.icount = 20;
    dls1_.meth   = 2;
    dls1_.miter  = dlsa_.jtyp;
    dlsa_.pdlast = 0.0;
    dls1_.nq     = nqm2;
    dls1_.l      = dls1_.nq + 1;
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
    exsm = 1.0 / static_cast<double>(dls1_.l);
    if (dlsa_.mxordn >= dls1_.nq) goto LABEL_484;
    nqm1  = dlsa_.mxordn;
    lm1   = dlsa_.mxordn + 1;
    exm1  = 1.0 / static_cast<double>(lm1);
    lm1p1 = lm1 + 1;
    dm1   = DMNORM(dls1_.n, &YH(1, lm1p1), ewt) / ARRAYF(dlsa_.cm1, dlsa_.mxordn);
    rh1   = 1.0 / (1.2 * std::pow(dm1, exm1) + 0.0000012);
    goto LABEL_486;
LABEL_484:
    dm1 = dsm * (ARRAYF(dlsa_.cm2, dls1_.nq) / ARRAYF(dlsa_.cm1, dls1_.nq));
    rh1 = 1.0 / (1.2 * std::pow(dm1, exsm) + 0.0000012);
    nqm1 = dls1_.nq;
    exm1 = exsm;
LABEL_486:
    rh1it = 2.0 * rh1;
    pdh = dlsa_.pdnorm * std::abs(dls1_.h);
    if (pdh * rh1 > 0.00001) rh1it = ARRAYF(sm1, nqm1) / pdh;
    rh1 = std::min(rh1, rh1it);
    rh2 = 1.0 / (1.2 * std::pow(dsm, exsm) + 0.0000012);
    if (rh1 * dlsa_.ratio < 5.0 * rh2) goto LABEL_488;
    alpha = std::max(0.001, rh1);
    dm1 = std::pow(alpha, exm1) * dm1;
    if (dm1 <= 1000.0 * dls1_.uround * pnorm) goto LABEL_488;
// The switch test passed.  Reset relevant quantities for Adams. --------
    rh      = rh1;
    dlsa_.icount = 20;
    dls1_.meth   = 1;
    dls1_.miter  = 0;
    dlsa_.pdlast = 0.0;
    dls1_.nq     = nqm1;
    dls1_.l      = dls1_.nq + 1;
    goto LABEL_170;
// No method switch is being made.  Do the usual step/order selection. --
LABEL_488:
    dls1_.ialth--;
    if (dls1_.ialth == 0) goto LABEL_520;
    if (dls1_.ialth > 1)  goto LABEL_700;
    if (dls1_.l == dls1_.lmax) goto LABEL_700;
    for (i = 1; i <= dls1_.n; ++i) {
        YH(i, dls1_.lmax) = ARRAYF(acor, i);
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
    dls1_.kflag--;
    dls1_.tn = told;
    i1 = dls1_.nqnyh + 1;
    for (jb = 1; jb <= dls1_.nq; ++jb) {
        i1 -= nyh;
        for (i = i1; i <= dls1_.nqnyh; ++i) {
            ARRAYF(yh1, i) -= ARRAYF(yh1, i+nyh);
        }
    }
    dls1_.rmax = 2.0;
    if (std::abs(dls1_.h) <= dls1_.hmin * 1.00001) goto LABEL_660;
    if (dls1_.kflag <= -3) goto LABEL_640;
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
    if (dls1_.l == dls1_.lmax) goto LABEL_540;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(savf, i) = ARRAYF(acor, i) - YH(i, dls1_.lmax);
    }
    dup  = DMNORM(dls1_.n, savf, ewt) / TESCO(3, dls1_.nq);
    exup = 1.0 / static_cast<double>(dls1_.l + 1);
    rhup = 1.0 / (1.4 * std::pow(dup, exup) + 0.0000014);
LABEL_540:
    exsm = 1.0 / static_cast<double>(dls1_.l);
    rhsm = 1.0 / (1.2 * std::pow(dsm, exsm) + 0.0000012);
    rhdn = 0.0;
    if (dls1_.nq == 1) goto LABEL_550;
    ddn = DMNORM(dls1_.n, &YH(1, dls1_.l), ewt) / TESCO(1, dls1_.nq);
    exdn = 1.0 / static_cast<double>(dls1_.nq);
    rhdn = 1.0 / (1.3 * std::pow(ddn, exdn) + 0.0000013);
// If METH = 1, limit RH according to the stability region also. --------
LABEL_550:
    if (dls1_.meth == 2) goto LABEL_560;
    pdh = std::max(std::abs(dls1_.h) * dlsa_.pdlast, 0.000001);
    if (dls1_.l < dls1_.lmax) rhup = std::min(rhup, ARRAYF(sm1, dls1_.l) / pdh);
    rhsm = std::min(rhsm, ARRAYF(sm1, dls1_.nq) / pdh);
    if (dls1_.nq > 1) rhdn = std::min(rhdn, ARRAYF(sm1, dls1_.nq-1) / pdh);
    dlsa_.pdest = 0.0;
LABEL_560:
    if (rhsm >= rhup) goto LABEL_570;
    if (rhup > rhdn)  goto LABEL_590;
    goto LABEL_580;
LABEL_570:
    if (rhsm < rhdn) goto LABEL_580;
    newq = dls1_.nq;
    rh = rhsm;
    goto LABEL_620;
LABEL_580:
    newq = dls1_.nq - 1;
    rh = rhdn;
    if (dls1_.kflag < 0 && rh > 1.0) rh = 1.0;
    goto LABEL_620;
LABEL_590:
    newq = dls1_.l;
    rh = rhup;
    if (rh < 1.1) goto LABEL_610;
    r = ARRAYF(dls1_.el, dls1_.l) / static_cast<double>(dls1_.l);
    for (i = 1; i <= dls1_.n; ++i) {
        YH(i, newq+1) = ARRAYF(acor, i) * r;
    }
    goto LABEL_630;
LABEL_610:
    dls1_.ialth = 3;
    goto LABEL_700;
// If METH = 1 and H is restricted by stability, bypass 10 percent test.
LABEL_620:
    if (dls1_.meth == 2) goto LABEL_622;
    if (rh * pdh * 1.00001 >= ARRAYF(sm1, newq)) goto LABEL_625;
LABEL_622:
    if (dls1_.kflag == 0 && rh < 1.1) goto LABEL_610;
LABEL_625:
    if (dls1_.kflag <= -2) rh = std::min(rh, 0.2);
//-----------------------------------------------------------------------
// If there is a change of order, reset NQ, L, and the coefficients.
// In any case H is reset according to RH and the YH array is rescaled.
// Then exit from 690 if the step was OK, or redo the step otherwise.
//-----------------------------------------------------------------------
    if (newq == dls1_.nq) goto LABEL_170;
LABEL_630:
    dls1_.nq = newq;
    dls1_.l = dls1_.nq + 1;
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
    if (dls1_.kflag == -10) goto LABEL_660;
    rh = 0.1;
    rh = std::max(dls1_.hmin / std::abs(dls1_.h), rh);
    dls1_.h *= rh;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = YH(i, 1);
    }
    f(neq, dls1_.tn, y, savf, user_data);
    dls1_.nfe++;
    for (i = 1; i <= dls1_.n; ++i) {
        YH(i, 2) = dls1_.h * ARRAYF(savf, i);
    }
    dls1_.ipup = dls1_.miter;
    dls1_.ialth = 5;
    if (dls1_.nq == 1) goto LABEL_200;
    dls1_.nq  = 1;
    dls1_.l   = 2;
    iret = 3;
    goto LABEL_150;
//-----------------------------------------------------------------------
// All returns are made through this section.  H is saved in HOLD
// to allow the caller to change H on the next step.
//-----------------------------------------------------------------------
LABEL_660:
    dls1_.kflag = -1;
    goto LABEL_720;
LABEL_670:
    dls1_.kflag = -2;
    goto LABEL_720;
LABEL_680:
    dls1_.kflag = -3;
    goto LABEL_720;
LABEL_690:
    dls1_.rmax = 10.0;
LABEL_700:
    r = 1.0 / TESCO(2, dls1_.nqu);
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(acor, i) *= r;
    }
LABEL_720:
    dls1_.hold = dls1_.h;
    dls1_.jstart = 1;
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

void Odepack::DPRJA(const int neq,  double *y, double *yh, const int nyh, const double *ewt,  
        double *ftem, double *savf, double *wm, int *iwm, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN1 jac, 
        void *user_data)
{
#ifndef YH
#define YH(i, j) MATF(yh, nyh, i, j)
#endif
//
    int i, i1, i2, ier, ii, j, j1, jj, lenp, mba, mband, meb1, meband,
        ml, ml3, mu, np1;
    double con, fac, hl0, r, r0, srur, yi, yj, yjj;
//
    dls1_.nje++;
    dls1_.ierpj = 0;
    dls1_.jcur  = 1;
    hl0    = dls1_.h * dls1_.el0;
    if (dls1_.miter == 1) {
        goto LABEL_100;
    } else if (dls1_.miter == 2) {
        goto LABEL_200;
    } else if (dls1_.miter == 3) {
        goto LABEL_300;
    } else if (dls1_.miter == 4) {
        goto LABEL_400;
    } else if (dls1_.miter == 5) {
        goto LABEL_500;
    }
// If MITER = 1, call JAC and multiply by scalar. -----------------------
LABEL_100:
    lenp = dls1_.n * dls1_.n;
    for (i = 1; i <= lenp; ++i) {
        ARRAYF(wm, i+2) = 0.0;
    }
    jac(neq, dls1_.tn, y, 0, 0, &ARRAYF(wm, 3), dls1_.n, user_data);
    con = -hl0;
    for (i = 1; i <= lenp; ++i) {
        ARRAYF(wm, i+2) *= con;
    }
    goto LABEL_240;
// If MITER = 2, make N calls to F to approximate J. --------------------
LABEL_200:
    fac = DMNORM(dls1_.n, savf, ewt);
    r0 = 1000.0 * std::abs(dls1_.h) * dls1_.uround * static_cast<double>(dls1_.n) * fac;
    if (r0 == 0.0) r0 = 1.0;
    srur = ARRAYF(wm, 1);
    j1 = 2;
    for (j = 1; j <= dls1_.n; ++j) {
        yj = ARRAYF(y, j);
        r = std::max(srur * std::abs(yj), r0 / ARRAYF(ewt, j));
        ARRAYF(y, j) += r;
        fac = -hl0 / r;
        f(neq, dls1_.tn, y, ftem, user_data);
        for (i = 1; i <= dls1_.n; ++i) {
            ARRAYF(wm, i+j1) = (ARRAYF(ftem, i) - ARRAYF(savf, i)) * fac;
        }
        ARRAYF(y, j) = yj;
        j1 += dls1_.n;
    }
    dls1_.nfe += dls1_.n;
LABEL_240:
// Compute norm of Jacobian. --------------------------------------------
    dlsa_.pdnorm = DFNORM(dls1_.n, &ARRAYF(wm, 3), ewt) / std::abs(hl0);
// Add identity matrix. -------------------------------------------------
    j = 3;
    np1 = dls1_.n + 1;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(wm, j) += 1.0;
        j += np1;
    }
// Do LU decomposition on P. --------------------------------------------
    DGEFA(&ARRAYF(wm, 3), dls1_.n, dls1_.n, &ARRAYF(iwm, 21), ier);
    if (ier != 0) dls1_.ierpj = 1;
    return;
// Dummy block only, since MITER is never 3 in this routine. ------------
LABEL_300:
    return;
// If MITER = 4, call JAC and multiply by scalar. -----------------------
LABEL_400:
    ml = ARRAYF(iwm, 1);
    mu = ARRAYF(iwm, 2);
    ml3 = ml + 3;
    mband = ml + mu + 1;
    meband = mband + ml;
    lenp = meband * dls1_.n;
    for (i = 1; i <= lenp; ++i) {
        ARRAYF(wm, i+2) = 0.0;
    }
    jac(neq, dls1_.tn, y, ml, mu, &ARRAYF(wm, ml3), meband, user_data);
    con = -hl0;
    for (i = 1; i <= lenp; ++i) {
        ARRAYF(wm, i+2) *= con;
    }
    goto LABEL_570;
// If MITER = 5, make MBAND calls to F to approximate J. ----------------
LABEL_500:
    ml = ARRAYF(iwm, 1);
    mu = ARRAYF(iwm, 2);
    mband = ml + mu + 1;
    mba = std::min(mband, dls1_.n);
    meband = mband + ml;
    meb1 = meband - 1;
    srur = ARRAYF(wm, 1);
    fac = DMNORM(dls1_.n, savf, ewt);
    r0 = 1000.0 * std::abs(dls1_.h) * dls1_.uround * static_cast<double>(dls1_.n) * fac;
    if (r0 == 0.0) r0 = 1.0;
    for (j = 1; j <= mba; ++j) {
        for (i = j; j <= dls1_.n; j += mband) {
            yi = ARRAYF(y, i);
            r = std::max(srur * std::abs(yi), r0 / ARRAYF(ewt, i));
            ARRAYF(y, i) += r;
        }
        f(neq, dls1_.tn, y, ftem, user_data);
        for (jj = j; jj <= dls1_.n; jj += mband) {
            ARRAYF(y, jj) = YH(jj, 1);
            yjj = ARRAYF(y, jj);
            r = std::max(srur * std::abs(yjj), r0 / ARRAYF(ewt, jj));
            fac = -hl0 / r;
            i1 = std::max(jj - mu, 1);
            i2 = std::min(jj + ml, dls1_.n);
            ii = jj * meb1 - ml + 2;
            for (i = i1; i <= i2; ++i) {
                ARRAYF(wm, ii+i) = (ARRAYF(ftem, i) - ARRAYF(savf, i)) * fac;
            }
        }
    }
    dls1_.nfe += mba;
LABEL_570:
// Compute norm of Jacobian. --------------------------------------------
    dlsa_.pdnorm = DBNORM(dls1_.n, &ARRAYF(wm, ml+3), meband, ml, mu, ewt) / std::abs(hl0);
// Add identity matrix. -------------------------------------------------
    ii = mband + 2;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(wm, ii) += 1.0;
        ii += mband;
    }
// Do LU decomposition of P. --------------------------------------------
    DGBFA(&ARRAYF(wm, 3), meband, dls1_.n, ml, mu, &ARRAYF(iwm, 3), ier);
    if (ier != 0) dls1_.ierpj = 1;
    return;
//
#ifdef YH
#undef YH
#endif
}

double Odepack::DMNORM(const int n, const double *v, const double *w)
{
    double vm = 0.0;
    for (int i = 0; i < n; ++i) {
        vm = std::max(vm, std::abs(v[i]) * w[i]);
    }
    return vm;
}


double Odepack::DFNORM(const int n, const double *a, const double *w)
{
#ifndef MATA
#define MATA(i, j) MATF(a, n, i, j)
#endif
//
    int i, j;
    double an, sum;
    an = 0.0;
    for (i = 1; i <= n; ++i) {
        sum = 0.0;
        for (j = 1; j <= n; ++j) {
            sum += std::abs(MATA(i, j)) / ARRAYF(w, j);
        }
        an = std::max(an, sum * ARRAYF(w, i));
    }
    return an;
//
#ifdef MATA
#undef MATA
#endif
}

double Odepack::DBNORM(const int n, const double *a, const int nra, const int ml, const int mu, const double *w)
{
#ifndef MATA
#define MATA(i, j) MATF(a, nra, i, j)
#endif
//
    int i, i1, jlo, jhi, j;
    double an, sum;
    an = 0.0;
    for (i = 1; i <= n; ++i) {
        sum = 0.0;
        i1 = i + mu + 1;
        jlo = std::max(i - ml, 1);
        jhi = std::min(i + mu, n);
        for (j = jlo; j <= jhi; ++j) {
            sum += std::abs(MATA(i1-j, j)) / ARRAYF(w, j);
        }
        an = std::max(an, sum * ARRAYF(w, j));
    }
    return an;
//
#ifdef MATA
#undef MATA
#endif
}

void Odepack::DSRCMA(double *rsav, int *isav, const int job)
{
    return;
}

void Odepack::DRCHEK(const int job, ODEPACK_CONSTRAINT g, const int neq, double *y, double *yh, const int nyh, double *g0, double *g1, double *gx, int *jroot, int &irt, void *user_data)
{
#ifndef YH
#define YH(i, j) MATF(yh, nyh, i, j)
#endif
//
    int i, iflag, jflag;
    double hming, t1, temp1, temp2, x;
    bool zroot;
//
    irt = 0;
    for (i = 1; i <= dlsr_.ngc; ++i) {
        ARRAYF(jroot, i) = 0;
    }
    hming = (std::abs(dls1_.tn) + std::abs(dls1_.h)) * dls1_.uround * 100.0;
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
    dlsr_.t0  = dls1_.tn;
    (*g)(neq, dlsr_.t0, y, dlsr_.ngc, g0, user_data);
    dlsr_.nge = 1;
    zroot     = false;
    for (i = 1; i <= dlsr_.ngc; ++i) {
        if (std::abs(ARRAYF(g0, i)) <= 0.0) zroot = true;
    }
    if (!zroot) goto LABEL_190;
// g has a zero at T.  Look at g at T + (small increment). --------------
    temp2    = std::max(hming / std::abs(dls1_.h), 0.1);
    temp1    = temp2 * dls1_.h;
    dlsr_.t0 = dlsr_.t0 + temp1;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = ARRAYF(y, i) + temp2 * YH(i, 2);
    }
    (*g)(neq, dlsr_.t0, y, dlsr_.ngc, g0, user_data);
    dlsr_.nge = dlsr_.nge + 1;
    zroot     = false;
    for (i = 1; i <= dlsr_.ngc; ++i) {
        if (std::abs(ARRAYF(g0, i)) <= 0.0) zroot = true;
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
    if (dlsr_.irfnd == 0) goto LABEL_260;
// If a root was found on the previous step, evaluate G0 = g(T0). -------
    DINTDY(dlsr_.t0, 0, yh, nyh, y, iflag);
    (*g)(neq, dlsr_.t0, y, dlsr_.ngc, g0, user_data);
    dlsr_.nge = dlsr_.nge + 1;
    zroot     = false;
    for (i = 1; i <= dlsr_.ngc; ++i) {
        if (std::abs(ARRAYF(g0, i)) <= 0.0) zroot = true;
    }
    if (!zroot) goto LABEL_260;
// g has a zero at T0.  Look at g at T + (small increment). -------------
    temp1    = std::abs(hming) * (dls1_.h >= 0 ? 1.0: -1.0);
    dlsr_.t0 = dlsr_.t0 + temp1;
    if ((dlsr_.t0 - dls1_.tn) * dls1_.h > 0.0) goto LABEL_230;
    temp2    = temp1 / dls1_.h;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = ARRAYF(y, i) + temp2 * YH(i, 2);
    }
    goto LABEL_240;
LABEL_230:
    DINTDY(dlsr_.t0, 0, yh, nyh, y, iflag);
LABEL_240:
    (*g)(neq, dlsr_.t0, y, dlsr_.ngc, g0, user_data);
    dlsr_.nge = dlsr_.nge + 1;
    zroot     = false;
    for (i = 1; i <= dlsr_.ngc; ++i) {
        if (std::abs(ARRAYF(g0, i)) > 0.0) continue;
        ARRAYF(jroot, i) = 1;
        zroot            = true;
    }
    if (!zroot) goto LABEL_260;
// g has a zero at T0 and also close to T0.  Return root. ---------------
    irt = 1;
    return;
// G0 has no zero components.  Proceed to check relevant interval. ------
LABEL_260:
    if (dls1_.tn == dlsr_.tlast) goto LABEL_390;
//
LABEL_300:
// Set T1 to TN or TOUTC, whichever comes first, and get g at T1. -------
    if (dlsr_.itaskc == 2 || dlsr_.itaskc == 3 || dlsr_.itaskc || 5) goto LABEL_310;
    if ((dlsr_.toutc - dls1_.tn) * dls1_.h >= 0.0) goto LABEL_310;
    t1 = dlsr_.toutc;
    if ((t1 - dlsr_.t0) * dls1_.h <= 0.0) goto LABEL_390;
    DINTDY(t1, 0, yh, nyh, y, iflag);
    goto LABEL_330;
LABEL_310:
    t1 = dls1_.tn;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = YH(i, 1);
    }
LABEL_330:
    (*g)(neq, t1, y, dlsr_.ngc, g1, user_data);
    dlsr_.nge = dlsr_.nge + 1;
// Call DROOTS to search for root in interval from T0 to T1. ------------
    jflag     = 0;
LABEL_350:
    DROOTS(dlsr_.ngc, hming, jflag, dlsr_.t0, t1, g0, g1, gx, x, jroot);
    if (jflag > 1) goto LABEL_360;
    DINTDY(x, 0, yh, nyh, y, iflag);
    (*g)(neq, x, y, dlsr_.ngc, gx, user_data);
    dlsr_.nge = dlsr_.nge + 1;
    goto LABEL_350;
LABEL_360:
    dlsr_.t0  = x;
    DCOPY(dlsr_.ngc, gx, 1, g0, 1);
    if (jflag == 4) goto LABEL_390;
// Found a root.  Interpolate to X and return. --------------------------
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

void Odepack::DROOTS(const int ng, const double hmin, int &jflag, double &x0, double &x1, double *g0, double *g1, double *gx, double &x, int *jroot)
{
    int i, imxold, nxlast;
    double t2, tmax, fracint, fracsub;
    bool zroot, sgnchg, xroot;
    const double zero  = 0.0;
    const double half  = 0.0;
    const double tenth = 0.1;
    const double five  = 5.0;
//
    if (jflag == 1) goto LABEL_200;
// JFLAG .ne. 1.  Check for change in sign of g or zero at X1. ----------
    dlsr_.imax = 0;
    tmax       = zero;
    zroot      = false;
    for (i = 1; i <= ng; ++i) {
        if (std::abs(ARRAYF(g1, i)) > zero) goto LABEL_110;
        zroot = true;
        continue;
// At this point, G0(i) has been checked and cannot be zero. ------------
LABEL_110:
        if (SIGN(ARRAYF(g0, i)) == SIGN(ARRAYF(g1, i))) continue;
        t2         = std::abs(ARRAYF(g1, i) / (ARRAYF(g1, i) - ARRAYF(g0, i)));
        if (t2 <= tmax) continue;
        tmax       = t2;
        dlsr_.imax = i;
    }
    if (dlsr_.imax > 0) goto LABEL_130;
    sgnchg = false;
    goto LABEL_140;
LABEL_130:
    sgnchg = true;
LABEL_140:
    if (!sgnchg) goto LABEL_400;
// There is a sign change.  Find the first root in the interval. --------
    xroot      = false;
    nxlast     = 0;
    dlsr_.last = 1;
//
// Repeat until the first root in the interval is found.  Loop point. ---
LABEL_150:
    if (xroot) goto LABEL_300;
    if (nxlast == dlsr_.last) goto LABEL_160;
    dlsr_.alpha = 1.0;
    goto LABEL_180;
LABEL_160:
    if (dlsr_.last == 0) goto LABEL_170;
    dlsr_.alpha = 0.5 * dlsr_.alpha;
    goto LABEL_180;
LABEL_170:
    dlsr_.alpha = 0.2 * dlsr_.alpha;
LABEL_180:
    dlsr_.x2    = x1 - (x1 - x0) * ARRAYF(g1, dlsr_.imax) / (ARRAYF(g1, dlsr_.imax) - dlsr_.alpha * ARRAYF(g0, dlsr_.imax));
// If X2 is too close to X0 or X1, adjust it inward, by a fractional ----
// distance that is between 0.1 and 0.5. --------------------------------
    if (std::abs(dlsr_.x2 - x0) < half * dls1_.hmin) {
        fracint  = std::abs(x1 - x0) / dls1_.hmin;
        fracsub  = tenth;
        if (fracint <= five) fracsub = half / fracint;
        dlsr_.x2 = x0 + fracsub * (x1 - x0);
    }
    if (std::abs(x1 - dlsr_.x2) < half * dls1_.hmin) {
        fracint  = std::abs(x1 - x0) / dls1_.hmin;
        fracsub  = tenth;
        if (fracint <= five) fracsub = half / fracint;
        dlsr_.x2 = x1 - fracsub * (x1 - x0);
    }
    jflag = 1;
    x     = dlsr_.x2;
// Return to the calling routine to get a value of GX = g(X). -----------
    return;
// Check to see in which interval g changes sign. -----------------------
LABEL_200:
    imxold     = dlsr_.imax;
    dlsr_.imax = 0;
    tmax       = zero;
    zroot      = false;
    for (i = 1; i <= ng; ++i) {
        if (std::abs(ARRAYF(gx, i)) > zero) goto LABEL_210;
        zroot = true;
        continue;
LABEL_210:
        if (SIGN(ARRAYF(g0, i)) == SIGN(ARRAYF(gx, i))) continue;
        t2         = std::abs(ARRAYF(gx, i) / (ARRAYF(gx, i) - ARRAYF(g0, i)));
        if (t2 <= tmax) continue;
        tmax       = t2;
        dlsr_.imax = i;
    }
    if (dlsr_.imax > 0) goto LABEL_230;
    sgnchg     = false;
    dlsr_.imax = imxold;
    goto LABEL_240;
LABEL_230:
    sgnchg = true;
LABEL_240:
    nxlast = dlsr_.last;
    if (!sgnchg) goto LABEL_250;
// Sign change between X0 and X2, so replace X1 with X2. ----------------
    x1         = dlsr_.x2;
    DCOPY(ng, gx, 1, g1, 1);
    dlsr_.last = 1;
    xroot      = false;
    goto LABEL_270;
LABEL_250:
    if (!zroot) goto LABEL_260;
// Zero value at X2 and no sign change in (X0,X2), so X2 is a root. -----
    x1    = dlsr_.x2;
    DCOPY(ng, gx, 1, g1, 1);
    xroot = true;
    goto LABEL_270;
// No sign change between X0 and X2.  Replace X0 with X2. ---------------
LABEL_260:
    DCOPY(ng, gx, 1, g0, 1);
    x0         = dlsr_.x2;
    dlsr_.last = 0;
    xroot      = false;
LABEL_270:
    if (std::abs(x1 - x0) <= dls1_.hmin) xroot = true;
    goto LABEL_150;
//
// Return with X1 as the root.  Set JROOT.  Set X = X1 and GX = G1. -----
LABEL_300:
    jflag = 2;
    x     = x1;
    DCOPY(ng, g1, 1, gx, 1);
    for (i = 1; i <= ng; ++i) {
        ARRAYF(jroot, i) = 0;
        if (std::abs(ARRAYF(g1, i)) > zero) goto LABEL_310;
        ARRAYF(jroot, i) = 1;
        continue;
LABEL_310:
        if (SIGN(ARRAYF(g0, i)) != SIGN(ARRAYF(g1, i))) ARRAYF(jroot, i) = 1;
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
        ARRAYF(jroot, i) = 0;
        if (std::abs(ARRAYF(g1, i)) <= zero) ARRAYF(jroot, i) = 1;
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

void Odepack::DSRCAR(double *rsav, int *isav, const int job)
{
    return;
}

void Odepack::DSTODPK(const int neq, double *y,double *yh, const int nyh, double *yh1, double *ewt, double *savf,
    double *savx, double *acor, double *wm, int *iwm, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN3 jac, 
    ODEPACK_PSOL psol, void *user_data)
{
#ifndef YH
#define YH(i, j) MATF(yh, dls1_.nyh, i, j)
#endif
#ifndef ELCO
#define ELCO(i, j) MATF(dls1_.elco, 13, i, j)
#endif
#ifndef TESCO
#define TESCO(i, j) MATF(dls1_.tesco, 3, i, j)
#endif
//
    int i, i1, iredo, iret, j, jb, m, ncf, newq;
    double dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup,
           r, rh, rhdn, rhsm, rhup, told;
//
    dls1_.kflag = 0;
    told        = dls1_.tn;
    ncf         = 0;
    dls1_.ierpj = 0;
    dls1_.iersl = 0;
    dls1_.jcur  = 0;
    dls1_.icf   = 0;
    delp        = 0.0;
    if (dls1_.jstart  >  0) goto LABEL_200;
    if (dls1_.jstart == -1) goto LABEL_100;
    if (dls1_.jstart == -2) goto LABEL_160;
//-----------------------------------------------------------------------
// On the first call, the order is set to 1, and other variables are
// initialized.  RMAX is the maximum ratio by which H can be increased
// in a single step.  It is initially 1.E4 to compensate for the small
// initial H, but then is normally equal to 10.  If a failure
// occurs (in corrector convergence or error test), RMAX is set at 2
// for the next increase.
//-----------------------------------------------------------------------
    dls1_.lmax  = dls1_.maxord + 1;
    dls1_.nq    = 1;
    dls1_.l     = 2;
    dls1_.ialth = 2;
    dls1_.rmax  = 10000.0;
    dls1_.rc    = 0.0;
    dls1_.el0   = 1.0;
    dls1_.crate = 0.7;
    dls1_.hold  = dls1_.h;
    dls1_.meo   = dls1_.meth;
    dls1_.nslp  = 0;
    dls1_.ipup  = dls1_.miter;
    iret        = 3;
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
    dls1_.ipup = dls1_.miter;
    dls1_.lmax = dls1_.maxord + 1;
    if (dls1_.ialth == 1) dls1_.ialth = 2;
    if (dls1_.meth == dls1_.meo) goto LABEL_110;
    DCFODE (dls1_.meth, dls1_.elco, dls1_.tesco);
    dls1_.meo   = dls1_.meth;
    if (dls1_.nq > dls1_.maxord) goto LABEL_120;
    dls1_.ialth = dls1_.l;
    iret        = 1;
    goto LABEL_150;
LABEL_110:
    if (dls1_.nq <= dls1_.maxord) goto LABEL_160;
LABEL_120:
    dls1_.nq = dls1_.maxord;
    dls1_.l  = dls1_.lmax;
    for (i = 1; i <= dls1_.l; ++i) {
        ARRAYF(dls1_.el, i) = ELCO(i, dls1_.nq);
    }
    dls1_.nqnyh = dls1_.nq * dls1_.nyh;
    dls1_.rc    = dls1_.rc * ARRAYF(dls1_.el, 1) / dls1_.el0;
    dls1_.el0   = ARRAYF(dls1_.el, 1);
    dls1_.coint = 0.5 / (dls1_.nq + 2);
    dlpk_.epcon = dls1_.coint * TESCO(2, dls1_.nq);
    ddn         = DVNORM(dls1_.n, savf, ewt) / TESCO(1, dls1_.l);
    exdn        = 1.0 / dls1_.l;
    rhdn        = 1.0 / (1.3 * std::pow(ddn, exdn) + 0.0000013);
    rh          = std::min(rhdn, 1.0);
    iredo       = 3;
    if (dls1_.h == dls1_.hold) goto LABEL_170;
    rh          = std::min(rh,std::abs(dls1_.h / dls1_.hold));
    dls1_.h     = dls1_.hold;
    goto LABEL_175;
//-----------------------------------------------------------------------
// DCFODE is called to get all the integration coefficients for the
// current METH.  Then the EL vector and related constants are reset
// whenever the order NQ is changed, or at the start of the problem.
//-----------------------------------------------------------------------
LABEL_140:
    DCFODE(dls1_.meth, dls1_.elco, dls1_.tesco);
LABEL_150:
    for (i = 1; i <= dls1_.l; ++i) {
        ARRAYF(dls1_.el, i) = ELCO(i, dls1_.nq);
    }
    dls1_.nqnyh = dls1_.nq * dls1_.nyh;
    dls1_.rc    = dls1_.rc * ARRAYF(dls1_.el, 1) / dls1_.el0;
    dls1_.el0   = ARRAYF(dls1_.el, 1);
    dls1_.coint = 0.5 / (dls1_.nq + 2);
    dlpk_.epcon = dls1_.coint * TESCO(2, dls1_.nq);
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
    if (dls1_.h == dls1_.hold) goto LABEL_200;
    rh      = dls1_.h / dls1_.hold;
    dls1_.h = dls1_.hold;
    iredo   = 3;
    goto LABEL_175;
LABEL_170:
    rh = std::max(rh, dls1_.hmin / std::abs(dls1_.h));
LABEL_175:
    rh = std::min(rh,dls1_.rmax);
    rh = rh / std::max(1.0, std::abs(dls1_.h) * dls1_.hmxi * rh);
    r  = 1.0;
    for (j = 2; j <= dls1_.l; ++j) {
        r = r * rh;
        for (i = 1; i <= dls1_.n; ++i) {
            YH(i, j) = YH(i, j) * r;
        }
    }
    dls1_.h     = dls1_.h * rh;
    dls1_.rc    = dls1_.rc * rh;
    dls1_.ialth = dls1_.l;
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
    if (dlpk_.jacflg != 0) goto LABEL_202;
    dls1_.ipup  = 0;
    dls1_.crate = 0.7;
    goto LABEL_205;
LABEL_202:
    if (std::abs(dls1_.rc-1.0) > dls1_.ccmax) dls1_.ipup = dls1_.miter;
    if (dls1_.nst >= dls1_.nslp + dls1_.msbp) dls1_.ipup = dls1_.miter;
LABEL_205:
    dls1_.tn = dls1_.tn + dls1_.h;
    i1       = dls1_.nqnyh + 1;
    for (jb = 1; jb <= dls1_.nq; ++jb) {
        i1 = i1 - dls1_.nyh;
        for (i = i1; i <= dls1_.nqnyh; ++i) {
            ARRAYF(yh1, i) = ARRAYF(yh, i) + ARRAYF(yh1, i + dls1_.nyh);
        }
    }
//-----------------------------------------------------------------------
// Up to MAXCOR corrector iterations are taken.  A convergence test is
// made on the RMS-norm of each correction, weighted by the error
// weight vector EWT.  The sum of the corrections is accumulated in the
// vector ACOR(i).  The YH array is not altered in the corrector loop.
//-----------------------------------------------------------------------
LABEL_220:  
    m           = 0;
    dlpk_.mnewt = 0;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = YH(i, 1);
    }
    (*f)(neq, dls1_.tn, y, savf, user_data);
    dls1_.nfe = dls1_.nfe + 1;
    if (dls1_.ipup <= 0) goto LABEL_250;
//-----------------------------------------------------------------------
// If indicated, DPKSET is called to update any matrix data needed,
// before starting the corrector iteration.
// IPUP is set to 0 as an indicator that this has been done.
//-----------------------------------------------------------------------
    DPKSET(neq, y, yh1, ewt, acor, savf, wm, iwm, f, jac, user_data);
    dls1_.ipup  = 0;
    dls1_.rc    = 1.0;
    dls1_.nslp  = dls1_.nst;
    dls1_.crate = 0.7;
    if (dls1_.ierpj != 0) goto LABEL_430;
LABEL_250:
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(acor, i) = 0.0;
    }
LABEL_270:
    if (dls1_.miter != 0) goto LABEL_350;
//-----------------------------------------------------------------------
// In the case of functional iteration, update Y directly from
// the result of the last function evaluation.
//-----------------------------------------------------------------------
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(savf, i) = dls1_.h * ARRAYF(savf, i) - YH(i, 2);
        ARRAYF(   y, i) = ARRAYF(savf, i) - ARRAYF(acor, i);
    }
    del = DVNORM(dls1_.n, y, ewt);
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(   y, i) = YH(i, 1) + ARRAYF(dls1_.el, 1) * ARRAYF(savf, i);
        ARRAYF(acor, i) = ARRAYF(savf, i);
    }
    goto LABEL_400;
//-----------------------------------------------------------------------
// In the case of the chord method, compute the corrector error,
// and solve the linear system with that as right-hand side and
// P as coefficient matrix.
//-----------------------------------------------------------------------
LABEL_350:  
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(savx, i) = dls1_.h * ARRAYF(savf, i) - (YH(i, 2) + ARRAYF(acor, i));
    }
    DSOLPK(neq, y, savf, savx, ewt, wm, iwm, f, psol, user_data);
    if (dls1_.iersl < 0) goto LABEL_430;
    if (dls1_.iersl > 0) goto LABEL_410;
    del = DVNORM (dls1_.n, savx, ewt);
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(acor, i) = ARRAYF(acor, i) + ARRAYF(savx, i);
        ARRAYF(   y, i) = YH(i, 1) + ARRAYF(dls1_.el, 1) * ARRAYF(acor, i);
    }
//-----------------------------------------------------------------------
// Test for convergence.  If M .gt. 0, an estimate of the convergence
// rate constant is stored in CRATE, and this is used in the test.
//-----------------------------------------------------------------------
LABEL_400:  
    if (m != 0) dls1_.crate = std::max(0.2 * dls1_.crate, del / delp);
    dcon = del * std::min(1.0, 1.5 * dls1_.crate) / dlpk_.epcon;
    if (dcon <= 1.0) goto LABEL_450;
    m = m + 1;
    if (m == dls1_.maxcor) goto LABEL_410;
    if (m >= 2 && del > 2.0 * delp) goto LABEL_410;
    dlpk_.mnewt = m;
    delp       = del;
    (*f)(neq, dls1_.tn, y, savf, user_data);
    dls1_.nfe = dls1_.nfe + 1;
    goto LABEL_270;
//-----------------------------------------------------------------------
// The corrector iteration failed to converge.
// If MITER .ne. 0 and the Jacobian is out of date, DPKSET is called for
// the next try.  Otherwise the YH array is retracted to its values
// before prediction, and H is reduced, if possible.  If H cannot be
// reduced or MXNCF failures have occurred, exit with KFLAG = -2.
//-----------------------------------------------------------------------
LABEL_410:  
    if (dls1_.miter == 0 || dls1_.jcur == 1 || dlpk_.jacflg == 0) goto LABEL_430;
    dls1_.icf  = 1;
    dls1_.ipup = dls1_.miter;
    goto LABEL_220;
LABEL_430:  
    dls1_.icf  = 2;
    ncf        = ncf + 1;
    dlpk_.ncfn = dlpk_.ncfn + 1;;
    dls1_.rmax = 2.0;
    dls1_.tn   = told;
    i1         = dls1_.nqnyh + 1;
    for (jb = 1; jb <= dls1_.nq; ++jb) {
        i1 = i1 - dls1_.nyh;
        for (i = i1; i <= dls1_.nqnyh; ++i) {
            ARRAYF(yh, i) = ARRAYF(yh, i) - ARRAYF(yh, i + dls1_.nyh);
        }
    }
    if (dls1_.ierpj < 0 || dls1_.iersl < 0) goto LABEL_680;
    if (std::abs(dls1_.h) <= dls1_.hmin * 1.00001) goto LABEL_670;
    if (ncf == dls1_.mxncf) goto LABEL_670;
    rh         = 0.5;
    dls1_.ipup = dls1_.miter;
    iredo      = 1;
    goto LABEL_170;
//-----------------------------------------------------------------------
// The corrector has converged.  JCUR is set to 0
// to signal that the Jacobian involved may need updating later.
// The local error test is made and control passes to statement 500
// if it fails.
//-----------------------------------------------------------------------
LABEL_450:  
    dls1_.jcur = 0;
    if (m == 0) dsm = del / TESCO(2, dls1_.nq);
    if (m  > 0) dsm = DVNORM(dls1_.n, acor, ewt) / TESCO(2, dls1_.nq);
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
    dls1_.kflag = 0;
    iredo       = 0;
    dls1_.nst   = dls1_.nst + 1;
    dls1_.hu    = dls1_.h;
    dls1_.nqu   = dls1_.nq;
    for (j = 1; j <= dls1_.l; ++j) {
        for (i = 1; i <= dls1_.n; ++i) {
            YH(i, j) = YH(i, j) + ARRAYF(dls1_.el, j) * ARRAYF(acor, i);
        }
    }
    dls1_.ialth = dls1_.ialth - 1;
    if (dls1_.ialth == 0) goto LABEL_520;
    if (dls1_.ialth  > 1) goto LABEL_700;
    if (dls1_.l == dls1_.lmax) goto LABEL_700;
    for (i = 1; i <= dls1_.n; ++i) {
        YH(i, dls1_.lmax) = ARRAYF(acor, i);
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
    dls1_.kflag = dls1_.kflag - 1;
    dls1_.tn    = told;
    i1          = dls1_.nqnyh + 1;
    for (jb = 1; jb <= dls1_.nq; ++jb) {
        i1 = i1 - dls1_.nyh;
        for (i = i1; i <= dls1_.nqnyh; ++i) {
            ARRAYF(yh, i) = ARRAYF(yh, i) - ARRAYF(yh, i + dls1_.nyh);
        }
    }
    dls1_.rmax = 2.0;
    if (std::abs(dls1_.h) <= dls1_.hmin * 1.00001) goto LABEL_660;
    if (dls1_.kflag <= -3) goto LABEL_640;
    iredo = 2;
    rhup  = 0.0;
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
    if (dls1_.l == dls1_.lmax) goto LABEL_540;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(savf, i) = ARRAYF(acor, i) - YH(i, dls1_.lmax);
    }
    dup  = DVNORM(dls1_.n, savf, ewt) / TESCO(3, dls1_.nq);
    exup = 1.0 / (static_cast<double>(dls1_.l + 1));
    rhup = 1.0 / (1.4 * std::pow(dup, exup) + 0.0000014);
LABEL_540:
    exsm = 1.0 / static_cast<double>(dls1_.l);
    rhsm = 1.0 / (1.2 * std::pow(dsm, exsm) + 0.0000012);
    rhdn = 0.0;
    if (dls1_.nq == 1) goto LABEL_560;
    ddn  = DVNORM(dls1_.n, &YH(1, dls1_.l), ewt ) /TESCO(1, dls1_.nq);
    exdn = 1.0 / static_cast<double>(dls1_.nq);
    rhdn = 1.0 / (1.3 * std::pow(ddn, exdn) + 0.0000013);
LABEL_560:
    if (rhsm >= rhup) goto LABEL_570;
    if (rhup  > rhdn) goto LABEL_590;
    goto LABEL_580;
LABEL_570:
    if (rhsm < rhdn) goto LABEL_580;
    newq = dls1_.nq;
    rh   = rhsm;
    goto LABEL_620;
LABEL_580:
    newq = dls1_.nq - 1;
    rh   = rhdn;
    if (dls1_.kflag < 0 && rh > 1.0) rh = 1.0;
    goto LABEL_620;
LABEL_590:
    newq = dls1_.l;
    rh   = rhup;
    if (rh < 1.1) goto LABEL_610;
    r = ARRAYF(dls1_.el, dls1_.l) / static_cast<double>(dls1_.l);
    for (i = 1; i <= dls1_.n; ++i) {
        YH(i, newq + 1) = ARRAYF(acor, i) * r;
    }
    goto LABEL_630;
LABEL_610:  
    dls1_.ialth = 3;
    goto LABEL_700;
LABEL_620:
    if ((dls1_.kflag == 0) && (rh < 1.1)) goto LABEL_610;
    if (dls1_.kflag <= -2) rh = std::min(rh, 0.2);
//-----------------------------------------------------------------------
// If there is a change of order, reset NQ, L, and the coefficients.
// In any case H is reset according to RH and the YH array is rescaled.
// Then exit from 690 if the step was OK, or redo the step otherwise.
//-----------------------------------------------------------------------
    if (newq == dls1_.nq) goto LABEL_170;
LABEL_630:
    dls1_.nq = newq;
    dls1_.l  = dls1_.nq + 1;
    iret     = 2;
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
    if (dls1_.kflag == -10) goto LABEL_660;
    rh       = 0.1;
    rh       = std::max(dls1_.hmin/std::abs(dls1_.h), rh);
    dls1_.h  = dls1_.h * rh;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = YH(i, 1);
    }
    (*f)(neq, dls1_.tn, y, savf, user_data);
    dls1_.nfe = dls1_.nfe + 1;
    for (i = 1; i <= dls1_.n; ++i) {
        YH(i, 2) = dls1_.h * ARRAYF(savf, i);
    }
    dls1_.ipup  = dls1_.miter;
    dls1_.ialth = 5;
    if (dls1_.nq == 1) goto LABEL_200;
    dls1_.nq    = 1;
    dls1_.l     = 2;
    iret        = 3;
    goto LABEL_150;
//-----------------------------------------------------------------------
// All returns are made through this section.  H is saved in HOLD
// to allow the caller to change H on the next step.
//-----------------------------------------------------------------------
LABEL_660:
    dls1_.kflag = -1;
    goto LABEL_720;
LABEL_670:
    dls1_.kflag = -2;
    goto LABEL_720;
LABEL_680:
    dls1_.kflag = -3;
    goto LABEL_720;
LABEL_690:
    dls1_.rmax = 10.0;
LABEL_700:
    r = 1.0 / TESCO(2, dls1_.nqu);
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(acor, i) = ARRAYF(acor, i) * r;
    }
LABEL_720:  
    dls1_.hold  = dls1_.h;
    dls1_.jstart = 1;
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

void Odepack::DPKSET(const int neq, double *y, double *ysv, double *ewt, double *ftem, double *savf, double *wm, int *iwm,
    ODEPACK_FUNCTION f, ODEPACK_JACOBIAN3 jac, void *user_data)
{
    int ier;
    double hl0;
//
    dls1_.ierpj = 0;
    dls1_.jcur  = 1;
    hl0         = dls1_.el0 * dls1_.h;
    (*jac)(f, neq, dls1_.tn, y, ysv, ewt, savf, ftem, hl0, &ARRAYF(wm, dlpk_.locwp), 
        &ARRAYF(iwm, dlpk_.lociwp), ier, user_data);
    dls1_.nje = dls1_.nje + 1;
    if (ier == 0) return;
    dls1_.ierpj = 1;
    return;
}

void Odepack::DSOLPK(const int neq, double *y, double *savf, double *x, double *ewt, double *wm, int *iwm,
        ODEPACK_FUNCTION f, ODEPACK_PSOL psol, void *user_data)
{
    int iflag, lb, ldl, lhes, liom, lgmr, lpcg, lp, lq, lr, lv, lw, lwk, lz, maxlp1, npsl;
    double delta, hl0;
//
    dls1_.iersl = 0;
    hl0         = dls1_.h * dls1_.el0;
    delta       = dlpk_.delt * dlpk_.epcon;
    if (dls1_.miter == 1) {
        goto LABEL_100;
    } else if (dls1_.miter == 2) {
        goto LABEL_200;
    } else if (dls1_.miter == 3) {
        goto LABEL_300;
    } else if (dls1_.miter == 4) {
        goto LABEL_400;
    } else {
        goto LABEL_900;
    }
//-----------------------------------------------------------------------
// Use the SPIOM algorithm to solve the linear system P*x = -f.
//-----------------------------------------------------------------------
LABEL_100:
    lv   = 1;
    lb   = lv + dls1_.n * dlpk_.maxl;
    lhes = lb + dls1_.n;
    lwk  = lhes + dlpk_.maxl * dlpk_.maxl;
    DCOPY(dls1_.n, x, 1, &ARRAYF(wm, lb), 1);
    DSCAL(dls1_.n, dlpk_.rsqrtn, ewt, 1);
    DSPIOM (neq, dls1_.tn, y, savf, &ARRAYF(wm, lb), ewt, dls1_.n, dlpk_.maxl, dlpk_.kmp, delta,
        hl0, dlpk_.jpre, dlpk_.mnewt, f, psol, npsl, x, &ARRAYF(wm, lv), &ARRAYF(wm, lhes), iwm,
        liom, &ARRAYF(wm, dlpk_.locwp), &ARRAYF(iwm, dlpk_.lociwp), &ARRAYF(wm, lwk), iflag, user_data);
    dlpk_.nni = dlpk_.nni + 1;
    dlpk_.nli = dlpk_.nli + liom;
    dlpk_.nps = dlpk_.nps + npsl;
    DSCAL(dls1_.n, dlpk_.sqrtn, ewt, 1);
    if (iflag != 0) dlpk_.ncfl = dlpk_.ncfl + 1;
    if (iflag >= 2) dls1_.iersl = 1;
    if (iflag <  0) dls1_.iersl = -1;
    return;
//-----------------------------------------------------------------------
// Use the SPIGMR algorithm to solve the linear system P*x = -f.
//-----------------------------------------------------------------------
LABEL_200:
    maxlp1 = dlpk_.maxl + 1;
    lv     = 1;
    lb     = lv + dls1_.n * dlpk_.maxl;
    lhes   = lb + dls1_.n + 1;
    lq     = lhes + dlpk_.maxl * maxlp1;
    lwk    = lq + 2 * dlpk_.maxl;
    ldl    = lwk + std::min(1, dlpk_.maxl - dlpk_.kmp) * dls1_.n;
    DCOPY(dls1_.n, x, 1, &ARRAYF(wm, lb), 1);
    DSCAL(dls1_.n, dlpk_.rsqrtn, ewt, 1);
    DSPIGMR(neq, dls1_.tn, y, savf, &ARRAYF(wm, lb), ewt, dls1_.n, dlpk_.maxl, maxlp1, dlpk_.kmp,
        delta, hl0, dlpk_.jpre, dlpk_.mnewt, f, psol, npsl, x, &ARRAYF(wm, lv), &ARRAYF(wm, lhes),
        &ARRAYF(wm, lq), lgmr, &ARRAYF(wm, dlpk_.locwp), &ARRAYF(iwm, dlpk_.lociwp), &ARRAYF(wm, lwk), 
        &ARRAYF(wm, ldl), iflag, user_data);
    dlpk_.nni = dlpk_.nni + 1;
    dlpk_.nli = dlpk_.nli + lgmr;
    dlpk_.nps = dlpk_.nps + npsl;
    DSCAL(dls1_.n, dlpk_.sqrtn, ewt, 1);
    if (iflag != 0) dlpk_.ncfl  = dlpk_.ncfl + 1;
    if (iflag >= 2) dls1_.iersl = 1;
    if (iflag <  0) dls1_.iersl = -1;
    return;
//-----------------------------------------------------------------------
//Use DPCG to solve the linear system P*x = -f
//-----------------------------------------------------------------------
LABEL_300:
    lr  = 1;
    lp  = lr + dls1_.n;
    lw  = lp + dls1_.n;
    lz  = lw + dls1_.n;
    lwk = lz + dls1_.n;
    DCOPY(dls1_.n, x, 1, &ARRAYF(wm, lr), 1);
    DPCG(neq, dls1_.tn, y, savf, &ARRAYF(wm, lr), ewt, dls1_.n, dlpk_.maxl, delta, hl0,
        dlpk_.jpre, dlpk_.mnewt, f, psol, npsl, x, &ARRAYF(wm, lp), &ARRAYF(wm, lw), &ARRAYF(wm, lz),
        lpcg, &ARRAYF(wm, dlpk_.locwp), &ARRAYF(iwm, dlpk_.lociwp), &ARRAYF(wm, lwk), iflag, user_data);
    dlpk_.nni = dlpk_.nni + 1;
    dlpk_.nli = dlpk_.nli + lpcg;
    dlpk_.nps = dlpk_.nps + npsl;
    if (iflag != 0) dlpk_.ncfl  = dlpk_.ncfl + 1;
    if (iflag >= 2) dls1_.iersl = 1;
    if (iflag <  0) dls1_.iersl = -1;
    return;
//-----------------------------------------------------------------------
// Use DPCGS to solve the linear system P*x = -f
//-----------------------------------------------------------------------
LABEL_400:
    lr  = 1;
    lp  = lr + dls1_.n;
    lw  = lp + dls1_.n;
    lz  = lw + dls1_.n;
    lwk = lz + dls1_.n;
    DCOPY(dls1_.n, x, 1, &ARRAYF(wm, lr), 1);
    DPCGS(neq, dls1_.tn, y, savf, &ARRAYF(wm, lr), ewt, dls1_.n, dlpk_.maxl,delta, hl0, 
        dlpk_.jpre, dlpk_.mnewt, f, psol, npsl, x, &ARRAYF(wm, lp), &ARRAYF(wm, lw), &ARRAYF(wm, lz),
        lpcg, &ARRAYF(wm, dlpk_.locwp), &ARRAYF(iwm, dlpk_.lociwp), &ARRAYF(wm, lwk), iflag, user_data);
    dlpk_.nni = dlpk_.nni + 1;
    dlpk_.nli = dlpk_.nli + lpcg;
    dlpk_.nps = dlpk_.nps + npsl;
    if (iflag != 0) dlpk_.ncfl  = dlpk_.ncfl + 1;
    if (iflag >= 2) dls1_.iersl = 1;
    if (iflag <  0) dls1_.iersl = -1;
    return;
//-----------------------------------------------------------------------
// Use DUSOL, which interfaces to psol, to solve the linear system
// (no Krylov iteration).
//-----------------------------------------------------------------------
LABEL_900:
    lb  = 1;
    lwk = lb + dls1_.n;
    DCOPY(dls1_.n, x, 1, &ARRAYF(wm, lb), 1);
    DUSOL(neq, dls1_.tn, y, savf, &ARRAYF(wm, lb), ewt, dls1_.n, delta, hl0, dlpk_.mnewt,
        psol, npsl, x, &ARRAYF(wm, dlpk_.locwp), &ARRAYF(iwm, dlpk_.lociwp), &ARRAYF(wm, lwk), iflag, user_data);
    dlpk_.nni = dlpk_.nni + 1;
    dlpk_.nps = dlpk_.nps + npsl;
    if (iflag != 0) dlpk_.ncfl  = dlpk_.ncfl + 1;
    if (iflag == 3) dls1_.iersl = 1;
    if (iflag <  0) dls1_.iersl = -1;
    return;
}

void Odepack::DSPIOM(const int neq, double tn, double *y, double *savf, double *b, double *wght, const int n, const int maxl, int kmp,
    double &delta, const double hl0, int jpre, int &mnewt, ODEPACK_FUNCTION f, ODEPACK_PSOL psol, int &npsl, double *x,
    double *v, double *hes, int *ipvt, int &liom, double *wp, int *iwp, double *wk, int &iflag, void *user_data)
{
#ifndef V
#define V(i, j) MATF(v, n, i, j)
#endif
#ifndef HES
#define HES(i, j) MATF(hes, maxl, i, j)
#endif
//
    int i, ier, info, j, k, ll, lm1;
    double bnrm, bnrm0, prod, rho, snormw, tem;
//
    iflag = 0;
    liom  = 0;
    npsl  = 0;
//-----------------------------------------------------------------------
// The initial residual is the vector b.  Apply scaling to b, and test
// for an immediate return with X = 0 or X = b.
//-----------------------------------------------------------------------
    for (i = 1; i <= n; ++i) {
        V(i, 1) = ARRAYF(b, i) * ARRAYF(wght, i);
    }
    bnrm0 = DNRM2(n, v, 1);
    bnrm  = bnrm0;
    if (bnrm0 > delta) goto LABEL_30;
    if (mnewt > 0) goto LABEL_20;
    DCOPY(n, b, 1, x, 1);
    return;
LABEL_20:
    for (i = 1; i <= n; ++i) {
        ARRAYF(x, i) = 0.0;
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
        V(i, 1) = ARRAYF(b, i) * ARRAYF(wght, i);
    }
    bnrm  = DNRM2(n, v, 1);
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
        if (ll > 1 && ARRAYF(ipvt, lm1) == lm1) prod = prod * HES(ll, lm1);
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
        ARRAYF(b, k) = 0.0;
    }
    ARRAYF(b, 1) = bnrm;
    DHESL(hes, maxl, ll, ipvt, b);
    for (k = 1; k <= n; ++k) {
        ARRAYF(x, k) = 0.0;
    }
    for (i = 1; i <= ll; ++i) {
        DAXPY(n, ARRAYF(b, i), &V(1, i), 1, x, 1);
    }
    for (i = 1; i <= n; ++i) {
        ARRAYF(x, i) = ARRAYF(x, i) / ARRAYF(wght, i);
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

void Odepack::DATV(const int neq, double *y, double *savf, double *v, double *wght, double *ftem, ODEPACK_FUNCTION f,
    ODEPACK_PSOL psol, double *z, double *vtem, double *wp, int *iwp, double hl0, int &jpre, int &ier, int &npsl,
    void *user_data)
{
    int i;
    double fac, rnorm, tempn;
//
// Set VTEM = D * V.
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(vtem, i) = ARRAYF(v, i) / ARRAYF(wght, i);
    }
    ier = 0;
    if (jpre >= 2) goto LABEL_30;
//
// JPRE = 0 or 1.  Save Y in Z and increment Y by VTEM.
    DCOPY(dls1_.n, y, 1, z, 1);
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = ARRAYF(z, i) + ARRAYF(vtem, i);
    }
    fac = hl0;
    goto LABEL_60;
//
// JPRE = 2 or 3.  Apply inverse of right preconditioner to VTEM.
LABEL_30:
    (*psol)(neq, dls1_.tn, y, savf, ftem, hl0, wp, iwp, vtem, 2, ier, user_data);
    npsl = npsl + 1;
    if (ier != 0) return;
// Calculate L-2 norm of (D-inverse) * VTEM.
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(z, i) = ARRAYF(vtem, i) * ARRAYF(wght,i);
    }
    tempn = DNRM2(dls1_.n, z, 1);
    rnorm = 1.0 / tempn;
// Save Y in Z and increment Y by VTEM/norm.
    DCOPY(dls1_.n, y, 1, z, 1);
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = ARRAYF(z, i) + ARRAYF(vtem, i) * rnorm;
    }
    fac = hl0 * tempn;
//
// For all JPRE, call F with incremented Y argument, and restore Y.
LABEL_60:
    (*f)(neq, dls1_.tn, y, ftem, user_data);
    dls1_.nfe = dls1_.nfe + 1;
    DCOPY(dls1_.n, z, 1, y, 1);
// Set Z = (identity - hl0*Jacobian) * VTEM, using difference quotient.
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(z, i) = ARRAYF(ftem, i) - ARRAYF(savf, i);
    }
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(z, i) = ARRAYF(vtem, i) - fac * ARRAYF(z, i);
    }
// Apply inverse of left preconditioner to Z, if nontrivial.
    if (jpre == 0 || jpre == 2) goto LABEL_85;
    (*psol)(neq, dls1_.n, y, savf, ftem, hl0, wp, iwp, z, 1, ier, user_data);
    npsl = npsl + 1;
    if (ier != 0) return;
LABEL_85:
// Apply D-inverse to Z and return.
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(z, i) = ARRAYF(z, i) * ARRAYF(wght, i);
    }
    return;
}

void Odepack::DORTHOG(double *vnew, double *v, double *hes, const int n, const int ll, const int ldhes, const int kmp, double &snormw)
{
#ifndef V
#define V(i, j) MATF(v, n, i, j)
#endif
#ifndef HES
#define HES(i, j) MATF(hes, ldhes, i, j)
#endif
//
    int i, i0;
    double arg, sumdsq, tem, vnrm;
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

void Odepack::DSPIGMR(const int neq, double tn, double *y, double *savf, double *b, double *wght, int n, const int maxl,
    int maxlp1, int kmp, double &delta, double hl0, int jpre, const int mnewt, ODEPACK_FUNCTION f, ODEPACK_PSOL psol,
    int &npsl, double *x, double *v, double *hes, double *q, int &lgmr, double *wp, int *iwp, double *wk, double *dl,
    int &iflag, void *user_data)
{
#ifndef V
#define V(i, j) MATF(v, n, i, j)
#endif
#ifndef HES
#define HES(i, j) MATF(hes, maxlp1, i, j)
#endif
//
    int i, ier, info, ip1, i2, j, k, ll, llp1;
    double bnrm, bnrm0, c, dlnrm, prod, rho, s, snormw, tem;
//
    iflag = 0;
    lgmr  = 0;
    npsl  = 0;
//-----------------------------------------------------------------------
// The initial residual is the vector b.  Apply scaling to b, and test
// for an immediate return with X = 0 or X = b.
//-----------------------------------------------------------------------
    for (i = 1; i <= n; ++i) {
        V(i, 1) = ARRAYF(b, i) * ARRAYF(wght, i);
    }
    bnrm0 = DNRM2(n, v, 1);
    bnrm  = bnrm0;
    if (bnrm0 > delta) goto LABEL_30;
    if (mnewt > 0) goto LABEL_20;
    DCOPY(n, b, 1, x, 1);
    return;
LABEL_20:
    for (i = 1; i <= n; ++i) {
        ARRAYF(x, i) = 0.0;
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
        V(i, 1) = ARRAYF(b, i) * ARRAYF(wght, i);
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
        prod = prod * ARRAYF(q, 2 * ll);
        rho = std::abs(prod * bnrm);
        if (ll > kmp && kmp < maxl) {
            if (ll == (kmp + 1)) {
                DCOPY(n, &V(1, 1), 1, dl, 1);
                for (i = 1; i <= kmp; ++i) {
                    ip1 = i + 1;
                    i2  = i * 2;
                    s   = ARRAYF(q, i2);
                    c   = ARRAYF(q, i2 - 1);
                    for (k = 1; k <= n; ++k) {
                        ARRAYF(dl, k) = s * ARRAYF(dl, k) + c * V(k, ip1);
                    }
                }
            }
            s = ARRAYF(q, 2 * ll);
            c = ARRAYF(q, 2 * ll - 1) / snormw;
            llp1 = ll + 1;
            for (k = 1; k <= n; ++k) {
                ARRAYF(dl, k) = s * ARRAYF(dl, k) + c * V(k, llp1);
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
        ARRAYF(b, k) = 0.0;
    }
    ARRAYF(b, 1) = bnrm;
    DHELS(hes, maxlp1, ll, q, b);
    for (k = 1; k <= n; ++k) {
        ARRAYF(x, k) = 0.0;
    }
    for (i = 1; i <= ll; ++i) {
        DAXPY(n, ARRAYF(b, i), &V(1, i), 1, x, 1);
    }
    for (i = 1; i <= n; ++i) {
        ARRAYF(x, i) = ARRAYF(x, i) / ARRAYF(wght, i);
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

void Odepack::DPCG(const int neq, double tn, double *y, double *savf, double *r, double *wght, int n, const int maxl,
    const double delta, double hl0, const int &jpre, const int mnewt, ODEPACK_FUNCTION f, ODEPACK_PSOL psol,
    int &npsl, double *x, double *p, double *w, double *z, int &lpcg, double *wp,int *iwp, double *wk, int &iflag,
    void *user_data)
{

}

void Odepack::DPCGS(const int neq, double tn, double *y, double *savf, double *r, double *wght, int n, const int maxl,
    const double delta, double hl0, const int jpre, const int mnewt, ODEPACK_FUNCTION f, ODEPACK_PSOL psol,
    int &npsl, double *x, double *p, double *w, double *z, int &lpcg, double *wp, int *iwp, double *wk, int &iflag,
    void *user_data)
{

}

void Odepack::DATP(const int neq, double *y, const double *savf, double *p, double *wght, const double hl0, double *wk,
    ODEPACK_FUNCTION f, double *w, void *user_data)
{

}

void Odepack::DUSOL(const int neq, double tn, double *y, double *savf, double *b, double *wght, int n, const double delta,
    double hl0, const int mnewt, ODEPACK_PSOL psol, int &npsl, double *x, double *wp, int *iwp, double *wk, int &iflag,
    void *user_data)
{

}

void Odepack::DSRCPK(double *rsav, int *isav, const int job)
{

}

void Odepack::DHEFA(double *a, const int lda, const int n, int *ipvt, int &info, const int job)
{

}

void Odepack::DHESL(double *a, const int lda, const int n, const int *ipvt, double *b)
{

}

void Odepack::DHEQR(double *a, const int lda, const int n, double *q, int &info, const int ijob)
{

}

void Odepack::DHELS(double *A, const int lda, const int n, const double *q, double *b)
{

}

void Odepack::DLHIN(const int neq, int n, const double t0, double *y0, const double *ydot, ODEPACK_FUNCTION f, const double tout,
    const double *uround, double *ewt, const int itol, const double *atol, double *y, double *temp, double &h0, 
    int &niter, int &ier)
{

}

} // end namespace odepack_cpp