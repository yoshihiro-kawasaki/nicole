/**
 * @fn opkda2.cpp
 * @note
 * 2025/01/21 kawasaki
*/

#include "odepack.hpp"

namespace odepack_cpp
{

/**
 * @fn DLSODES
 * @brief DLSODE solves the initial-value problem for stiff or
 * nonstiff systems of first-order ODE's,
 * dy/dt = f(t,y),   or, in component form,
 * dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(N)),  i=1,...,N.
 * @date 2025/02/08 kawasaki
 */
void Odepack::DLSODE(ODEPACK_FUNCTION f, const int neq, double *y, double &t, const double tout,
            const int itol, double *rtol, double *atol, const int itask, int &istate,
            const int iopt, double *rwork, const int lrw, int *iwork, const int liw,
            ODEPACK_JACOBIAN1 jac, const int mf, void *user_data)
{
    std::string msg;
    int i, i1, i2, iflag, kgo, leniw, lenrw, lenwm, ml, mu, imxer, lf0;
    double atoli, ayi, h0, hmax, hmx, rh, rtoli, tcrit, tdist, tnext, tol, tolsf, tp, sum, w0, big, size, ewti;
    bool ihit = false;
    const int mxstp0 = 500;
    const int mxhnl0 = 10;
    const int mord[2] = {12, 5};
//-----------------------------------------------------------------------
// Block a.
// This code block is executed on every call.
// It tests *istate and itask for legality and branches appropriately.
// If *istate > 1 but the flag init shows that initialization has not
// yet been done, an error return occurs.
// If istate = 1 and tout = t, return immediately.
//-----------------------------------------------------------------------
//
//***FIRST EXECUTABLE STATEMENT  DLSODE
    if (istate < 1 || istate > 3) goto LABEL_601;
    if (itask < 1 || itask > 5)   goto LABEL_602;
    if (istate == 1) goto LABEL_10;
    if (dls1_.init == 0)  goto LABEL_603;
    if (istate == 2) goto LABEL_200;
    goto LABEL_20;
LABEL_10:
    dls1_.init = 0;
    if (tout == t) return;
//-----------------------------------------------------------------------
// Block B.
// The next code block is executed for the initial call (ISTATE = 1),
// or for a continuation call with parameter changes (ISTATE = 3).
// It contains checking of all inputs and various initializations.
//
// First check legality of the non-optional inputs NEQ, ITOL, IOPT,
// MF, ML, and MU.
//-----------------------------------------------------------------------
LABEL_20:
    if (neq <= 0) goto LABEL_604;
    if (istate == 1) goto LABEL_25;
    if (neq > dls1_.n) goto LABEL_605;
LABEL_25:
    dls1_.n = neq;
    if (itol < 1 || itol > 4) goto LABEL_606;
    if (iopt < 0 || iopt > 1) goto LABEL_607;
    dls1_.meth = mf / 10;
    dls1_.miter = mf - 10 * dls1_.meth;
    if (dls1_.meth < 1 || dls1_.meth > 2)  goto LABEL_608;
    if (dls1_.miter < 0 || dls1_.miter > 5) goto LABEL_608;
    if (dls1_.miter <= 3) goto LABEL_30;
    ml = ARRAYF(iwork, 1);
    mu = ARRAYF(iwork, 2);
    if (ml < 0 || ml >= dls1_.n) goto LABEL_609;
    if (mu < 0 || mu >= dls1_.n) goto LABEL_610;
LABEL_30:
// Next process and check the optional inputs. --------------------------
    if (iopt == 1) goto LABEL_40;
    dls1_.maxord = ARRAYF(mord, dls1_.meth);
    dls1_.mxstep = mxstp0;
    dls1_.mxhnil = mxhnl0;
    if (istate == 1) h0 = 0.0;
    dls1_.hmxi = 0.0;
    dls1_.hmin = 0.0;
    goto LABEL_60;
LABEL_40:
    dls1_.maxord = ARRAYF(iwork, 5);
    if (dls1_.maxord < 0) goto LABEL_611;
    if (dls1_.maxord == 0) dls1_.maxord = 100;
    dls1_.maxord = std::min(dls1_.maxord, ARRAYF(mord, dls1_.meth));
    dls1_.mxstep = ARRAYF(iwork, 6);
    if (dls1_.mxstep < 0) goto LABEL_612;
    if (dls1_.mxstep == 0) dls1_.mxstep = mxstp0;
    dls1_.mxhnil = ARRAYF(iwork, 7);
    if (dls1_.mxhnil < 0) goto LABEL_613;
    if (dls1_.mxhnil == 0) dls1_.mxhnil = mxhnl0;
    if (istate != 1) goto LABEL_50;
    h0 = ARRAYF(rwork, 5);
    if ((tout - t) * h0 < 0.0) goto LABEL_614;
LABEL_50:
    hmax = ARRAYF(rwork, 6);
    if (hmax < 0.0) goto LABEL_615;
    dls1_.hmxi = 0.0;
    if (hmax > 0.0) dls1_.hmxi = 1.0 / hmax;
    dls1_.hmin = ARRAYF(rwork, 7);
    if (dls1_.hmin < 0.0) goto LABEL_616;
//-----------------------------------------------------------------------
// Set work array pointers and check lengths LRW and LIW.
// Pointers to segments of RWORK and IWORK are named by prefixing L to
// the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
// Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
//-----------------------------------------------------------------------
LABEL_60:
    dls1_.lyh = 21;
    if (istate == 1) dls1_.nyh = dls1_.n;
    dls1_.lwm = dls1_.lyh + (dls1_.maxord + 1) * dls1_.nyh;
    if (dls1_.miter == 0) lenwm = 0;
    if (dls1_.miter == 1 || dls1_.miter == 2) lenwm = dls1_.n * dls1_.n + 2;
    if (dls1_.miter == 3) lenwm = dls1_.n + 2;
    if (dls1_.miter >= 4) lenwm = (2*ml + mu + 1) * dls1_.n + 2;
    dls1_.lewt = dls1_.lwm + lenwm;
    dls1_.lsavf = dls1_.lewt + dls1_.n;
    dls1_.lacor = dls1_.lsavf + dls1_.n;
    lenrw  = dls1_.lacor + dls1_.n - 1;
    ARRAYF(iwork, 17) = lenrw;
    dls1_.liwm = 1;
    leniw = 20 + dls1_.n;
    if (dls1_.miter == 0 || dls1_.miter == 3) leniw = 20;
    ARRAYF(iwork, 18) = leniw;
    if (lenrw > lrw) goto LABEL_617;
    if (leniw > liw) goto LABEL_618;
// Check RTOL and ATOL for legality. ------------------------------------
    rtoli = ARRAYF(rtol, 1);
    atoli = ARRAYF(atol, 1);
    for (i = 1; i <= dls1_.n; ++i) {
        if (itol >= 3) rtoli = ARRAYF(rtol, i);
        if (itol == 2 || itol == 4) atoli = ARRAYF(atol, i);
        if (rtoli < 0.0) goto LABEL_619;
        if (atoli < 0.0) goto LABEL_620;
    }
    if (istate == 1) goto LABEL_100;
// If ISTATE = 3, set flag to signal parameter changes to DSTODE. -------
    dls1_.jstart = -1;
    if (dls1_.nq <= dls1_.maxord) goto LABEL_90;
// MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(rwork, i+dls1_.lsavf-1) = ARRAYF(rwork, i + dls1_.lwm - 1);
    }
// Reload WM(1) = RWORK(LWM), since LWM may have changed. ---------------
LABEL_90:
    if (dls1_.miter > 0) ARRAYF(rwork, dls1_.lwm) = std::sqrt(dls1_.uround);
    if (dls1_.n == dls1_.nyh) goto LABEL_200;
// NEQ was reduced.  Zero part of YH to avoid undefined references. -----
    i1 = dls1_.lyh + dls1_.l * dls1_.nyh;
    i2 = dls1_.lyh + (dls1_.maxord + 1) * dls1_.nyh - 1;
    if (i1 > i2) goto LABEL_200;
    for (i = i1; i <= i2; ++i) {
        ARRAYF(rwork, i) = 0.0;
    }
    goto LABEL_200;
//-----------------------------------------------------------------------
// Block C.
// The next block is for the initial call only (ISTATE = 1).
// It contains all remaining initializations, the initial call to F,
// and the calculation of the initial step size.
// The error weights in EWT are inverted after being loaded.
//-----------------------------------------------------------------------
LABEL_100:
    dls1_.uround = DUMACH();
    dls1_.tn = t;
    if (itask != 4 && itask != 5) goto LABEL_110;
    tcrit = ARRAYF(rwork, 1);
    if ((tcrit - tout) * (tout - t) < 0.0) goto LABEL_625;
    if (h0 != 0.0 && (t + h0 - tcrit) * h0 > 0.0) h0 = tcrit - t;
LABEL_110:
    dls1_.jstart = 0;
    if (dls1_.miter > 0) ARRAYF(rwork, dls1_.lwm) = std::sqrt(dls1_.uround);
    dls1_.nhnil  = 0;
    dls1_.nst    = 0;
    dls1_.nje    = 0;
    dls1_.nslast = 0;
    dls1_.hu     = 0.0;
    dls1_.nqu    = 0;
    dls1_.ccmax  = 0.3;
    dls1_.maxcor = 3;
    dls1_.msbp   = 20;
    dls1_.mxncf  = 10;
// Initial call to F.  (LF0 points to YH(*,2).) -------------------------
    lf0 = dls1_.lyh + dls1_.nyh;
    f(neq, t, y, &ARRAYF(rwork, lf0), user_data);
    dls1_.nfe = 1;
// Load the initial value vector in ydls1_.h. --------------------------------
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(rwork, i + dls1_.lyh - 1) = ARRAYF(y, i);
    }
// Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
    dls1_.nq = 1;
    dls1_.h = 1.0;
    DEWSET(dls1_.n, itol, rtol, atol, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lewt));
    for (i = 1; i <= dls1_.n; ++i) {
        if (ARRAYF(rwork, i + dls1_.lewt - 1) <= 0.0) goto LABEL_621;
        ARRAYF(rwork, i + dls1_.lewt - 1) = 1.0 / ARRAYF(rwork, i + dls1_.lewt - 1);
    }
//-----------------------------------------------------------------------
// The coding below computes the step size, H0, to be attempted on the
// first step, unless the user has supplied a value for this.
// First check that TOUT - T differs significantly from zero.
// A scalar tolerance quantity TOL is computed, as MAX(RTOL(I))
// if this is positive, or MAX(ATOL(I)/ABS(Y(I))) otherwise, adjusted
// so as to be between 100*UROUND and 1.0E-3.
// Then the computed value H0 is given by..
//                                     NEQ
// H0**2 = TOL / ( w0**-2 + (1/NEQ) * SUM ( f(i)/ywt(i) )**2  )
//                                     1
// where   w0     = MAX ( ABS(T), ABS(TOUT) ),
//         f(i)   = i-th component of initial value of f,
//         ywt(i) = EWT(i)/TOL  (a weight for y(i)).
// The sign of H0 is inferred from the initial values of TOUT and T.
//-----------------------------------------------------------------------
    if (h0 != 0.0) goto LABEL_180;
    tdist = std::abs(tout - t);
    w0 = std::max(std::abs(t), std::abs(tout));
    if (tdist < 2.0 * dls1_.uround * w0) goto LABEL_622;
    tol = ARRAYF(rtol, 1);
    if (itol <= 2) goto LABEL_140;
    for (i = 1; i <= dls1_.n; ++i) {
        tol = std::max(tol, ARRAYF(rtol, i));
    }
LABEL_140:
    if (tol > 0.0) goto LABEL_160;
    atoli = ARRAYF(atol, 1);
    for (i = 1; i <= dls1_.n; ++i) {
        if (itol == 2 || itol == 4) atoli = ARRAYF(atol, i);
        ayi = std::abs(ARRAYF(y, i));
        if (ayi != 0.0) tol = std::max(tol, atoli / ayi);
    }
LABEL_160:
    tol = std::max(tol, 100.0 * dls1_.uround);
    tol = std::min(tol, 0.001);
    sum = DVNORM(dls1_.n, &ARRAYF(rwork, lf0), &ARRAYF(rwork, dls1_.lewt));
    sum = 1.0 / (tol * w0 * w0) + tol * sum * sum;
    h0  = 1.0 / std::sqrt(sum);
    h0  = std::min(h0, tdist);
    h0  = h0 * ((tout - t >= 0.0) ? 1.0 : -1.0);
// Adjust H0 if necessary to meet HMAX bound. ---------------------------
LABEL_180:
    rh = std::abs(h0) * dls1_.hmxi;
    if (rh > 1.0) h0 /= rh;
// Load H with H0 and scale YH(*,2) by H0. ------------------------------
    dls1_.h = h0;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(rwork, i + lf0 - 1) *= h0;
    }
    goto LABEL_270;
//-----------------------------------------------------------------------
// Block D.
// The next code block is for continuation calls only (ISTATE = 2 or 3)
// and is to check stop conditions before taking a step.
//-----------------------------------------------------------------------
LABEL_200:
    dls1_.nslast = dls1_.nst;
    if (itask == 1) {
        goto LABEL_210;
    } else if (itask == 2) {
        goto LABEL_250;
    } else if (itask == 3) {
        goto LABEL_220;
    } else if (itask == 4) {
        goto LABEL_230;
    } else if (itask == 5) {
        goto LABEL_240;
    }
LABEL_210:
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_250;
    DINTDY(tout, 0, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, y, iflag);
    if (iflag != 0) goto LABEL_627;
    t = tout;
    goto LABEL_420;
LABEL_220:
    tp = dls1_.tn - dls1_.hu * (1.0 + 100.0 * dls1_.uround);
    if ((tp  - tout) * dls1_.h > 0.0) goto LABEL_623;
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_250;
    goto LABEL_400;
LABEL_230:
    tcrit = ARRAYF(rwork, 1);
    if ((dls1_.tn - tcrit) * dls1_.h > 0.0)  goto LABEL_624;
    if ((tcrit - tout) * dls1_.h < 0.0) goto LABEL_625;
    if ((dls1_.tn - tout) * dls1_.h < 0.0)   goto LABEL_245;
    DINTDY(tout, 0, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, y, iflag);
    if (iflag != 0) goto LABEL_627;
    t = tout;
    goto LABEL_420;
LABEL_240:
    tcrit = ARRAYF(rwork, 1);
    if ((dls1_.tn - tcrit) * dls1_.h > 0.0) goto LABEL_624;    
LABEL_245:
    hmx  = std::abs(dls1_.tn) + std::abs(dls1_.h);
    ihit = std::abs(dls1_.tn - tcrit) <= (100.0 * dls1_.uround * hmx);
    if (ihit) goto LABEL_400;
    tnext = dls1_.tn + dls1_.h * (1.0 + 4.0 * dls1_.uround);
    if ((tnext - tcrit) * dls1_.h <= 0.0) goto LABEL_250;
    dls1_.h = (tcrit - dls1_.tn) * (1.0 - 4.0 * dls1_.uround);
    if (istate == 2) dls1_.jstart = -2;
//-----------------------------------------------------------------------
// Block E.
// The next block is normally executed for all calls and contains
// the call to the one-step core integrator DSTODE.

// This is a looping point for the integration steps.

// First check for too many steps being taken, update EWT (if not at
// start of problem), check for too much accuracy being requested, and
// check for H below the roundoff level in T.
//-----------------------------------------------------------------------
LABEL_250:
    if ((dls1_.nst - dls1_.nslast) >= dls1_.mxstep) goto LABEL_500;
    DEWSET(dls1_.n, itol, rtol, atol, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lewt));
    for (i = 1; i <= dls1_.n; ++i) {
        if (ARRAYF(rwork, i + dls1_.lewt - 1) <= 0.0) goto LABEL_510;
        ARRAYF(rwork, i + dls1_.lewt - 1) = 1.0 / ARRAYF(rwork, i + dls1_.lewt - 1);
    }
LABEL_270:
    tolsf = dls1_.uround * DVNORM(dls1_.n, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lewt));
    if (tolsf <= 1.0) goto LABEL_280;
    tolsf *= 2.0;
    if (dls1_.nst == 0) goto LABEL_626;
    goto LABEL_520;
LABEL_280:
    if ((dls1_.tn + dls1_.h) != dls1_.tn) goto LABEL_290;
    dls1_.nhnil++;
    if (dls1_.nhnil > dls1_.mxhnil) goto LABEL_290;
    msg = "DLSODE-  Warning..internal T (=R1) and H (=R2) are";
    XERRWD(msg, 50, 101, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      such that in the machine, T + H = T on the next step  ";
    XERRWD(msg, 60, 101, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      (H = step size). Solver will continue anyway";
    XERRWD(msg, 50, 101, 0, 0, 0, 0, 2, dls1_.tn, dls1_.h);
    if (dls1_.nhnil < dls1_.mxhnil) goto LABEL_290;
    msg = "DLSODE-  Above warning has been issued I1 times.  ";
    XERRWD(msg, 50, 102, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      It will not be issued again for this problem";
    XERRWD(msg, 50, 102, 0, 1, dls1_.mxhnil, 0, 0, 0.0, 0.0);
LABEL_290:
//-----------------------------------------------------------------------
//  CALL DSTODE(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,DPREPJ,DSOLSY)
//-----------------------------------------------------------------------
    DSTODE(neq, y, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, &ARRAYF(rwork, dls1_.lyh), 
        &ARRAYF(rwork, dls1_.lewt), &ARRAYF(rwork, dls1_.lsavf), &ARRAYF(rwork, dls1_.lacor), 
        &ARRAYF(rwork, dls1_.lwm), &ARRAYF(iwork, dls1_.liwm), f, jac, &Odepack::DPREPJ, &Odepack::DSOLSY, 
        user_data);
    kgo = 1 - dls1_.kflag;
    if (kgo == 1) {
        goto LABEL_300;
    } else if (kgo == 2) {
        goto LABEL_530;
    } else if (kgo == 3) {
        goto LABEL_540;
    }
//-----------------------------------------------------------------------
// Block F.
// The following block handles the case of a successful return from the
// core integrator (KFLAG = 0).  Test for stop conditions.
//-----------------------------------------------------------------------
LABEL_300:
    dls1_.init = 1;
    if (itask == 1) {
        goto LABEL_310;
    } else if (itask == 2) {
        goto LABEL_400;
    } else if (itask == 3) {
        goto LABEL_330;
    } else if (itask == 4) {
        goto LABEL_340;
    } else if (itask == 5) {
        goto LABEL_350;
    }
// ITASK = 1.  If TOUT has been reached, interpolate. -------------------
LABEL_310:
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_250;
    DINTDY(tout, 0, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, y, iflag);
    t = tout;
    goto LABEL_420;
// ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
LABEL_330:
    if ((dls1_.tn - tout) * dls1_.h >= 0.0) goto LABEL_400;
    goto LABEL_250;
// ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
LABEL_340:
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_345;
    DINTDY(tout, 0, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, y, iflag);
    t = tout;
    goto LABEL_420;
LABEL_345:
    hmx = std::abs(dls1_.tn) + std::abs(dls1_.h);
    ihit = std::abs(dls1_.tn - tcrit) <= (100.0 * dls1_.uround * hmx);
    if (ihit) goto LABEL_400;
    tnext = dls1_.tn + dls1_.h * (1.0 + 4.0 * dls1_.uround);
    if ((tnext - tcrit) * dls1_.h <= 0.0) goto LABEL_250;
    dls1_.h = (tcrit - dls1_.tn) * (1.0 - 4.0 * dls1_.uround);
    dls1_.jstart = -2;
    goto LABEL_250;
// ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
LABEL_350:
    hmx  = std::abs(dls1_.tn) + std::abs(dls1_.h);
    ihit = std::abs(dls1_.tn - tcrit) <= (100.0 * dls1_.uround * hmx);
//-----------------------------------------------------------------------
// Block G.
// The following block handles all successful returns from DLSODE.
// If ITASK .NE. 1, Y is loaded from YH and T is set accordingly.
// ISTATE is set to 2, and the optional outputs are loaded into the
// work arrays before returning.
//-----------------------------------------------------------------------
LABEL_400:
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = ARRAYF(rwork, i + dls1_.lyh - 1);
    }
    t = dls1_.tn;
    if (itask != 4 && itask != 5) goto LABEL_420;
    if (ihit) t = tcrit;
LABEL_420:
    istate            = 2;
    ARRAYF(rwork, 11) = dls1_.hu;
    ARRAYF(rwork, 12) = dls1_.h;
    ARRAYF(rwork, 13) = dls1_.tn;
    ARRAYF(iwork, 11) = dls1_.nst;
    ARRAYF(iwork, 12) = dls1_.nfe;
    ARRAYF(iwork, 13) = dls1_.nje;
    ARRAYF(iwork, 14) = dls1_.nqu;
    ARRAYF(iwork, 15) = dls1_.nq;
    return;
//-----------------------------------------------------------------------
// Block H.
// The following block handles all unsuccessful returns other than
// those for illegal input.  First the error message routine is called.
// If there was an error test or convergence test failure, IMXER is set.
// Then Y is loaded from YH and T is set to TN.  The optional outputs
// are loaded into the work arrays before returning.
//-----------------------------------------------------------------------
// The maximum number of steps was taken before reaching TOUT. ----------
LABEL_500:
    msg = "DLSODE-  At current T (=R1), MXSTEP (=I1) steps   ";
    XERRWD(msg, 50, 201, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      taken on this call before reaching TOUT     ";
    XERRWD(msg, 50, 201, 0, 1, dls1_.mxstep, 0, 1, dls1_.tn, 0.0);
    istate = -1;
    goto LABEL_580;
// EWT(I) .LE. 0.0 for some I (not at start of problem). ----------------
LABEL_510:
    ewti = ARRAYF(rwork, i+dls1_.lewt-1);
    msg = "DLSODE-  At T (=R1), EWT(I1) has become R2 .LE. 0.";
    XERRWD(msg, 50, 202, 0, 1, i, 0, 2, dls1_.tn, ewti);
    istate = -6;
    goto LABEL_580;
// Too much accuracy requested for machine precision. -------------------
LABEL_520:
    msg = "DLSODE-  At T (=R1), too much accuracy requested  ";
    XERRWD(msg, 50, 203, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      for precision of machine..  see TOLSF (=R2) ";
    XERRWD(msg, 50, 203, 0, 0, 0, 0, 2, dls1_.tn, tolsf);
    ARRAYF(rwork, 14) = tolsf;
    istate = -2;
    goto LABEL_580;
// KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
LABEL_530:
    msg = "DLSODE-  At T(=R1) and step size H(=R2), the error";
    XERRWD(msg, 50, 204, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      test failed repeatedly or with ABS(H) = HMIN'";
    XERRWD(msg, 50, 204, 0, 0, 0, 0, 2, dls1_.tn, dls1_.h);
    istate = -4;
    goto LABEL_560;
// KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
LABEL_540:
    msg = "DLSODE-  At T (=R1) and step size H (=R2), the    ";
    XERRWD(msg, 50, 205, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      corrector convergence failed repeatedly     ";
    XERRWD(msg, 50, 205, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      or with ABS(H) = HMIN   ";
    XERRWD(msg, 30, 205, 0, 0, 0, 0, 2, dls1_.tn, dls1_.h);
    istate = -5;
// Compute IMXER if relevant. -------------------------------------------
LABEL_560:
    big = 0.0;
    imxer = 1;
    for (i = 1; i <= dls1_.n; ++i) {
        size = std::abs(ARRAYF(rwork, i + dls1_.lacor - 1) * ARRAYF(rwork, i + dls1_.lewt - 1));
        if (big >= size) continue;
        big = size;
        imxer = i;
    }
    ARRAYF(iwork, 16) = imxer;
// Set Y vector, T, and optional outputs. -------------------------------
LABEL_580:
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = ARRAYF(rwork, i + dls1_.lyh - 1);
    }
    t = dls1_.tn;
    ARRAYF(rwork, 11) = dls1_.hu;
    ARRAYF(rwork, 12) = dls1_.h;
    ARRAYF(rwork, 13) = dls1_.tn;
    ARRAYF(iwork, 11) = dls1_.nst;
    ARRAYF(iwork, 12) = dls1_.nfe;
    ARRAYF(iwork, 13) = dls1_.nje;
    ARRAYF(iwork, 14) = dls1_.nqu;
    ARRAYF(iwork, 15) = dls1_.nq;
    return;
//-----------------------------------------------------------------------
// Block I.
// The following block handles all error returns due to illegal input
// (ISTATE = -3), as detected before calling the core integrator.
// First the error message routine is called.  If the illegal input 
// is a negative ISTATE, the run is aborted (apparent infinite loop).
//-----------------------------------------------------------------------
LABEL_601:
    msg = "DLSODE-  ISTATE (=I1) illegal";
    XERRWD(msg, 30, 1, 0, 1, istate, 0, 0, 0.0, 0.0);
    if (istate < 0) goto LABEL_800;
    goto LABEL_700;
LABEL_602:
    msg = "DLSODE-  ITASK (=I1) illegal  ";
    XERRWD(msg, 30, 2, 0, 1, itask, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_603:
    msg = "DLSODE-  ISTATE .GT. 1 but DLSODE not initialized ";
    XERRWD(msg, 50, 3, 0, 0, 0, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_604:
    msg = "DLSODE-  NEQ (=I1) .LT. 1     ";
    XERRWD(msg, 30, 4, 0, 1, neq, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_605:
    msg = "DLSODE-  ISTATE = 3 and NEQ increased (I1 to I2)  ";
    XERRWD(msg, 50, 5, 0, 2, dls1_.n, neq, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_606:
    msg = "DLSODE-  ITOL (=I1) illegal   ";
    XERRWD(msg, 30, 6, 0, 1, itol, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_607:
    msg = "DLSODE-  IOPT (=I1) illegal   ";
    XERRWD(msg, 30, 7, 0, 1, iopt, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_608:
    msg = "DLSODE-  MF (=I1) illegal     ";
    XERRWD(msg, 30, 8, 0, 1, mf, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_609:
    msg = "DLSODE-  ML (=I1) illegal.. .LT.0 or .GE.NEQ (=I2)";
    XERRWD(msg, 50, 9, 0, 2, ml, neq, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_610:
    msg = "DLSODE-  MU (=I1) illegal.. .LT.0 or .GE.NEQ (=I2)";
    XERRWD(msg, 50, 10, 0, 2, mu, neq, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_611:
    msg = "DLSODE-  MAXORD (=I1) .LT. 0  ";
    XERRWD(msg, 30, 11, 0, 1, dls1_.maxord, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_612:
    msg = "DLSODE-  MXSTEP (=I1) .LT. 0  ";
    XERRWD(msg, 30, 12, 0, 1, dls1_.mxstep, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_613:
    msg = "DLSODE-  MXHNIL (=I1) .LT. 0  ";
    XERRWD(msg, 30, 13, 0, 1, dls1_.mxhnil, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_614:
    msg = "DLSODE-  TOUT (=R1) behind T (=R2)      ";
    XERRWD(msg, 40, 14, 0, 0, 0, 0, 2, tout, t);
    msg = "      Integration direction is given by H0 (=R1)  ";
    XERRWD(msg, 50, 14, 0, 0, 0, 0, 1, h0, 0.0);
    goto LABEL_700;
LABEL_615:
    msg = "DLSODE-  HMAX (=R1) .LT. 0.0  ";
    XERRWD(msg, 30, 15, 0, 0, 0, 0, 1, hmax, 0.0);
    goto LABEL_700;
LABEL_616:
    msg = "DLSODE-  HMIN (=R1) .LT. 0.0  ";
    XERRWD(msg, 30, 16, 0, 0, 0, 0, 1, dls1_.hmin, 0.0);
    goto LABEL_700;
LABEL_617:
    msg = "DLSODE-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)";
    XERRWD(msg, 60, 17, 0, 2, lenrw, lrw, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_618:
    msg = "DLSODE-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)";
    XERRWD(msg, 60, 18, 0, 2, leniw, liw, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_619:
    msg = "DLSODE-  RTOL(I1) is R1 .LT. 0.0        ";
    XERRWD(msg, 40, 19, 0, 1, i, 0, 1, rtoli, 0.0);
    goto LABEL_700;
LABEL_620:
    msg = "DLSODE-  ATOL(I1) is R1 .LT. 0.0        ";
    XERRWD(msg, 40, 20, 0, 1, i, 0, 1, atoli, 0.0);
    goto LABEL_700;
LABEL_621:
    ewti = ARRAYF(rwork, dls1_.lewt+i);
    msg = "DLSODE-  EWT(I1) is R1 .LE. 0.0         ";
    XERRWD(msg, 40, 21, 0, 1, i, 0, 1, ewti, 0.0);
    goto LABEL_700;
LABEL_622:
    msg = "DLSODE-  TOUT (=R1) too close to T(=R2) to start integration";
    XERRWD(msg, 60, 22, 0, 0, 0, 0, 2, tout, t);
    goto LABEL_700;
LABEL_623:
    msg = "DLSODE-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  ";
    XERRWD(msg, 60, 23, 0, 1, itask, 0, 2, tout, tp);
    goto LABEL_700;
LABEL_624:
    msg = "DLSODE-  ITASK = 4 OR 5 and TCRIT (=R1) behind TCUR (=R2)   ";
    XERRWD(msg, 60, 24, 0, 0, 0, 0, 2, tcrit, dls1_.tn);
    goto LABEL_700;
LABEL_625:
    msg = "DLSODE-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   ";
    XERRWD(msg, 60, 25, 0, 0, 0, 0, 2, tcrit, tout);
    goto LABEL_700;
LABEL_626:
    msg = "DLSODE-  At start of problem, too much accuracy   ";
    XERRWD(msg, 50, 26, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      requested for precision of machine..  See TOLSF (=R1) ";
    XERRWD(msg, 60, 26, 0, 0, 0, 0, 1, tolsf, 0.0);
    goto LABEL_700;
LABEL_627:
    msg = "DLSODE-  Trouble in DINTDY.  ITASK = I1, TOUT = R1";
    XERRWD(msg, 50, 27, 0, 1, itask, 0, 1, tout, 0.0);
LABEL_700:
    istate = -3;
    return;
LABEL_800:
    msg = "DLSODE-  Run aborted.. apparent infinite loop     ";
    XERRWD(msg, 50, 303, 2, 0, 0, 0, 0, 0.0, 0.0);
    return;
}


void Odepack::DLSODES(ODEPACK_FUNCTION f, const int neq, double *y, double &t, const double tout, 
            const int itol, double *rtol, double *atol, const int itask, int &istate, 
            const int iopt, double *rwork, const int lrw, int *iwork, const int liw, 
            ODEPACK_JACOBIAN2 jac, const int mf, void *user_data)
{
std::string msg;
    int i, j, i1, i2, iflag, imax, imul, imxer, ipflag, ipgo, irem,
        kgo, lenyht, leniw, lenrw, lf0, lia, lja,
        lrtem, lwtem, lyhd, lyhn, mf1, ncolm;
    double atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli,
           tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0;
    bool ihit;
    const int mord[2] = {12, 5};
    const int mxstp0  = 500;
    const int mxhnl0  = 10;
//-----------------------------------------------------------------------
// In the Data statement below, set LENRAT equal to the ratio of
// the wordlength for a real number to that for an integer.  Usually,
// LENRAT = 1 for single precision and 2 for double precision.  If the
// true ratio is not an integer, use the next smaller integer (.ge. 1).
//-----------------------------------------------------------------------
   const int lenrat = 2;
//-----------------------------------------------------------------------
// Block A.
// This code block is executed on every call.
// It tests ISTATE and ITASK for legality and branches appropriately.
// If ISTATE .gt. 1 but the flag INIT shows that initialization has
// not yet been done, an error return occurs.
// If ISTATE = 1 and TOUT = T, return immediately.
//-----------------------------------------------------------------------
    if (istate < 1 || istate > 3) goto LABEL_601;
    if (itask < 1 || itask > 5) goto LABEL_602;
    if (istate == 1) goto LABEL_10;
    if (dls1_.init == 0) goto LABEL_603;
    if (istate == 2) goto LABEL_200;
    goto LABEL_20;
LABEL_10:
    dls1_.init = 0;
    if (tout == t) return;
//-----------------------------------------------------------------------
// Block B.
// The next code block is executed for the initial call (ISTATE = 1),
// or for a continuation call with parameter changes (ISTATE = 3).
// It contains checking of all inputs and various initializations.
// If ISTATE = 1, the final setting of work space pointers, the matrix
// preprocessing, and other initializations are done in Block C.
//
// First check legality of the non-optional inputs NEQ, ITOL, IOPT,
// MF, ML, and MU.
//-----------------------------------------------------------------------
LABEL_20:
    if (neq <= 0) goto LABEL_604;
    if (istate == 1) goto LABEL_25;
    if (neq > dls1_.n) goto LABEL_605;
LABEL_25:
    dls1_.n = neq;
    if (itol < 1 || itol > 4) goto LABEL_606;
    if (iopt < 0 || iopt > 1) goto LABEL_607;
    dlss_.moss = mf / 100;
    mf1 = mf - 100 * dlss_.moss;
    dls1_.meth = mf1 / 10;
    dls1_.miter = mf1 - 10 * dls1_.meth;
    if (dlss_.moss < 0 || dlss_.moss > 2) goto LABEL_608;
    if (dls1_.meth < 1 || dls1_.meth > 2) goto LABEL_608;
    if (dls1_.miter < 0 || dls1_.miter > 3) goto LABEL_608;
    if (dls1_.miter == 0 || dls1_.miter == 3) dlss_.moss = 0;
// Next process and check the optional inputs. --------------------------
    if (iopt == 1) goto LABEL_40;
    dls1_.maxord = ARRAYF(mord, dls1_.meth);
    dls1_.mxstep = mxstp0;
    dls1_.mxhnil = mxhnl0;
    if (istate == 1) h0 = 0.0;
    dls1_.hmxi = 0.0;
    dls1_.hmin = 0.0;
    dlss_.seth = 0.0;
    goto LABEL_60;
LABEL_40:
    dls1_.maxord = ARRAYF(iwork, 5);
    if (dls1_.maxord < 0) goto LABEL_611;
    if (dls1_.maxord == 0) dls1_.maxord = 100;
    dls1_.maxord = std::min(dls1_.maxord, ARRAYF(mord, dls1_.meth));
    dls1_.mxstep = ARRAYF(iwork, 6);
    if (dls1_.mxstep < 0) goto LABEL_612;
    if (dls1_.mxstep == 0) dls1_.mxstep = mxstp0;
    dls1_.mxhnil = ARRAYF(iwork, 7);
    if (dls1_.mxhnil < 0) goto LABEL_613;
    if (dls1_.mxhnil == 0) dls1_.mxhnil = mxhnl0;
    if (istate != 1) goto LABEL_50;
    h0 = ARRAYF(rwork, 5);
    if ((tout - t) * h0 < 0.0) goto LABEL_614;
LABEL_50:
    hmax = ARRAYF(rwork, 6);
    if (hmax < 0.0) goto LABEL_615;
    dls1_.hmxi = 0.0;
    if (hmax > 0.0) dls1_.hmxi = 1.0 / hmax;
    dls1_.hmin = ARRAYF(rwork, 7);
    if (dls1_.hmin < 0.0) goto LABEL_616;
    dlss_.seth = ARRAYF(rwork, 8);
    if (dlss_.seth < 0.0) goto LABEL_609;
// Check RTOL and ATOL for legality. ------------------------------------
LABEL_60:
    rtoli = ARRAYF(rtol, 1);
    atoli = ARRAYF(atol, 1);
    for (i = 1; i <= dls1_.n; ++i) {
        if (itol >= 3) rtoli = ARRAYF(rtol, i);
        if (itol == 2 || itol == 4) atoli = ARRAYF(atol, i);
        if (rtoli < 0.0) goto LABEL_619;
        if (atoli < 0.0) goto LABEL_620;
    }
//-----------------------------------------------------------------------
// Compute required work array lengths, as far as possible, and test
// these against LRW and LIW.  Then set tentative pointers for work
// arrays.  Pointers to RWORK/IWORK segments are named by prefixing L to
// the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
// Segments of RWORK (in order) are denoted  WM, YH, SAVF, EWT, ACOR.
// If MITER = 1 or 2, the required length of the matrix work space WM
// is not yet known, and so a crude minimum value is used for the
// initial tests of LRW and LIW, and YH is temporarily stored as far
// to the right in RWORK as possible, to leave the maximum amount
// of space for WM for matrix preprocessing.  Thus if MITER = 1 or 2
// and MOSS .ne. 2, some of the segments of RWORK are temporarily
// omitted, as they are not needed in the preprocessing.  These
// omitted segments are: ACOR if ISTATE = 1, EWT and ACOR if ISTATE = 3
// and MOSS = 1, and SAVF, EWT, and ACOR if ISTATE = 3 and MOSS = 0.
//-----------------------------------------------------------------------
    dlss_.lrat = lenrat;
    if (istate == 1) dls1_.nyh = dls1_.n;
    dlss_.lwmin = 0;
    if (dls1_.miter == 1) dlss_.lwmin = 4 * dls1_.n + 10 * dls1_.n / dlss_.lrat;
    if (dls1_.miter == 2) dlss_.lwmin = 4 * dls1_.n + 11 * dls1_.n / dlss_.lrat;
    if (dls1_.miter == 3) dlss_.lwmin = dls1_.n + 2;
    dlss_.lenyh = (dls1_.maxord + 1) * dls1_.nyh;
    dlss_.lrest = dlss_.lenyh + 3 * dls1_.n;
    lenrw = 20 + dlss_.lwmin + dlss_.lrest;
    ARRAYF(iwork, 17) = lenrw;
    leniw = 30;
    if (dlss_.moss == 0 && dls1_.miter != 0 && dls1_.miter != 3) leniw = leniw + dls1_.n + 1;
    ARRAYF(iwork, 18) = leniw;
    if (lenrw > lrw) goto LABEL_617;
    if (leniw > liw) goto LABEL_618;
    lia = 31;
    if (dlss_.moss == 0 && dls1_.miter != 0 && dls1_.miter != 3) leniw = leniw + ARRAYF(iwork, lia + dls1_.n) - 1;
    ARRAYF(iwork, 18) = leniw;
    if (leniw > liw) goto LABEL_618;
    lja = lia + dls1_.n + 1;
    lia = std::min(lia, liw);
    lja = std::min(lja, liw);
    dls1_.lwm = 21;
    if (istate == 1) dls1_.nq = 1;
    ncolm = std::min(dls1_.nq + 1, dls1_.maxord + 2);
    dlss_.lenyhm = ncolm * dls1_.nyh;
    lenyht = dlss_.lenyh;
    if (dls1_.miter == 1 || dls1_.miter == 2) lenyht = dlss_.lenyhm;
    imul = 2;
    if (istate == 3) imul = dlss_.moss;
    if (dlss_.moss == 2) imul = 3;
    lrtem = lenyht + imul * dls1_.n;
    lwtem = dlss_.lwmin;
    if (dls1_.miter == 1 || dls1_.miter == 2) lwtem = lrw - 20 - lrtem;
    dlss_.lenwk = lwtem;
    lyhn = dls1_.lwm + lwtem;
    dls1_.lsavf = lyhn + lenyht;
    dls1_.lewt = dls1_.lsavf + dls1_.n;
    dls1_.lacor = dls1_.lewt  + dls1_.n;
    dlss_.istatc = istate;
    if (istate == 1) goto LABEL_100;
//-----------------------------------------------------------------------
// ISTATE = 3.  Move YH to its new location.
// Note that only the part of YH needed for the next step, namely
// MIN(NQ+1,MAXORD+2) columns, is actually moved.
// A temporary error weight array EWT is loaded if MOSS = 2.
// Sparse matrix processing is done in DIPREP/DPREP if MITER = 1 or 2.
// If MAXORD was reduced below NQ, then the pointers are finally set
// so that SAVF is identical to YH(*,MAXORD+2).
//-----------------------------------------------------------------------
    lyhd = dls1_.lyh - lyhn;
    imax = lyhn - 1 + dlss_.lenyhm;
// Move YH.  Move right if LYHD < 0; move left if LYHD > 0. -------------
    if (lyhd < 0) {
        for (i = lyhn; i <= imax; ++i) {
            j = imax + lyhn - i;
            ARRAYF(rwork, j) = ARRAYF(rwork, j + lyhd);
        }
    }
    if (lyhd > 0) {
        for (i = lyhn; i <= imax; ++i) {
            ARRAYF(rwork, i) = ARRAYF(rwork, i + lyhd);
        }
    }
LABEL_80:
    dls1_.lyh = lyhn;
    ARRAYF(iwork, 22) = dls1_.lyh;
    if (dls1_.miter == 0 || dls1_.miter == 3) goto LABEL_92;
    if (dlss_.moss  != 2) goto LABEL_85;
// Temporarily load EWT if MITER = 1 or 2 and MOSS = 2. -----------------
    DEWSET(dls1_.n, itol, rtol, atol, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lewt));
    for (i = 1; i <= dls1_.n; ++i) {
        if (ARRAYF(rwork, i + dls1_.lewt - 1) <= 0.0) goto LABEL_621;
        ARRAYF(rwork, i + dls1_.lewt - 1) = 1.0 / ARRAYF(rwork, i + dls1_.lewt - 1);
    }
LABEL_85:
// DIPREP and DPREP do sparse matrix preprocessing if MITER = 1 or 2. ---
    dls1_.lsavf = std::min(dls1_.lsavf, lrw);
    dls1_.lewt  = std::min(dls1_.lewt,  lrw);
    dls1_.lacor = std::min(dls1_.lacor, lrw);
    DIPREP(neq, y, rwork, &ARRAYF(iwork, lia), &ARRAYF(iwork, lja), ipflag, f, jac, user_data);
    lenrw = dls1_.lwm - 1 + dlss_.lenwk + dlss_.lrest;
    ARRAYF(iwork, 17) = lenrw;
    if (ipflag != -1) ARRAYF(iwork, 23) = dlss_.ipian;
    if (ipflag != -1) ARRAYF(iwork, 24) = dlss_.ipjan;
    ipgo = -ipflag + 1;
    if (ipgo == 1) {
        goto LABEL_90;
    } else if (ipgo == 2) {
        goto LABEL_628;
    } else if (ipgo == 3) {
        goto LABEL_629;
    } else if (ipgo == 4) {
        goto LABEL_630;
    } else if (ipgo == 5) {
        goto LABEL_631;
    } else if (ipgo == 6) {
        goto LABEL_632;
    } else if (ipgo == 7) {
        goto LABEL_633;
    }
LABEL_90:
    ARRAYF(iwork, 22) = dls1_.lyh;
    if (lenrw > lrw) goto LABEL_617;
// Set flag to signal parameter changes to DSTODE. ----------------------
LABEL_92:
    dls1_.jstart = -1;
    if (dls1_.n == dls1_.nyh) goto LABEL_200;
// NEQ was reduced.  Zero part of YH to avoid undefined references. -----
    i1 = dls1_.lyh + dls1_.l * dls1_.nyh;
    i2 = dls1_.lyh + (dls1_.maxord + 1) * dls1_.nyh - 1;
    if (i1 > i2) goto LABEL_200;
    for (i = i1; i <= i2; ++i) {
        ARRAYF(rwork, i) = 0.0;
    }
//-----------------------------------------------------------------------
// Block C.
// The next block is for the initial call only (ISTATE = 1).
// It contains all remaining initializations, the initial call to F,
// the sparse matrix preprocessing (MITER = 1 or 2), and the
// calculation of the initial step size.
// The error weights in EWT are inverted after being loaded.
//-----------------------------------------------------------------------
LABEL_100:
    dls1_.lyh = lyhn;
    ARRAYF(iwork, 22) = dls1_.lyh;
    dls1_.tn  = t;
    dls1_.nst = 0;
    dls1_.h   = 1.0;
    dlss_.nnz = 0;
    dlss_.ngp = 0;
    dlss_.nzl = 0;
    dlss_.nzu = 0;
// Load the initial value vector in YH. ---------------------------------
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(rwork, i + dls1_.lyh - 1) = ARRAYF(y, i);
    }
// Initial call to F.  (LF0 points to YH(*,2).) -------------------------
    lf0 = dls1_.lyh + dls1_.nyh;
    (*f)(neq, t, y, &ARRAYF(rwork, lf0), user_data);
    dls1_.nfe = 1;
// Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
    DEWSET(dls1_.n, itol, rtol, atol, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lewt));
    for (i = 1; i <= dls1_.n; ++i) {
        if (ARRAYF(rwork, i + dls1_.lewt - 1) <= 0.0) goto LABEL_621;
        ARRAYF(rwork, i + dls1_.lewt - 1) = 1.0 / ARRAYF(rwork, i + dls1_.lewt - 1);
    }
    if (dls1_.miter == 0 || dls1_.miter == 3) goto LABEL_120;
// DIPREP and DPREP do sparse matrix preprocessing if MITER = 1 or 2. ---
    dls1_.lacor = std::min(dls1_.lacor, lrw);
    DIPREP(neq, y, rwork,  &ARRAYF(iwork, lia), &ARRAYF(iwork, lja), ipflag, f, jac, user_data);
    lenrw = dls1_.lwm - 1 + dlss_.lenwk + dlss_.lrest;
    ARRAYF(iwork, 17) = lenrw;
    if (ipflag != -1) ARRAYF(iwork, 23) = dlss_.ipian;
    if (ipflag != -1) ARRAYF(iwork, 24) = dlss_.ipjan;
    ipgo = - ipflag + 1;
    if (ipgo == 1) {
        goto LABEL_115;
    } else if (ipgo == 2) {
        goto LABEL_628;
    } else if (ipgo == 3) {
        goto LABEL_629;
    } else if (ipgo == 4) {
        goto LABEL_630;
    } else if (ipgo == 5) {
        goto LABEL_631;
    } else if (ipgo == 6) {
        goto LABEL_632;
    } else if (ipgo == 7) {
        goto LABEL_633;
    }
LABEL_115:
    ARRAYF(iwork, 22) = dls1_.lyh;
    if (lenrw > lrw) goto LABEL_617;
// Check TCRIT for legality (ITASK = 4 or 5). ---------------------------
LABEL_120:
    if (itask != 4 && itask != 5) goto LABEL_125;
    tcrit = ARRAYF(rwork, 1);
    if ((tcrit - tout) * (tout - t) < 0.0) goto LABEL_625;
    if (h0 != 0.0 && (t + h0 - tcrit) * h0 > 0.0) h0 = tcrit - t;
// Initialize all remaining parameters. ---------------------------------
LABEL_125:
    dls1_.uround = DUMACH();
    dls1_.jstart = 0;
    if (dls1_.miter != 0) ARRAYF(rwork, dls1_.lwm) = std::sqrt(dls1_.uround);
    dlss_.msbj   = 50;
    dlss_.nslj   = 0;
    dlss_.ccmxj  = 0.2;
    dlss_.psmall = 1000.0 * dls1_.uround;
    dlss_.rbig   = 0.01 / dlss_.psmall;
    dls1_.nhnil  = 0;
    dls1_.nje    = 0;
    dlss_.nlu    = 0;
    dls1_.nslast = 0;
    dls1_.hu     = 0.0;
    dls1_.nqu    = 0;
    dls1_.ccmax  = 0.3;
    dls1_.maxcor = 3;
    dls1_.msbp   = 20;
    dls1_.mxncf  = 10;
//-----------------------------------------------------------------------
// The coding below computes the step size, H0, to be attempted on the
// first step, unless the user has supplied a value for this.
// First check that TOUT - T differs significantly from zero.
// A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
// if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
// so as to be between 100*UROUND and 1.0E-3.
// Then the computed value H0 is given by..
//                                        NEQ
//     H0**2 = TOL / ( w0**-2 + (1/NEQ) * Sum ( f(i)/ywt(i) )**2  )
//                                         1
// where   w0     = MAX ( ABS(T), ABS(TOUT) ),
//         f(i)   = i-th component of initial value of f,
//         ywt(i) = EWT(i)/TOL  (a weight for y(i)).
// The sign of H0 is inferred from the initial values of TOUT and T.
// ABS(H0) is made .le. ABS(TOUT-T) in any case.
//-----------------------------------------------------------------------
    lf0 = dls1_.lyh + dls1_.nyh;
    if (h0 != 0.0) goto LABEL_180;
    tdist = std::abs(tout - t);
    w0 = std::max(std::abs(t), std::abs(tout));
    if (tdist < 2.0 * dls1_.uround * w0) goto LABEL_622;
    tol = ARRAYF(rtol, 1);
    if (itol <= 2) goto LABEL_140;
    for (i = 1; i <= dls1_.n; ++i) {
        tol = std::max(tol, ARRAYF(rtol, i));
    }
LABEL_140:
    if (tol > 0.0) goto LABEL_160;
    atoli = ARRAYF(atol, 1);
    for (i = 1; i <= dls1_.n; ++i) {
        if (itol == 2 || itol == 4) atoli = ARRAYF(atol, i);
        ayi = std::abs(ARRAYF(y, i));
        if (ayi != 0.0) tol = std::max(tol, atoli / ayi);
    }
LABEL_160:
    tol = std::max(tol, 100.0 * dls1_.uround);
    tol = std::min(tol, 0.001);
    sum = DVNORM(dls1_.n, &ARRAYF(rwork, lf0), &ARRAYF(rwork, dls1_.lewt));
    sum = 1.0 / (tol * w0 * w0) + tol * sum * sum;
    h0  = 1.0 / std::sqrt(sum);
    h0  = std::min(h0, tdist);
    h0  = h0 * ((tout - t) >= 0.0 ? 1.0 : -1.0);
// Adjust H0 if necessary to meet HMAX bound. ---------------------------
LABEL_180:
    rh = std::abs(h0) * dls1_.hmxi;
    if (rh > 1.0) h0 = h0 / rh;
// Load H with H0 and scale YH(*,2) by H0. ------------------------------
    dls1_.h = h0;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(rwork, i + lf0 - 1) = h0 * ARRAYF(rwork, i + lf0 - 1);
    }
    goto LABEL_270;
//-----------------------------------------------------------------------
// Block D.
// The next code block is for continuation calls only (ISTATE = 2 or 3)
// and is to check stop conditions before taking a step.
//-----------------------------------------------------------------------
LABEL_200:
    dls1_.nslast = dls1_.nst;
    if (itask == 1) {
        goto LABEL_210;
    } else if (itask == 2) {
        goto LABEL_250;
    } else if (itask == 3) {
        goto LABEL_220;
    } else if (itask == 4) {
        goto LABEL_230;
    } else if (itask == 5) {
        goto LABEL_240;
    }
LABEL_210:
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_250;
    DINTDY(tout, 0, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, y, iflag);
    if (iflag != 0) goto LABEL_627;
    t = tout;
    goto LABEL_420;
LABEL_220:
    tp = dls1_.tn - dls1_.hu * (1.0 + 100.0 * dls1_.uround);
    if ((tp - tout) * dls1_.h > 0.0)  goto LABEL_623;
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_250;
    goto LABEL_400;
LABEL_230:
    tcrit = ARRAYF(rwork, 1);
    if ((dls1_.tn - tcrit) * dls1_.h > 0.0) goto LABEL_624;
    if ((tcrit - tout) * dls1_.h < 0.0) goto LABEL_625;
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_245;
    DINTDY(tout, 0, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, y, iflag);
    if (iflag != 0) goto LABEL_627;
    t = tout;
    goto LABEL_420;
LABEL_240:
    tcrit = ARRAYF(rwork, 1);
    if ((dls1_.tn - tcrit) * dls1_.h > 0.0) goto LABEL_624;    
LABEL_245:
    hmx  = std::abs(dls1_.tn) + std::abs(dls1_.h);
    ihit = std::abs(dls1_.tn - tcrit) <= (100.0 * dls1_.uround * hmx);
    if (ihit) goto LABEL_400;
    tnext = dls1_.tn + dls1_.h * (1.0 + 4.0 * dls1_.uround);
    if ((tnext - tcrit) * dls1_.h <= 0.0) goto LABEL_250;
    dls1_.h = (tcrit - dls1_.tn) * (1.0 - 4.0 * dls1_.uround);
    if (istate == 2) dls1_.jstart = -2;
//-----------------------------------------------------------------------
// Block E.
// The next block is normally executed for all calls and contains
// C the call to the one-step core integrator DSTODE.

// This is a looping point for the integration steps.

// First check for too many steps being taken, update EWT (if not at
// start of problem), check for too much accuracy being requested, and
// check for H below the roundoff level in T.
//-----------------------------------------------------------------------
LABEL_250:
    if ((dls1_.nst - dls1_.nslast) >= dls1_.mxstep) goto LABEL_500;
    DEWSET(dls1_.n, itol, rtol, atol, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lewt));
    for (i = 1; i <= dls1_.n; ++i) {
        if (ARRAYF(rwork, i + dls1_.lewt - 1) <= 0.0) goto LABEL_510;
        ARRAYF(rwork, i + dls1_.lewt - 1) = 1.0 / ARRAYF(rwork, i + dls1_.lewt - 1);
    }
LABEL_270:
    tolsf = dls1_.uround * DVNORM(dls1_.n, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lewt));
    if (tolsf <= 1.0) goto LABEL_280;
    tolsf = tolsf * 2.0;
    if (dls1_.nst == 0) goto LABEL_626;
    goto LABEL_520;
LABEL_280:
    if ((dls1_.tn + dls1_.h) != dls1_.tn) goto LABEL_290;
    dls1_.nhnil++;
    if (dls1_.nhnil > dls1_.mxhnil) goto LABEL_290;
    msg = "DLSODES- Warning..internal T (=R1) and H (=R2) are";
    XERRWD(msg, 50, 101, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      such that in the machine, T + H = T on the next step  ";
    XERRWD(msg, 60, 101, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      (H = step size). Solver will continue anyway";
    XERRWD(msg, 50, 101, 0, 0, 0, 0, 2, dls1_.tn, dls1_.h);
    if (dls1_.nhnil < dls1_.mxhnil) goto LABEL_290;
    msg = "DLSODES- Above warning has been issued I1 times.  ";
    XERRWD(msg, 50, 102, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      It will not be issued again for this problem";
    XERRWD(msg, 50, 102, 0, 1, dls1_.mxhnil, 0, 0, 0.0, 0.0);
LABEL_290:
//-----------------------------------------------------------------------
// CALL DSTODE(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,WM,F,JAC,DPRJS,DSOLSS)
//-----------------------------------------------------------------------
    DSTODE(neq, y, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lewt), 
        &ARRAYF(rwork, dls1_.lsavf), &ARRAYF(rwork, dls1_.lacor), &ARRAYF(rwork, dls1_.lwm), &ARRAYF(rwork, dls1_.lwm), 
        f, jac, &Odepack::DPRJS, &Odepack::DSOLSS, user_data);
    kgo = 1 - dls1_.kflag;
    if (kgo == 1) {
        goto LABEL_300;
    } else if (kgo == 2) {
        goto LABEL_530;
    } else if (kgo == 3) {
        goto LABEL_540;
    } else if (kgo == 4) {
        goto LABEL_550;
    }
//-----------------------------------------------------------------------
// Block F.
// The following block handles the case of a successful return from the
// core integrator (KFLAG = 0).  Test for stop conditions.
//-----------------------------------------------------------------------
LABEL_300:
    dls1_.init = 1;
    if (itask == 1) {
        goto LABEL_310;
    } else if (itask == 2) {
        goto LABEL_400;
    } else if (itask == 3) {
        goto LABEL_330;
    } else if (itask == 4) {
        goto LABEL_340;
    } else if (itask == 5) {
        goto LABEL_350;
    }
// ITASK = 1.  if TOUT has been reached, interpolate. -------------------
LABEL_310:
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_250;
    DINTDY(tout, 0, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, y, iflag);
    t = tout;
    goto LABEL_420;
// ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
LABEL_330:
    if ((dls1_.tn - tout) * dls1_.h >= 0.0) goto LABEL_400;
    goto LABEL_250;
// ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
LABEL_340:
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_345;
    DINTDY(tout, 0, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, y, iflag);
    t = tout;
    goto LABEL_420;
LABEL_345:
    hmx = std::abs(dls1_.tn) + std::abs(dls1_.h);
    ihit = std::abs(dls1_.tn - tcrit) <= (100.0 * dls1_.uround * hmx);
    if (ihit) goto LABEL_400;
    tnext = dls1_.tn + dls1_.h * (1.0 + 4.0 * dls1_.uround);
    if ((tnext - tcrit) * dls1_.h <= 0.0) goto LABEL_250;
    dls1_.h = (tcrit - dls1_.tn) * (1.0 - 4.0 * dls1_.uround);
    dls1_.jstart = -2;
    goto LABEL_250;
// ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
LABEL_350:
    hmx = std::abs(dls1_.tn) + std::abs(dls1_.h);
    ihit = std::abs(dls1_.tn - tcrit) <= (100.0 * dls1_.uround * hmx);
//-----------------------------------------------------------------------
// Block G.
// The following block handles all successful returns from DLSODES.
// If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
// ISTATE is set to 2, and the optional outputs are loaded into the
// work arrays before returning.
//-----------------------------------------------------------------------
LABEL_400:
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = ARRAYF(rwork, i + dls1_.lyh - 1);
    }
    t = dls1_.tn;
    if (itask != 4 && itask != 5) goto LABEL_420;
    if (ihit) t = tcrit;
LABEL_420:
    istate = 2;
    ARRAYF(rwork, 11) = dls1_.hu;
    ARRAYF(rwork, 12) = dls1_.h;
    ARRAYF(rwork, 13) = dls1_.tn;
    ARRAYF(iwork, 11) = dls1_.nst;
    ARRAYF(iwork, 12) = dls1_.nfe;
    ARRAYF(iwork, 13) = dls1_.nje;
    ARRAYF(iwork, 14) = dls1_.nqu;
    ARRAYF(iwork, 15) = dls1_.nq;
    ARRAYF(iwork, 19) = dlss_.nnz;
    ARRAYF(iwork, 20) = dlss_.ngp;
    ARRAYF(iwork, 21) = dlss_.nlu;
    ARRAYF(iwork, 25) = dlss_.nzl;
    ARRAYF(iwork, 26) = dlss_.nzu;
    return;
//-----------------------------------------------------------------------
// Block H.
// The following block handles all unsuccessful returns other than
// those for illegal input.  First the error message routine is called.
// If there was an error test or convergence test failure, IMXER is set.
// Then Y is loaded from YH and T is set to TN.
// The optional outputs are loaded into the work arrays before returning.
//-----------------------------------------------------------------------
// The maximum number of steps was taken before reaching TOUT.
LABEL_500:
    msg = "DLSODES- At current T (=R1), MXSTEP (=I1) steps   ";
    XERRWD(msg, 50, 201, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      taken on this call before reaching TOUT     ";
    XERRWD(msg, 50, 201, 0, 1, dls1_.mxstep, 0, 1, dls1_.tn, 0.0);
    istate = -1;
    goto LABEL_580;
// EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
LABEL_510:
    ewti = ARRAYF(rwork, dls1_.lewt + i - 1);
    msg = "DLSODES- At T (=R1), EWT(I1) has become R2 .LE. 0.";
    XERRWD(msg, 50, 202, 0, 1, i, 0, 2, dls1_.tn, ewti);
    istate = -6;
    goto LABEL_580;
// Too much accuracy requested for machine precision. -------------------
LABEL_520:
    msg = "DLSODES- At T (=R1), too much accuracy requested  ";
    XERRWD(msg, 50, 203, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      for precision of machine..  see TOLSF (=R2) ";
    XERRWD(msg, 50, 203, 0, 0, 0, 0, 2, dls1_.tn, tolsf);
    ARRAYF(rwork, 14) = tolsf;
    istate = -2;
    goto LABEL_580;
// KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
LABEL_530:
    msg = "DLSODES- At T(=R1) and step size H(=R2), the error";
    XERRWD(msg, 50, 204, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      test failed repeatedly or with ABS(H) = HMIN'";
    XERRWD(msg, 50, 204, 0, 0, 0, 0, 2, dls1_.tn, dls1_.h);
    istate = -4;
    goto LABEL_560;
// KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
LABEL_540:
    msg = "DLSODES- At T (=R1) and step size H (=R2), the    ";
    XERRWD(msg, 50, 205, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      corrector convergence failed repeatedly     ";
    XERRWD(msg, 50, 205, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      or with ABS(H) = HMIN   ";
    XERRWD(msg, 30, 205, 0, 0, 0, 0, 2, dls1_.tn, dls1_.h);
    istate = -5;
    goto LABEL_560;
// KFLAG = -3.  Fatal error flag returned by DPRJS or DSOLSS (CDRV). ----
LABEL_550:
    msg = "DLSODES- At T (=R1) and step size H (=R2), a fatal";
    XERRWD(msg, 50, 207, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      error flag was returned by CDRV (by way of  ";
    XERRWD(msg, 50, 207, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      Subroutine DPRJS or DSOLSS)       ";
    XERRWD(msg, 40, 207, 0, 0, 0, 0, 2, dls1_.tn, dls1_.h);
    istate = -7;
    goto LABEL_580;
// Compute IMXER if relevant.　-------------------------------------------
LABEL_560:
    big = 0.0;
    imxer = 1;
    for (i = 1; i <= dls1_.n; ++i) {
        size = std::abs(ARRAYF(rwork, i + dls1_.lacor - 1) * ARRAYF(rwork, i + dls1_.lewt - 1));
        if (big >= size) continue;
        big = size;
        imxer = i;
    }
    ARRAYF(iwork, 16) = imxer;
// Set Y vector, T, and optional outputs. -------------------------------
    LABEL_580:
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = ARRAYF(rwork, i + dls1_.lyh - 1);
    }
    t = dls1_.tn;
    ARRAYF(rwork, 11) = dls1_.hu;
    ARRAYF(rwork, 12) = dls1_.h;
    ARRAYF(rwork, 13) = dls1_.tn;
    ARRAYF(iwork, 11) = dls1_.nst;
    ARRAYF(iwork, 12) = dls1_.nfe;
    ARRAYF(iwork, 13) = dls1_.nje;
    ARRAYF(iwork, 14) = dls1_.nqu;
    ARRAYF(iwork, 15) = dls1_.nq;
    ARRAYF(iwork, 19) = dlss_.nnz;
    ARRAYF(iwork, 20) = dlss_.ngp;
    ARRAYF(iwork, 21) = dlss_.nlu;
    ARRAYF(iwork, 25) = dlss_.nzl;
    ARRAYF(iwork, 26) = dlss_.nzu;
    return;
//-----------------------------------------------------------------------
// Block I.
// The following block handles all error returns due to illegal input
// (ISTATE = -3), as detected before calling the core integrator.
// First the error message routine is called.  If the illegal input
// is a negative ISTATE, the run is aborted (apparent infinite loop).
//-----------------------------------------------------------------------
LABEL_601:
    msg = "DLSODES- ISTATE (=I1) illegal.";
    XERRWD(msg, 30, 1, 0, 1, istate, 0, 0, 0.0, 0.0);
    if (istate < 0) goto LABEL_800;
    goto LABEL_700;
LABEL_602:
    msg = "DLSODES- ITASK (=I1) illegal. ";
    XERRWD(msg, 30, 2, 0, 1, itask, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_603:
    msg = "DLSODES- ISTATE.gt.1 but DLSODES not initialized. ";
    XERRWD(msg, 50, 3, 0, 0, 0, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_604:
    msg = "DLSODES- NEQ (=I1) .lt. 1     ";
    XERRWD(msg, 30, 4, 0, 1, neq, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_605:
    msg = "DLSODES- ISTATE = 3 and NEQ increased (I1 to I2). ";
    XERRWD(msg, 50, 5, 0, 2, dls1_.n, neq, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_606:
    msg = "DLSODES- ITOL (=I1) illegal.  ";
    XERRWD(msg, 30, 6, 0, 1, itol, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_607:
    msg = "DLSODES- IOPT (=I1) illegal.  ";
    XERRWD(msg, 30, 7, 0, 1, iopt, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_608:
    msg = "DLSODES- MF (=I1) illegal.    ";
    XERRWD(msg, 30, 8, 0, 1, mf, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_609:
    msg = "DLSODES- SETH (=R1) .lt. 0.0  ";
    XERRWD(msg, 30, 9, 0, 0, 0, 0, 1, dlss_.seth, 0.0);
    goto LABEL_700;
LABEL_611:
    msg = "DLSODES- MAXORD (=I1) .lt. 0  ";
    XERRWD(msg, 30, 11, 0, 1, dls1_.maxord, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_612:
    msg = "DLSODES- MXSTEP (=I1) .lt. 0  ";
    XERRWD(msg, 30, 12, 0, 1, dls1_.mxstep, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_613:
    msg = "DLSODES- MXHNIL (=I1) .lt. 0  ";
    XERRWD(msg, 30, 13, 0, 1, dls1_.mxhnil, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_614:
    msg = "DLSODES- TOUT (=R1) behind T (=R2)      ";
    XERRWD(msg, 40, 14, 0, 0, 0, 0, 2, tout, t);
    msg = "      Integration direction is given by H0 (=R1)  ";
    XERRWD(msg, 50, 14, 0, 0, 0, 0, 1, h0, 0.0);
    goto LABEL_700;
LABEL_615:
    msg = "DLSODES- HMAX (=R1) .LT. 0.0  ";
    XERRWD(msg, 30, 15, 0, 0, 0, 0, 1, hmax, 0.0);
    goto LABEL_700;
LABEL_616:
    msg = "DLSODES- HMIN (=R1) .LT. 0.0  ";
    XERRWD(msg, 30, 16, 0, 0, 0, 0, 1, dls1_.hmin, 0.0);
    goto LABEL_700;
LABEL_617:
    msg = "DLSODES- RWORK length is insufficient to proceed. ";
    XERRWD(msg, 50, 17, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)";
    XERRWD(msg, 60, 17, 0, 2, lenrw, lrw, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_618:
    msg = "DLSODES- IWORK length is insufficient to proceed. ";
    XERRWD(msg, 50, 18, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "        Length needed is .ge. LENIW (=I1), exceeds LIW (=I2)";
    XERRWD(msg, 60, 18, 0, 2, leniw, liw, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_619:
    msg = "DLSODES- RTOL(I1) is R1 .LT. 0.0        ";
    XERRWD(msg, 40, 19, 0, 1, i, 0, 1, rtoli, 0.0);
    goto LABEL_700;
LABEL_620:
    msg = "DLSODES- ATOL(I1) is R1 .LT. 0.0        ";
    XERRWD(msg, 40, 20, 0, 1, i, 0, 1, atoli, 0.0);
    goto LABEL_700;
LABEL_621:
    ewti = ARRAYF(rwork, dls1_.lewt + i - 1);
    msg = "DLSODES-  EWT(I1) is R1 .le. 0.0         ";
    XERRWD(msg, 40, 21, 0, 1, i, 0, 1, ewti, 0.0);
    goto LABEL_700;
LABEL_622:
    msg = "DLSODES- TOUT (=R1) too close to T(=R2) to start integration";
    XERRWD(msg, 60, 22, 0, 0, 0, 0, 2, tout, t);
    goto LABEL_700;
LABEL_623:
    msg = "DLSODES- ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  ";
    XERRWD(msg, 60, 23, 0, 1, itask, 0, 2, tout, tp);
    goto LABEL_700;
LABEL_624:
    msg = "DLSODES- ITASK = 4 OR 5 and TCRIT (=R1) behind TCUR (=R2)   ";
    XERRWD(msg, 60, 24, 0, 0, 0, 0, 2, tcrit, dls1_.tn);
    goto LABEL_700;
LABEL_625:
    msg = "DLSODES- ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   ";
    XERRWD(msg, 60, 25, 0, 0, 0, 0, 2, tcrit, tout);
    goto LABEL_700;
LABEL_626:
    msg = "DLSODES- At start of problem, too much accuracy   ";
    XERRWD(msg, 50, 26, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      requested for precision of machine..  See TOLSF (=R1) ";
    XERRWD(msg, 60, 26, 0, 0, 0, 0, 1, tolsf, 0.0);
    goto LABEL_700;
LABEL_627:
    msg = "DLSODES- Trouble in DINTDY.  ITASK = I1, TOUT = R1";
    XERRWD(msg, 50, 27, 0, 1, itask, 0, 1, tout, 0.0);
    goto LABEL_700;
LABEL_628:
    msg = "DLSODES- RWORK length insufficient (for Subroutine DPREP).  ";
    XERRWD(msg, 60, 28, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)";
    XERRWD(msg, 60, 28, 0, 2, lenrw, lrw, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_629:
    msg = "DLSODES- RWORK length insufficient (for Subroutine JGROUP). ";
    XERRWD(msg, 60, 29, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)";
    XERRWD(msg, 60, 29, 0, 2, lenrw, lrw, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_630:
    msg = "DLSODES- RWORK length insufficient (for Subroutine ODRV).   ";
    XERRWD(msg, 60, 30, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)";
    XERRWD(msg, 60, 30, 0, 2, lenrw, lrw, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_631:
    msg = "DLSODES- Error from ODRV in Yale Sparse Matrix Package.     ";
    XERRWD(msg, 60, 31, 0, 0, 0, 0, 0, 0.0, 0.0);
    imul = (dlss_.iys - 1) / dls1_.n;
    irem = dlss_.iys - imul * dls1_.n;
    msg = "      At T (=R1), ODRV returned error flag = I1*NEQ + I2.   ";
    XERRWD(msg, 60, 31, 0, 2, imul, irem, 1, dls1_.tn, 0.0);
    goto LABEL_700;
LABEL_632:
    msg = "DLSODES- RWORK length insufficient (for Subroutine CDRV).   ";
    XERRWD(msg, 60, 32, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)";
    XERRWD(msg, 60, 32, 0, 2, lenrw, lrw, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_633:
    msg = "DLSODES- Error from CDRV in Yale Sparse Matrix Package.     ";
    XERRWD(msg, 60, 33, 0, 0, 0, 0, 0, 0.0, 0.0);
    imul = (dlss_.iys - 1) / dls1_.n;
    irem = dlss_.iys - imul * dls1_.n;
    msg = "      At T (=R1), CDRV returned error flag = I1*NEQ + I2.   ";
    XERRWD(msg, 60, 33, 0, 2, imul, irem, 1, dls1_.tn, 0.0);
    if (imul == 2) {
        msg = "        Duplicate entry in sparsity structure descriptors.  ";
        XERRWD(msg, 60, 33, 0, 0, 0, 0, 0, 0.0, 0.0);
    }
    if (imul == 3 || imul == 6) {
        msg = "        Insufficient storage for NSFC (called by CDRV).     ";
        XERRWD(msg, 60, 33, 0, 0, 0, 0, 0, 0.0, 0.0);
    }
//
LABEL_700:
    istate = -3;
    return;
//
LABEL_800:
    msg = "DLSODES-  Run aborted.. apparent infinite loop     ";
    XERRWD(msg, 50, 303, 2, 0, 0, 0, 0, 0.0, 0.0);
    return;
}

void Odepack::DLSODA(ODEPACK_FUNCTION f, const int neq, double *y, double &t, const double tout,
            const int itol, double *rtol, double *atol, const int itask, int &istate,
            const int iopt, double *rwork, const int lrw, int *iwork, const int liw,
            ODEPACK_JACOBIAN1 jac, const int jt, void *user_data)
{
    int i, i1, i2, iflag, imxer, kgo, lf0, leniw, lenrw, lenwm, ml, mu;
    int len1, len1c, len1n, len1s, len2, leniwc, lenrwc;
    double atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli,
           tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0;
    bool ihit = false;
    const int mxstp0 = 500;
    const int mxhnl0 = 10;
    const int mord[2] = {12, 5};
    std::string msg;
//-----------------------------------------------------------------------
// Block A.
// This code block is executed on every call.
// It tests ISTATE and ITASK for legality and branches appropriately.
// If ISTATE .gt. 1 but the flag INIT shows that initialization has
// not yet been done, an error return occurs.
// If ISTATE = 1 and TOUT = T, return immediately.
//-----------------------------------------------------------------------
    if (istate < 1 || istate > 3) goto LABEL_601;
    if (itask  < 1 || itask  > 5) goto LABEL_602;
    if (istate == 1) goto LABEL_10;
    if (dls1_.init  == 0) goto LABEL_603;
    if (istate == 2) goto LABEL_200;
    goto LABEL_20;
LABEL_10:
    dls1_.init = 0;
    if (tout == t) return;
//-----------------------------------------------------------------------
// Block B.
// The next code block is executed for the initial call (ISTATE = 1),
// or for a continuation call with parameter changes (ISTATE = 3).
// It contains checking of all inputs and various initializations.

// First check legality of the non-optional inputs NEQ, ITOL, IOPT,
// JT, ML, and MU.
//-----------------------------------------------------------------------
LABEL_20:
    if (neq <= 0)    goto LABEL_604;
    if (istate == 1) goto LABEL_25;
    if (neq > dls1_.n)    goto LABEL_605;
LABEL_25:
    dls1_.n = neq;
    if (itol < 1 || itol > 4) goto LABEL_606;
    if (iopt < 0 || iopt > 1) goto LABEL_607;
    if (jt == 3 || jt < 1 || jt > 5) goto LABEL_608;
    dlsa_.jtyp = jt;
    if (jt <= 2) goto LABEL_30;
    ml = ARRAYF(iwork, 1);
    mu = ARRAYF(iwork, 2);
    if (ml < 0 || ml >= dls1_.n) goto LABEL_609;
    if (mu < 0 || mu >= dls1_.n) goto LABEL_610;
LABEL_30:
// Next process and check the optional inputs. --------------------------
    if (iopt == 1) goto LABEL_40;
    dlsa_.ixpr   = 0;
    dls1_.mxstep = mxstp0;
    dls1_.mxhnil = mxhnl0;
    dls1_.hmxi   = 0.0;
    dls1_.hmin   = 0.0;
    if (istate != 1) goto LABEL_60;
    h0      = 0.0;
    dlsa_.mxordn = ARRAYF(mord, 1);
    dlsa_.mxords = ARRAYF(mord, 2);
    goto LABEL_60;
LABEL_40:
    dlsa_.ixpr   = ARRAYF(iwork, 5);
    if (dlsa_.ixpr < 0 || dlsa_.ixpr > 1) goto LABEL_611;
    dls1_.mxstep = ARRAYF(iwork, 6);
    if (dls1_.mxstep < 0) goto LABEL_612;
    if (dls1_.mxstep == 0) dls1_.mxstep = mxstp0;
    dls1_.mxhnil = ARRAYF(iwork, 7);
    if (dls1_.mxhnil < 0) goto LABEL_613;
    if (dls1_.mxhnil == 0) dls1_.mxhnil = mxhnl0;
    if (istate != 1) goto LABEL_50;
    h0      = ARRAYF(rwork, 5);
    dlsa_.mxordn = ARRAYF(iwork, 8);
    if (dlsa_.mxordn < 0) goto LABEL_628;
    if (dlsa_.mxordn == 0) dlsa_.mxordn = 100;
    dlsa_.mxordn = std::min(dlsa_.mxordn, ARRAYF(mord, 1));
    dlsa_.mxords = ARRAYF(iwork, 9);
    if (dlsa_.mxords < 0) goto LABEL_629;
    if (dlsa_.mxords == 0) dlsa_.mxords = 100;
    dlsa_.mxords = std::min(dlsa_.mxords, ARRAYF(mord, 2));
    if ((tout - t) * h0 < 0.0) goto LABEL_614;
LABEL_50:
    hmax    = ARRAYF(rwork, 6);
    if (hmax < 0.0) goto LABEL_615;
    dls1_.hmxi   = 0.0;
    if (hmax > 0.0) dls1_.hmxi = 1.0 / hmax;
    dls1_.hmin   = ARRAYF(rwork, 7);
    if (dls1_.hmin < 0.0) goto LABEL_616;
//-----------------------------------------------------------------------
// Set work array pointers and check lengths LRW and LIW.
// If ISTATE = 1, METH is initialized to 1 here to facilitate the
// checking of work space lengths.
// Pointers to segments of RWORK and IWORK are named by prefixing L to
// the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
// Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
// If the lengths provided are insufficient for the current method,
// an error return occurs.  This is treated as illegal input on the
// first call, but as a problem interruption with ISTATE = -7 on a
// continuation call.  If the lengths are sufficient for the current
// method but not for both methods, a warning message is sent.
//-----------------------------------------------------------------------
LABEL_60:
    if (istate == 1) dls1_.meth = 1;
    if (istate == 1) dls1_.nyh = dls1_.n;
    dls1_.lyh  = 21;
    len1n = 20 + (dlsa_.mxordn + 1) * dls1_.nyh;
    len1s = 20 + (dlsa_.mxords + 1) * dls1_.nyh;
    dls1_.lwm  = len1s + 1;
    if (jt <= 2) lenwm = dls1_.n * dls1_.n + 2;
    if (jt >= 4) lenwm = (2*ml + mu + 1) * dls1_.n + 2;
    len1s = len1s + lenwm;
    len1c = len1n;
    if (dls1_.meth == 2) len1c = len1s;
    len1   = std::max(len1n, len1s);
    len2   = 3 * dls1_.n;
    lenrw  = len1  + len2;
    lenrwc = len1c + len2;
    ARRAYF(iwork, 17) = lenrw;
    dls1_.liwm  = 1;
    leniw  = 20 + dls1_.n;
    leniwc = 20;
    if (dls1_.meth == 2) leniwc = leniw;
    ARRAYF(iwork, 18) = leniw;
    if (istate == 1 && lrw < lenrwc) goto LABEL_617;
    if (istate == 1 && liw < leniwc) goto LABEL_618;
    if (istate == 3 && lrw < lenrwc) goto LABEL_550;
    if (istate == 3 && liw < leniwc) goto LABEL_555;
    dls1_.lewt = len1 + 1;
    dlsa_.insufr = 0;
    if (lrw >= lenrw) goto LABEL_65;
    dlsa_.insufr = 2;
    dls1_.lewt = len1c + 1;
    msg = "DLSODA-  Warning.. RWORK length is sufficient for now, but  ";
    XERRWD(msg, 60, 103, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      may not be later.  Integration will proceed anyway.   ";
    XERRWD(msg, 60, 103, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      Length needed is LENRW = I1, while LRW = I2.";
    XERRWD(msg, 50, 103, 0, 2, lenrw, lrw, 0, 0.0, 0.0);
LABEL_65:
    dls1_.lsavf = dls1_.lewt + dls1_.n;
    dls1_.lacor = dls1_.lsavf + dls1_.n;
    dlsa_.insufi = 0;
    if (liw >= leniw) goto LABEL_70;
    dlsa_.insufi = 2;
    msg = "DLSODA-  Warning.. IWORK length is sufficient for now, but  '";
    XERRWD(msg, 60, 104, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      may not be later.  Integration will proceed anyway.   ";
    XERRWD(msg, 60, 104, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      Length needed is LENIW = I1, while LIW = I2.";
    XERRWD(msg, 50, 104, 0, 2, leniw, liw, 0, 0.0, 0.0);
LABEL_70:
// Check RTOL and ATOL for legality. ------------------------------------
    rtoli = ARRAYF(rtol, 1);
    atoli = ARRAYF(atol, 1);
    for (i = 1; i <= dls1_.n; ++i) {
        if (itol >= 3) rtoli = ARRAYF(rtol, i);
        if (itol == 2 || itol == 4) atoli = ARRAYF(atol, i);
        if (rtoli < 0.0) goto LABEL_619;
        if (atoli < 0.0) goto LABEL_620;
    }
    if (istate == 1) goto LABEL_100;
// If ISTATE = 3, set flag to signal parameter changes to DSTODA. -------
    dls1_.jstart = -1;
    if (dls1_.n == dls1_.nyh) goto LABEL_200;
// NEQ was reduced.  Zero part of YH to avoid undefined references. -----
    i1 = dls1_.lyh + dls1_.l * dls1_.nyh;
    i2 = dls1_.lyh + (dls1_.maxord + 1) * dls1_.nyh - 1;
    if (i1 > i2) goto LABEL_200;
    for (i = i1; i <= i2; ++i) {
        ARRAYF(rwork, i) = 0.0;
    }
    goto LABEL_200;
//-----------------------------------------------------------------------
// Block C.
// The next block is for the initial call only (ISTATE = 1).
// It contains all remaining initializations, the initial call to F,
// and the calculation of the initial step size.
// The error weights in EWT are inverted after being loaded.
//-----------------------------------------------------------------------
LABEL_100:
    dls1_.uround = DUMACH();
    dls1_.tn  = t;
    dlsa_.tsw = t;
    dls1_.maxord = dlsa_.mxordn;
    if (itask != 4 && itask != 5) goto LABEL_110;
    tcrit = ARRAYF(rwork, 1);
    if ((tcrit - tout) * (tout - t) < 0.0) goto LABEL_625;
    if (h0 != 0.0 && (t + h0 - tcrit) * h0 > 0.0) h0 = tcrit - t;
LABEL_110:
    dls1_.jstart = 0;
    dls1_.nhnil  = 0;
    dls1_.nst    = 0;
    dls1_.nje    = 0;
    dls1_.nslast = 0;
    dls1_.hu     = 0.0;
    dls1_.nqu    = 0;
    dlsa_.mused  = 0;
    dls1_.miter  = 0;
    dls1_.ccmax  = 0.30;
    dls1_.maxcor = 3;
    dls1_.msbp   = 20;
    dls1_.mxncf  = 10;
// Initial call to F.  (LF0 points to YH(*,2).) -------------------------
    lf0 = dls1_.lyh + dls1_.nyh;
    f(neq, t, y, &ARRAYF(rwork, lf0), user_data);
    dls1_.nfe = 1;
// Load the initial value vector in YH. ---------------------------------
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(rwork, i+dls1_.lyh-1) = ARRAYF(y, i);
    }
// Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
    dls1_.nq = 1;
    dls1_.h = 1.0;
    DEWSET(dls1_.n, itol, rtol, atol, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lewt));
    for (i = 1; i <= dls1_.n; ++i) {
        if (ARRAYF(rwork, i+dls1_.lewt-1) <= 0.0) goto LABEL_621;
        ARRAYF(rwork, i+dls1_.lewt-1) = 1.0 / ARRAYF(rwork, i+dls1_.lewt-1);
    }
//-----------------------------------------------------------------------
// The coding below computes the step size, H0, to be attempted on the
// first step, unless the user has supplied a value for this.
// First check that TOUT - T differs significantly from zero.
// A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
// if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
// so as to be between 100*UROUND and 1.0E-3.
// Then the computed value H0 is given by:

// H0**(-2)  =  1./(TOL * w0**2)  +  TOL * (norm(F))**2

// where   w0     = MAX ( ABS(T), ABS(TOUT) ),
//         F      = the initial value of the vector f(t,y), and
//         norm() = the weighted vector norm used throughout, given by
//                 the DMNORM function routine, and weighted by the
//                 tolerances initially loaded into the EWT array.
// The sign of H0 is inferred from the initial values of TOUT and T.
// ABS(H0) is made .le. ABS(TOUT-T) in any case.
//-----------------------------------------------------------------------
    if (h0 != 0.0) goto LABEL_180;
    tdist = std::abs(tout - t);
    w0 = std::max(std::abs(t), std::abs(tout));
    if (tdist < 2.0 * dls1_.uround * w0) goto LABEL_622;
    tol = ARRAYF(rtol, 1);
    if (itol <= 2) goto LABEL_140;
    for (i = 1; i <= dls1_.n; ++i) {
        tol = std::max(tol, ARRAYF(rtol, i));
    }
LABEL_140:
    if (tol > 0.0) goto LABEL_160;
    atoli = ARRAYF(atol, 1);
    for (i = 1; i <= dls1_.n; ++i) {
        if (itol == 2 || itol == 4) atoli = ARRAYF(atol, i);
        ayi = std::abs(ARRAYF(y, i));
        if (ayi != 0.0) tol = std::max(tol, atoli / ayi);
    }
LABEL_160:
    tol = std::max(tol, 100.0 * dls1_.uround);
    tol = std::min(tol, 0.001);
    sum = DVNORM(dls1_.n, &ARRAYF(rwork, lf0), &ARRAYF(rwork, dls1_.lewt));
    sum = 1.0 / (tol * w0 * w0) + tol * sum * sum;
    h0  = 1.0 / std::sqrt(sum);
    h0  = std::min(h0, tdist);
    h0  = h0 * ((tout - t >= 0.0) ? 1.0 : -1.0);
// Adjust H0 if necessary to meet HMAX bound. ---------------------------
LABEL_180:
    rh = std::abs(h0) * dls1_.hmxi;
    if (rh > 1.0) h0 /= rh;
// Load H with H0 and scale YH(*,2) by H0. ------------------------------
    dls1_.h = h0;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(rwork, i + lf0 - 1) *= h0;
    }
    goto LABEL_270;
//-----------------------------------------------------------------------
// Block D.
// The next code block is for continuation calls only (ISTATE = 2 or 3)
// and is to check stop conditions before taking a step.
//-----------------------------------------------------------------------
LABEL_200:
    dls1_.nslast = dls1_.nst;
    if (itask == 1) {
        goto LABEL_210;
    } else if (itask == 2) {
        goto LABEL_250;
    } else if (itask == 3) {
        goto LABEL_220;
    } else if (itask == 4) {
        goto LABEL_230;
    } else if (itask == 5) {
        goto LABEL_240;
    }
LABEL_210:
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_250;
    DINTDY(tout, 0, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, y, iflag);
    if (iflag != 0) goto LABEL_627;
    t = tout;
    goto LABEL_420;
LABEL_220:
    tp = dls1_.tn - dls1_.hu * (1.0 + 100.0 * dls1_.uround);
    if ((tp  - tout) * dls1_.h > 0.0) goto LABEL_623;
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_250;
    t = dls1_.tn;
    goto LABEL_400;
LABEL_230:
    tcrit = ARRAYF(rwork, 1);
    if ((dls1_.tn - tcrit)  * dls1_.h > 0.0) goto LABEL_624;
    if ((tcrit - tout) * dls1_.h < 0.0) goto LABEL_625;
    if ((dls1_.tn - tout)   * dls1_.h < 0.0) goto LABEL_245;
    DINTDY(tout, 0, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, y, iflag);
    if (iflag != 0) goto LABEL_627;
    t = tout;
    goto LABEL_420;
LABEL_240:
    tcrit = ARRAYF(rwork, 1);
    if ((dls1_.tn - tcrit) * dls1_.h > 0.0) goto LABEL_624;    
LABEL_245:
    hmx  = std::abs(dls1_.tn) + std::abs(dls1_.h);
    ihit = std::abs(dls1_.tn - tcrit) <= (100.0 * dls1_.uround * hmx);
    if (ihit) t = tcrit;
    if (ihit) goto LABEL_400;
    tnext = dls1_.tn + dls1_.h * (1.0 + 4.0 * dls1_.uround);
    if ((tnext - tcrit) * dls1_.h <= 0.0) goto LABEL_250;
    dls1_.h = (tcrit - dls1_.tn) * (1.0 - 4.0 * dls1_.uround);
    if (istate == 2 && dls1_.jstart >= 0) dls1_.jstart = -2;
//-----------------------------------------------------------------------
// Block E.
// The next block is normally executed for all calls and contains
// the call to the one-step core integrator DSTODA.
//
// This is a looping point for the integration steps.
//
// First check for too many steps being taken, update EWT (if not at
// start of problem), check for too much accuracy being requested, and
// check for H below the roundoff level in T.
//-----------------------------------------------------------------------
LABEL_250:
    if (dls1_.meth == dlsa_.mused) goto LABEL_255;
    if (dlsa_.insufr == 1) goto LABEL_550;
    if (dlsa_.insufi == 1) goto LABEL_555;
LABEL_255:
    if ((dls1_.nst - dls1_.nslast) >= dls1_.mxstep) goto LABEL_500;
    DEWSET(dls1_.n, itol, rtol, atol, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lewt));
    for (i = 1; i <= dls1_.n; ++i) {
        if (ARRAYF(rwork, i + dls1_.lewt - 1) <= 0.0) goto LABEL_510;
        ARRAYF(rwork, i + dls1_.lewt - 1) = 1.0 / ARRAYF(rwork, i + dls1_.lewt - 1);
    }
LABEL_270:
    tolsf = dls1_.uround * DMNORM(dls1_.n, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lewt));
    if (tolsf <= 1.0) goto LABEL_280;
    tolsf *= 2.0;
    if (dls1_.nst == 0) goto LABEL_626;
    goto LABEL_520;
LABEL_280:
    if ((dls1_.tn + dls1_.h) != dls1_.tn) goto LABEL_290;
    dls1_.nhnil += 1;
    if (dls1_.nhnil > dls1_.mxhnil) goto LABEL_290;
    msg = "DLSODA-  Warning..internal T (=R1) and H (=R2) are";
    XERRWD(msg, 50, 101, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      such that in the machine, T + H = T on the next step  ";
    XERRWD(msg, 60, 101, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      (H = step size). Solver will continue anyway";
    XERRWD(msg, 50, 101, 0, 0, 0, 0, 2, dls1_.tn, dls1_.h);
    if (dls1_.nhnil < dls1_.mxhnil) goto LABEL_290;
    msg = "DLSODA-  Above warning has been issued I1 times.  ";
    XERRWD(msg, 50, 102, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      It will not be issued again for this problem";
    XERRWD(msg, 50, 102, 0, 1, dls1_.mxhnil, 0, 0, 0.0, 0.0);
LABEL_290:
//-----------------------------------------------------------------------
//   CALL DSTODA(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,DPRJA,DSOLSY)
//-----------------------------------------------------------------------
    DSTODA(neq, y, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lewt), 
        &ARRAYF(rwork, dls1_.lsavf), &ARRAYF(rwork, dls1_.lacor), &ARRAYF(rwork, dls1_.lwm), &ARRAYF(iwork, dls1_.liwm), 
        f, jac, &Odepack::DPRJA, &Odepack::DSOLSY, user_data);
    kgo = 1 - dls1_.kflag;
    if (kgo == 1) {
        goto LABEL_300;
    } else if (kgo == 2) {
        goto LABEL_530;
    } else if (kgo == 3) {
        goto LABEL_540;
    }
//-----------------------------------------------------------------------
// Block F.
// The following block handles the case of a successful return from the
// core integrator (KFLAG = 0).
// If a method switch was just made, record TSW, reset MAXORD,
// set JSTART to -1 to signal DSTODA to complete the switch,
// and do extra printing of data if IXPR = 1.
// Then, in any case, check for stop conditions.
//-----------------------------------------------------------------------
LABEL_300:
    dls1_.init = 1;
    if (dls1_.meth == dlsa_.mused) goto LABEL_310;
    dlsa_.tsw = dls1_.tn;
    dls1_.maxord = dlsa_.mxordn;
    if (dls1_.meth == 2) dls1_.maxord = dlsa_.mxords;
    if (dls1_.meth == 2) ARRAYF(rwork, dls1_.lwm) = std::sqrt(dls1_.uround);
    dlsa_.insufr = std::min(dlsa_.insufr, 1);
    dlsa_.insufi = std::min(dlsa_.insufi, 1);
    dls1_.jstart = -1;
    if (dlsa_.ixpr == 0) goto LABEL_310;
    if (dls1_.meth == 2) {
        msg = "DLSODA- A switch to the BDF (stiff) method has occurred     ";
        XERRWD(msg, 60, 105, 0, 0, 0, 0, 0, 0.0, 0.0);
    }
    if (dls1_.meth == 1) {
        msg = "DLSODA- A switch to the Adams (nonstiff) method has occurred";
        XERRWD(msg, 60, 106, 0, 0, 0, 0, 0, 0.0, 0.0);
    }
    msg = "     at T = R1,  tentative step size H = R2,  step NST = I1 ";
    XERRWD(msg, 60, 107, 0, 1, dls1_.nst, 0, 2, dls1_.tn, dls1_.h);
LABEL_310:
    if (itask == 1) {
        goto LABEL_320;
    } else if (itask == 2) {
        goto LABEL_400;
    } else if (itask == 3) {
        goto LABEL_330;
    } else if (itask == 4) {
        goto LABEL_340;
    } else if (itask == 5) {
        goto LABEL_350;
    }
// ITASK = 1.  If TOUT has been reached, interpolate. -------------------
LABEL_320:
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_250;
    DINTDY(tout, 0, &ARRAYF(rwork, dls1_.lyh) , dls1_.nyh, y, iflag);
    t = tout;
    goto LABEL_420;
// ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
LABEL_330:
    if ((dls1_.tn - tout) * dls1_.h >= 0.0) goto LABEL_400;
    goto LABEL_250;
// ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
LABEL_340:
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_345;
    DINTDY(tout, 0, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, y, iflag);
    t = tout;
    goto LABEL_420;
LABEL_345:
    hmx = std::abs(dls1_.tn) + std::abs(dls1_.h);
    ihit = std::abs(dls1_.tn - tcrit) <= (100.0 * dls1_.uround * hmx);
    if (ihit) goto LABEL_400;
    tnext = dls1_.tn + dls1_.h * (1.0 + 4.0 * dls1_.uround);
    if ((tnext - tcrit) * dls1_.h <= 0.0) goto LABEL_250;
    dls1_.h = (tcrit - dls1_.tn) * (1.0 - 4.0 * dls1_.uround);
    if (dls1_.jstart >= 0) dls1_.jstart = -2;
    goto LABEL_250;
// ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
LABEL_350:
    hmx  = std::abs(dls1_.tn) + std::abs(dls1_.h);
    ihit = std::abs(dls1_.tn - tcrit) <= (100.0 * dls1_.uround * hmx);
//-----------------------------------------------------------------------
// Block G.
// The following block handles all successful returns from DLSODA.
// If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
// ISTATE is set to 2, and the optional outputs are loaded into the
// work arrays before returning.
//-----------------------------------------------------------------------
LABEL_400:
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = ARRAYF(rwork, i + dls1_.lyh - 1);
    }
    t = dls1_.tn;
    if (itask != 4 && itask != 5) goto LABEL_420;
    if (ihit) t = tcrit;
LABEL_420:
    istate = 2;
    ARRAYF(rwork, 11) = dls1_.hu;
    ARRAYF(rwork, 12) = dls1_.h;
    ARRAYF(rwork, 13) = dls1_.tn;
    ARRAYF(rwork, 15) = dlsa_.tsw;
    ARRAYF(iwork, 11) = dls1_.nst;
    ARRAYF(iwork, 12) = dls1_.nfe;
    ARRAYF(iwork, 13) = dls1_.nje;
    ARRAYF(iwork, 14) = dls1_.nqu;
    ARRAYF(iwork, 15) = dls1_.nq;
    ARRAYF(iwork, 19) = dlsa_.mused;
    ARRAYF(iwork, 20) = dls1_.meth;
    return;
//-----------------------------------------------------------------------
// Block H.
// The following block handles all unsuccessful returns other than
// those for illegal input.  First the error message routine is called.
// If there was an error test or convergence test failure, IMXER is set.
// Then Y is loaded from YH and T is set to TN.
// The optional outputs are loaded into the work arrays before returning.
//-----------------------------------------------------------------------
// The maximum number of steps was taken before reaching TOUT. ----------
LABEL_500:
    msg = "DLSODA-  At current T (=R1), MXSTEP (=I1) steps   ";
    XERRWD(msg, 50, 201, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      taken on this call before reaching TOUT     ";
    XERRWD(msg, 50, 201, 0, 1, dls1_.mxstep, 0, 1, dls1_.tn, 0.0);
    istate = -1;
    goto LABEL_580;
// EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
LABEL_510:
    ewti = ARRAYF(rwork, dls1_.lewt+i-1);
    msg = "DLSODA-  At T (=R1), EWT(I1) has become R2 .le. 0.";
    XERRWD(msg, 50, 202, 0, 1, i, 0, 2, dls1_.tn, ewti);
    istate = -6;
    goto LABEL_580;
// Too much accuracy requested for machine precision. -------------------
LABEL_520:
    msg = "DLSODA-  At T (=R1), too much accuracy requested  ";
    XERRWD(msg, 50, 203, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      for precision of machine..  see TOLSF (=R2) ";
    XERRWD(msg, 50, 203, 0, 0, 0, 0, 2, dls1_.tn, tolsf);
    ARRAYF(rwork, 14) = tolsf;
    istate = -2;
    goto LABEL_580;
// KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. ------
LABEL_530:
    msg = "DLSODA-  At T(=R1) and step size H(=R2), the error";
    XERRWD(msg, 50, 204, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      test failed repeatedly or with ABS(H) = HMIN'";
    XERRWD(msg, 50, 204, 0, 0, 0, 0, 2, dls1_.tn, dls1_.h);
    istate = -4;
    goto LABEL_560;
// KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. -----
LABEL_540:
    msg = "DLSODA-  At T (=R1) and step size H (=R2), the    ";
    XERRWD(msg, 50, 205, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      corrector convergence failed repeatedly     ";
    XERRWD(msg, 50, 205, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      or with ABS(H) = HMIN   ";
    XERRWD(msg, 30, 205, 0, 0, 0, 0, 2, dls1_.tn, dls1_.h);
    istate = -5;
    goto LABEL_560;
// RWORK length too small to proceed. -----------------------------------
LABEL_550:
    msg = "DLSODA-  At current T(=R1), RWORK length too small";
    XERRWD(msg, 50, 206, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      to proceed.  The integration was otherwise successful.";
    XERRWD(msg, 60, 207, 0, 0, 0, 0, 1, dls1_.tn, 0.0);
    istate = -7;
    goto LABEL_580;
// IWORK length too small to proceed. -----------------------------------
LABEL_555:
    msg = "DLSODA-  At current T(=R1), IWORK length too small";
    XERRWD(msg, 50, 207, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      to proceed.  The integration was otherwise successful.";
    XERRWD(msg, 60, 207, 0, 0, 0, 0, 1, dls1_.tn, 0.0);
    istate = -7;
    goto LABEL_580;
// Compute IMXER if relevant. -------------------------------------------
LABEL_560:
    big = 0.0;
    imxer = 1;
    for (i = 1; i <= dls1_.n; ++i) {
        size = std::abs(ARRAYF(rwork, i + dls1_.lacor - 1) * ARRAYF(rwork, i + dls1_.lewt - 1));
        if (big >= size) continue;
        big = size;
        imxer = i;
    }
    ARRAYF(iwork, 16) = imxer;
// Set Y vector, T, and optional outputs. -------------------------------
LABEL_580:
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = ARRAYF(rwork, i + dls1_.lyh - 1);
    }
    t = dls1_.tn;
    ARRAYF(rwork, 11) = dls1_.hu;
    ARRAYF(rwork, 12) = dls1_.h;
    ARRAYF(rwork, 13) = dls1_.tn;
    ARRAYF(rwork, 15) = dlsa_.tsw;
    ARRAYF(iwork, 11) = dls1_.nst;
    ARRAYF(iwork, 12) = dls1_.nfe;
    ARRAYF(iwork, 13) = dls1_.nje;
    ARRAYF(iwork, 14) = dls1_.nqu;
    ARRAYF(iwork, 15) = dls1_.nq;
    ARRAYF(iwork, 19) = dlsa_.mused;
    ARRAYF(iwork, 20) = dls1_.meth;
    return;
//-----------------------------------------------------------------------
// Block I.
// The following block handles all error returns due to illegal input
// (ISTATE = -3), as detected before calling the core integrator.
// First the error message routine is called.  If the illegal input 
// is a negative ISTATE, the run is aborted (apparent infinite loop).
//-----------------------------------------------------------------------
LABEL_601:
    msg = "DLSODA-  ISTATE (=I1) illegal.";
    XERRWD(msg, 30, 1, 0, 1, istate, 0, 0, 0.0, 0.0);
    if (istate < 0) goto LABEL_800;
    goto LABEL_700;
LABEL_602:
    msg = "DLSODA-  ITASK (=I1) illegal. ";
    XERRWD(msg, 30, 2, 0, 1, itask, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_603:
    msg = "DLSODA-  ISTATE .gt. 1 but DLSODA not initialized.";
    XERRWD(msg, 50, 3, 0, 0, 0, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_604:
    msg = "DLSODA-  NEQ (=I1) .lt. 1     ";
    XERRWD(msg, 30, 4, 0, 1, neq, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_605:
    msg = "DLSODA-  ISTATE = 3 and NEQ increased (I1 to I2). ";
    XERRWD(msg, 50, 5, 0, 2, dls1_.n, neq, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_606:
    msg = "DLSODA-  ITOL (=I1) illegal.  ";
    XERRWD(msg, 30, 6, 0, 1, itol, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_607:
    msg = "DLSODA-  IOPT (=I1) illegal   ";
    XERRWD(msg, 30, 7, 0, 1, iopt, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_608:
    msg = "DLSODA-  JT (=I1) illegal.    ";
    XERRWD(msg, 30, 8, 0, 1, jt, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_609:
    msg = "DLSODA-  ML (=I1) illegal: .lt.0 or .ge.NEQ (=I2) ";
    XERRWD(msg, 50, 9, 0, 2, ml, neq, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_610:
    msg = "DLSODA-  MU (=I1) illegal: .lt.0 or .ge.NEQ (=I2) ";
    XERRWD(msg, 50, 10, 0, 2, mu, neq, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_611:
    msg = "DLSODA-  IXPR (=I1) illegal.  ";
    XERRWD(msg, 30, 11, 0, 1, dlsa_.ixpr, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_612:
    msg = "DLSODA-  MXSTEP (=I1) .lt. 0  ";
    XERRWD(msg, 30, 12, 0, 1, dls1_.mxstep, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_613:
    msg = "DLSODA-  MXHNIL (=I1) .lt. 0  ";
    XERRWD(msg, 30, 13, 0, 1, dls1_.mxhnil, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_614:
    msg = "DLSODA-  TOUT (=R1) behind T (=R2)      ";
    XERRWD(msg, 40, 14, 0, 0, 0, 0, 2, tout, t);
    msg = "      Integration direction is given by H0 (=R1)  ";
    XERRWD(msg, 50, 14, 0, 0, 0, 0, 1, h0, 0.0);
    goto LABEL_700;
LABEL_615:
    msg = "DLSODA-  HMAX (=R1) .lt. 0.0  ";
    XERRWD(msg, 30, 15, 0, 0, 0, 0, 1, hmax, 0.0);
    goto LABEL_700;
LABEL_616:
    msg = "DLSODA-  HMIN (=R1) .lt. 0.0  ";
    XERRWD(msg, 30, 16, 0, 0, 0, 0, 1, dls1_.hmin, 0.0);
    goto LABEL_700;
LABEL_617:
    msg = "DLSODA-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)";
    XERRWD(msg, 60, 17, 0, 2, lenrw, lrw, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_618:
    msg = "DLSODA-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)";
    XERRWD(msg, 60, 18, 0, 2, leniw, liw, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_619:
    msg = "DLSODA-  RTOL(I1) is R1 .lt. 0.0        ";
    XERRWD(msg, 40, 19, 0, 1, i, 0, 1, rtoli, 0.0);
    goto LABEL_700;
LABEL_620:
    msg = "DLSODA-  ATOL(I1) is R1 .lt. 0.0        ";
    XERRWD(msg, 40, 20, 0, 1, i, 0, 1, atoli, 0.0);
    goto LABEL_700;
LABEL_621:
    ewti = ARRAYF(rwork, dls1_.lewt+i-1);
    msg = "DLSODA-  EWT(I1) is R1 .le. 0.0         ";
    XERRWD(msg, 40, 21, 0, 1, i, 0, 1, ewti, 0.0);
    goto LABEL_700;
LABEL_622:
    msg = "DLSODA-  TOUT(=R1) too close to T(=R2) to start integration.";
    XERRWD(msg, 60, 22, 0, 0, 0, 0, 2, tout, t);
    goto LABEL_700;
LABEL_623:
    msg = "DLSODA-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  ";
    XERRWD(msg, 60, 23, 0, 1, itask, 0, 2, tout, tp);
    goto LABEL_700;
LABEL_624:
    msg = "DLSODA-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   ";
    XERRWD(msg, 60, 24, 0, 0, 0, 0, 2, tcrit, dls1_.tn);
    goto LABEL_700;
LABEL_625:
    msg = "DLSODA-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   ";
    XERRWD(msg, 60, 25, 0, 0, 0, 0, 2, tcrit, tout);
    goto LABEL_700;
LABEL_626:
    msg = "DLSODA-  At start of problem, too much accuracy   ";
    XERRWD(msg, 50, 26, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      requested for precision of machine..  See TOLSF (=R1) ";
    XERRWD(msg, 60, 26, 0, 0, 0, 0, 1, tolsf, 0.0);
    goto LABEL_700;
LABEL_627:
    msg = "DLSODA-  Trouble in DINTDY.  ITASK = I1, TOUT = R1";
    XERRWD(msg, 50, 27, 0, 1, itask, 0, 1, tout, 0.0);
    goto LABEL_700;
LABEL_628:
    msg = "DLSODA-  MXORDN (=I1) .lt. 0  ";
    XERRWD(msg, 30, 28, 0, 1, dlsa_.mxordn, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_629:
    msg = "DLSODA-  mxords (=I1) .lt. 0  ";
    XERRWD(msg, 30, 29, 0, 1, dlsa_.mxords, 0, 0, 0.0, 0.0);
//
LABEL_700:
    istate = -3;
    return;
//
LABEL_800:
    msg = "DLSODA-  Run aborted.. apparent infinite loop.    ";
    XERRWD(msg, 50, 303, 2, 0, 0, 0, 0, 0.0, 0.0);
    return;
}

void Odepack::DLSODAR(ODEPACK_FUNCTION f, const int neq, double *y, double &t, const double tout,
    const int itol, double *rtol, double *atol, const int itask, int &istate,
    const int iopt, double *rwork, const int lrw, int *iwork, const int liw,
    ODEPACK_JACOBIAN1 jac, const int jt, ODEPACK_CONSTRAINT g, const int ng, int *jroot,
    void *user_data)
{
    int i, i1, i2, iflag, imxer, kgo, leniw, lenrw, lenwm, lf0, ml, mu;
    int len1, len1c, len1n, len1s, len2, leniwc, lenrwc;
    int irfp, irt, lenyh, lyhnew;
    double atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli, tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0;
    bool ihit;
    const int mord[2] = {12, 5};
    const int mxstp0  = 500;
    const int mxhnl0  = 10;
    std::string msg;
// -----------------------------------------------------------------------
//  Block A.
//  This code block is executed on every call.
//  It tests ISTATE and ITASK for legality and branches appropriately.
//  If ISTATE > 1 but the flag INIT shows that initialization has
//  not yet been done, an error return occurs.
//  If ISTATE = 1 and TOUT = T, return immediately.
// -----------------------------------------------------------------------
    if (istate < 1 || istate > 3) goto LABEL_601;
    if (itask  < 1 || itask  > 5) goto LABEL_602;
    dlsr_.itaskc = itask;
    if (istate == 1)     goto LABEL_10;
    if (dls1_.init == 0) goto LABEL_603;
    if (istate == 2)     goto LABEL_200;
    goto LABEL_20;
LABEL_10:
    dls1_.init = 0;
    if (tout == t) return;
// -----------------------------------------------------------------------
//  Block B.
//  The next code block is executed for the initial (ISTATE = 1),
//  or for a continuation with parameter changes (ISTATE = 3).
//  It contains checking of all inputs and various initializations.
//
//  First check legality of the non-optional inputs NEQ, ITOL, IOPT,
//  JT, ML, MU, and NG.
// -----------------------------------------------------------------------
LABEL_20:
    if (neq <= 0)      goto LABEL_604;
    if (istate == 1)   goto LABEL_25;
    if (neq > dls1_.n) goto LABEL_605;
LABEL_25:
    dls1_.n = neq;
    if (itol < 1 || itol > 4) goto LABEL_606;
    if (iopt < 0 || iopt > 1) goto LABEL_607;
    if (jt == 3 || jt < 1 || jt > 5) goto LABEL_608;
    dlsa_.jtyp = jt;
    if (jt <= 2) goto LABEL_30;
    ml = ARRAYF(iwork, 1);
    mu = ARRAYF(iwork, 2);
    if (ml < 0 || ml >= dls1_.n) goto LABEL_609;
    if (mu < 0 || mu >= dls1_.n) goto LABEL_610;
LABEL_30:
    if (ng < 0) goto LABEL_630;
    if (istate == 1) goto LABEL_35;
    if (dlsr_.irfnd == 0 && ng != dlsr_.ngc) goto LABEL_631;
LABEL_35:
    dlsr_.ngc = ng;
// Next process and check the optional inputs. --------------------------
    if (iopt == 1) goto LABEL_40;
    dlsa_.ixpr   = 0;
    dls1_.mxstep = mxstp0;
    dls1_.mxhnil = mxhnl0;
    dls1_.hmxi   = 0.0;
    dls1_.hmin   = 0.0;
    if (istate != 1) goto LABEL_60;
    h0           = 0.0;
    dlsa_.mxordn = ARRAYF(mord, 1);
    dlsa_.mxords = ARRAYF(mord, 2);
    goto LABEL_60;
LABEL_40:
    dlsa_.ixpr   = ARRAYF(iwork, 5);
    if (dlsa_.ixpr < 0 || dlsa_.ixpr > 1) goto LABEL_611;
    dls1_.mxstep = ARRAYF(iwork, 6);
    if (dls1_.mxstep  < 0) goto LABEL_612;
    if (dls1_.mxstep == 0) dls1_.mxstep = mxstp0;
    dls1_.mxhnil = ARRAYF(iwork, 7);
    if (dls1_.mxhnil  < 0) goto LABEL_613;
    if (dls1_.mxhnil == 0) dls1_.mxhnil = mxhnl0;
    if (istate != 1) goto LABEL_50;
    h0           = ARRAYF(rwork, 5);
    dlsa_.mxordn = ARRAYF(iwork, 8);
    if (dlsa_.mxordn  < 0) goto LABEL_628;
    if (dlsa_.mxordn == 0) dlsa_.mxordn = 100;
    dlsa_.mxordn = std::min(dlsa_.mxordn, ARRAYF(mord, 1));
    dlsa_.mxords = ARRAYF(iwork, 9);
    if (dlsa_.mxords  < 0) goto LABEL_629;
    if (dlsa_.mxords == 0) dlsa_.mxords = 100;
    dlsa_.mxords = std::min(dlsa_.mxords, ARRAYF(mord, 2));
    if ((tout - t) * h0 < 0.0) goto LABEL_614;
LABEL_50:
    hmax       = ARRAYF(rwork, 6);
    if (hmax < 0.0) goto LABEL_615;
    dls1_.hmxi = 0.0;
    if (hmax > 0.0) dls1_.hmxi = 1.0 / hmax;
    dls1_.hmin = ARRAYF(rwork, 7);
    if (dls1_.hmin < 0.0) goto LABEL_616;
//-----------------------------------------------------------------------
// Set work array pointers and check lengths LRW and LIW.
// If ISTATE = 1, METH is initialized to 1 here to facilitate the
// checking of work space lengths.
// Pointers to segments of RWORK and IWORK are named by prefixing L to
// the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
// Segments of RWORK (in order) are denoted  G0, G1, GX, YH, WM,
// EWT, SAVF, ACOR.
// If the lengths provided are insufficient for the current method,
// an error return occurs.  This is treated as illegal input on the
// first call, but as a problem interruption with ISTATE = -7 on a
// continuation call.  If the lengths are sufficient for the current
// method but not for both methods, a warning message is sent.
// -----------------------------------------------------------------------
LABEL_60:
    if (istate == 1) dls1_.meth = 1;
    if (istate == 1) dls1_.nyh = dls1_.n;
    dlsr_.lg0 = 21;
    dlsr_.lg1 = dlsr_.lg0 + ng;
    dlsr_.lgx = dlsr_.lg1 + ng;
    lyhnew    = dlsr_.lgx + ng;
    if (istate == 1) dls1_.lyh = lyhnew;
    if (lyhnew == dls1_.lyh) goto LABEL_62;
// If ISTATE = 3 and NG was changed, shift YH to its new location. ------
    lenyh = dls1_.l * dls1_.nyh;
    if (lrw < (lyhnew - 1 + lenyh)) goto LABEL_62;
    i1 = 1;
    if (lyhnew > dls1_.lyh) i1 = -1;
    DCOPY(lenyh, &ARRAYF(rwork, dls1_.lyh), i1, &ARRAYF(rwork, lyhnew), i1);
    dls1_.lyh = lyhnew;
LABEL_62:
    len1n     = lyhnew - 1 + (dlsa_.mxordn + 1) * dls1_.nyh;
    len1s     = lyhnew - 1 + (dlsa_.mxords + 1) * dls1_.nyh;
    dls1_.lwm = len1s + 1;
    if (jt <= 2) lenwm = dls1_.n * dls1_.n + 2;
    if (jt >= 4) lenwm = (2 * ml + mu + 1) * dls1_.n + 2;
    len1s = len1s + lenwm;
    len1c = len1n;
    if (dls1_.meth == 2) len1c = len1s;
    len1              = std::max(len1n, len1s);
    len2              = 3 * dls1_.n;
    lenrw             = len1 + len2;
    lenrwc            = len1c + len2;
    ARRAYF(iwork, 17) = lenrw;
    dls1_.liwm        = 1;
    leniw             = 20 + dls1_.n;
    leniwc            = 20;
    if (dls1_.meth == 2) leniwc = leniw;
    ARRAYF(iwork, 18) = leniw;
    if (istate == 1 && lrw < lenrwc) goto LABEL_617;
    if (istate == 1 && liw < leniwc) goto LABEL_618;
    if (istate == 3 && lrw < lenrwc) goto LABEL_550;
    if (istate == 3 && liw < leniwc) goto LABEL_555;
    dls1_.lewt   = len1 + 1;
    dlsa_.insufr = 0;
    if (lrw >= lenrw) goto LABEL_65;
    dlsa_.insufr = 2;
    dls1_.lewt   = len1c + 1;
    msg = "DLSODAR-  Warning.. RWORK length is sufficient for now, but ";
    XERRWD(msg, 60, 103, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      may not be later.  Integration will proceed anyway.   ";
    XERRWD(msg, 60, 103, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      Length needed is LENRW = I1, while LRW = I2.";
    XERRWD(msg, 50, 103, 0, 2, lenrw, lrw, 0, 0.0, 0.0);
LABEL_65:
    dls1_.lsavf  = dls1_.lewt  + dls1_.n;
    dls1_.lacor  = dls1_.lsavf + dls1_.n;
    dlsa_.insufi = 0;
    if (liw >= leniw) goto LABEL_70;
    dlsa_.insufi = 2;
    msg = "DLSODAR-  Warning.. IWORK length is sufficient for now, but ";
    XERRWD(msg, 60, 104, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      may not be later.  Integration will proceed anyway.   ";
    XERRWD(msg, 60, 104, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      Length needed is LENIW = I1, while LIW = I2.";
    XERRWD(msg, 50, 104, 0, 2, leniw, liw, 0, 0.0, 0.0);
LABEL_70:
// Check RTOL and ATOL for legality. ------------------------------------
    rtoli = ARRAYF(rtol, 1);
    atoli = ARRAYF(atol, 1);
    for (i = 1; i <= dls1_.n; ++i) {
        if (itol >= 3) rtoli = ARRAYF(rtol, i);
        if (itol == 2 || itol == 4) atoli = ARRAYF(atol, i);
        if (rtoli < 0.0) goto LABEL_619;
        if (atoli < 0.0) goto LABEL_620;
    }
    if (istate == 1) goto LABEL_100;
// if ISTATE = 3, set flag to signal parameter changes to DSTODA. -------
    dls1_.jstart = -1;
    if (dls1_.n == dls1_.nyh) goto LABEL_200;
// NEQ was reduced.  zero part of yh to avoid undefined references. -----
    i1 = dls1_.lyh + dls1_.l * dls1_.nyh;
    i2 = dls1_.lyh + (dls1_.maxord + 1) * dls1_.nyh - 1;
    if (i1 > i2) goto LABEL_200;
    for (i = i1; i1 <= i2; ++i) {
        ARRAYF(rwork, i) = 0.0;
    }
    goto LABEL_200;
// -----------------------------------------------------------------------
//  Block C.
//  The next block is for the initial only (ISTATE = 1).
//  It contains all remaining initializations, the initial to F,
//  and the calculation of the initial step size.
//  The error weights in EWT are inverted after being loaded.
// -----------------------------------------------------------------------
LABEL_100:
    dls1_.uround = DUMACH();
    dls1_.tn     = t;
    dlsa_.tsw    = t;
    dls1_.maxord = dlsa_.mxordn;
    if (itask != 4 && itask != 5) goto LABEL_110;
    tcrit        = ARRAYF(rwork, 1);
    if ((tcrit - tout) * (tout - t) < 0.0) goto LABEL_625;
    if (h0 != 0.0 && (t + h0 - tcrit) * h0 > 0.0) h0 = tcrit - t;
LABEL_110:
    dls1_.jstart = 0;
    dls1_.nhnil  = 0;
    dls1_.nst    = 0;
    dls1_.nje    = 0;
    dls1_.nslast = 0;
    dls1_.hu     = 0.0;
    dls1_.nqu    = 0;
    dlsa_.mused  = 0;
    dls1_.miter  = 0;
    dls1_.ccmax  = 0.3;
    dls1_.maxcor = 3;
    dls1_.msbp   = 20;
    dls1_.mxncf  = 10;
// Initial call to F.  (LF0 points to YH(*,2).) -------------------------
    lf0 = dls1_.lyh + dls1_.nyh;
    (*f)(neq, t, y, &ARRAYF(rwork, lf0), user_data);
    dls1_.nfe = 1;
// Load the initial value vector in YH. ---------------------------------
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(rwork, i + dls1_.lyh - 1) = ARRAYF(y, i);
    }
// Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
    dls1_.nq = 1;
    dls1_.h  = 1.0;
    DEWSET(dls1_.n, itol, rtol, atol, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lewt));
    for (i = 1; i <= dls1_.n; ++i) {
        if (ARRAYF(rwork, i + dls1_.lewt - 1) <= 0.0) goto LABEL_621;
        ARRAYF(rwork, i + dls1_.lewt - 1) = 1.0 / ARRAYF(rwork, i + dls1_.lewt - 1);
    }
//-----------------------------------------------------------------------
// The coding below computes the step size, H0, to be attempted on the
// first step, unless the user has supplied a value for this.
// First check that TOUT - T differs significantly from zero.
// A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
// if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
// so as to be between 100*UROUND and 1.0E-3.
// Then the computed value H0 is given by:
//
//   H0**(-2)  =  1./(TOL * w0**2)  +  TOL * (norm(F))**2
//
// where   w0     = MAX ( ABS(T), ABS(TOUT) ),
//         F      = the initial value of the vector f(t,y), and
//         norm() = the weighted vector norm used throughout, given by
//                 the DMNORM function routine, and weighted by the
//                 tolerances initially loaded into the EWT array.
// The sign of H0 is inferred from the initial values of TOUT and T.
// ABS(H0) is made .le. ABS(TOUT-T) in any case.
// -----------------------------------------------------------------------
    if (h0 != 0.0) goto LABEL_180;
    tdist = std::abs(tout - t);
    w0    = std::max(std::abs(t), std::abs(tout));
    if (tdist < 2.0 * dls1_.uround * w0) goto LABEL_622;
    tol   = ARRAYF(rtol, 1);
    if (itol <= 2) goto LABEL_140;
    for (i = 1; i <= dls1_.n; ++i) {
        tol = std::max(tol, ARRAYF(rtol, i));
    }
LABEL_140:
    if (tol > 0.0) goto LABEL_160;
    atoli = ARRAYF(atol, 1);
    for (i = 1; i <= dls1_.n; ++i) {
        if (itol == 2 || itol == 4) atoli = ARRAYF(atol, i);
        ayi = std::abs(ARRAYF(y, i));
        if (ayi != 0.0) tol = std::max(tol, atoli / ayi);
    }
LABEL_160:
    tol = std::max(tol, 100.0 * dls1_.uround);
    tol = std::min(tol, 0.001);
    sum = DMNORM(dls1_.n, &ARRAYF(rwork, lf0), &ARRAYF(rwork, dls1_.lewt));
    sum = 1.0 / (tol * w0 * w0) + tol * sum * sum;
    h0  = 1.0 / std::sqrt(sum);
    h0  = std::min(h0, tdist);
    h0  = std::abs(h0) * ((tout - t) >= 0 ? 1: -1);
// Adjust H0 if necessary to meet HMAX bound. ---------------------------
LABEL_180:
    rh = std::abs(h0) * dls1_.hmxi;
    if (rh > 1.0) h0 = h0 / rh;
// Load H with H0 and scale YH(*,2) by H0. ------------------------------
    dls1_.h = h0;
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(rwork, i + lf0 - 1) = h0 * ARRAYF(rwork, i + lf0 - 1);
    }
//
// Check for a zero of g at T. ------------------------------------------
    dlsr_.irfnd = 0;
    dlsr_.toutc = tout;
    if (dlsr_.ngc == 0) goto LABEL_270;
    DRCHEK(1, g, neq, y, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh,
        &ARRAYF(rwork, dlsr_.lg0), &ARRAYF(rwork, dlsr_.lg1), &ARRAYF(rwork, dlsr_.lgx), jroot, irt, user_data);
    if (irt == 0) goto LABEL_270;
    goto LABEL_632;
// -----------------------------------------------------------------------
//  Block D.
//  The next code block is for continuation calls only (ISTATE = 2 or 3)
//  and is to check stop conditions before taking a step.
//  First, DRCHEK is called to check for a root within the last step
//  taken, other than the last root found there, if any.
//  If ITASK = 2 or 5, and y(TN) has not yet been returned to the user
//  because of an intervening root, return through Block G.
// -----------------------------------------------------------------------
LABEL_200:
    dls1_.nslast = dls1_.nst;
//
    irfp         = dlsr_.irfnd;
    if (dlsr_.ngc == 0) goto LABEL_205;
    if (itask == 1 || itask == 4) dlsr_.toutc = tout;
    DRCHEK(2, g, neq, y, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, 
        &ARRAYF(rwork, dlsr_.lg0), &ARRAYF(rwork, dlsr_.lg1), &ARRAYF(rwork, dlsr_.lgx), jroot, irt, user_data);
    if (irt != 1) goto LABEL_205;
    dlsr_.irfnd = 1;
    istate      = 3;
    t           = dlsr_.t0;
    goto LABEL_425;
LABEL_205:
    dlsr_.irfnd = 0;
    if (irfp == 1 && dlsr_.tlast != dls1_.tn && itask == 2) goto LABEL_400;
//
    if (itask == 1) {
        goto LABEL_210;
    } else if (itask == 2) {
        goto LABEL_250;
    } else if (itask == 3) {
        goto LABEL_220;
    } else if (itask == 4) {
        goto LABEL_230;
    } else if (itask == 5) {
        goto LABEL_240;
    }
LABEL_210:
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_250;
    DINTDY(tout, 0, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, y, iflag);
    if (iflag != 0) goto LABEL_627;
    t = tout;
    goto LABEL_420;
LABEL_220:
    tp = dls1_.tn - dls1_.hu * (1.0 + 100.0 * dls1_.uround);
    if ((tp - tout) * dls1_.h > 0.0) goto LABEL_623;
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_250;
    t = dls1_.tn;
    goto LABEL_400;
LABEL_230:
    tcrit = ARRAYF(rwork, 1);
    if ((dls1_.tn - tcrit) * dls1_.h > 0.0) goto LABEL_624;
    if ((tcrit - tout) * dls1_.h < 0.0) goto LABEL_625;
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_245;
    DINTDY(tout, 0, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, y, iflag);
    if (iflag != 0) goto LABEL_627;
    t = tout;
    goto LABEL_420;
LABEL_240:
    tcrit = ARRAYF(rwork, 1);
    if ((dls1_.tn - tcrit) * dls1_.h > 0.0) goto LABEL_624;
LABEL_245:
    hmx  = std::abs(dls1_.tn) + std::abs(dls1_.h);
    ihit = (std::abs(dls1_.tn - tcrit) <= 100.0 * dls1_.uround * hmx);
    if (ihit) t = tcrit;
    if (irfp == 1 && dlsr_.tlast != dls1_.tn && itask == 5) goto LABEL_400;
    if (ihit) goto LABEL_400;
    tnext = dls1_.tn + dls1_.h * (1.0 + 4.0 * dls1_.uround);
    if ((tnext - tcrit) * dls1_.h <= 0.0) goto LABEL_250;
    dls1_.h = (tcrit - dls1_.tn) * (1.0 - 4.0 * dls1_.uround);
    if (istate == 2 && dls1_.jstart >= 0) dls1_.jstart = -2;
// -----------------------------------------------------------------------
//  Block E.
//  The next block is normally executed for all calls and contains
//  the to the one-step core integrator DSTODA.
//
//  This is a looping point for the integration steps.
//
//  First check for too many steps being taken, update EWT (if not at
//  start of problem), check for too much accuracy being requested, and
//  check for H below the roundoff level in T.
// -----------------------------------------------------------------------
LABEL_250:
    if (dls1_.meth == dlsa_.mused) goto LABEL_255;
    if (dlsa_.insufr == 1) goto LABEL_550;
    if (dlsa_.insufi == 1) goto LABEL_555;
LABEL_255:
    if ((dls1_.nst - dls1_.nslast) >= dls1_.mxstep) goto LABEL_500;
    DEWSET(dls1_.n, itol, rtol, atol, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lewt));
    for (i = 1; i <= dls1_.n; ++i) {
        if (ARRAYF(rwork, i + dls1_.lewt - 1) <= 0.0) goto LABEL_510;
        ARRAYF(rwork, i + dls1_.lewt - 1) = 1.0 / ARRAYF(rwork, i + dls1_.lewt - 1);
    }
LABEL_270:
    tolsf = dls1_.uround * DMNORM(dls1_.n, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lewt));
    if (tolsf <= 1.0) goto LABEL_280;
    tolsf = tolsf * 2.0;
    if (dls1_.nst == 0) goto LABEL_626;
    goto LABEL_520;
LABEL_280:
    if ((dls1_.tn + dls1_.h) != dls1_.tn) goto LABEL_290;
    dls1_.nhnil = dls1_.nhnil + 1;
    if (dls1_.nhnil > dls1_.mxhnil) goto LABEL_290;
    msg = "DLSODAR-  Warning..Internal T(=R1) and H(=R2) are ";
    XERRWD(msg, 50, 101, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      such that in the machine, T + H = T on the next step  ";
    XERRWD (msg, 60, 101, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "     (H = step size). Solver will continue anyway.";
    XERRWD(msg, 50, 101, 0, 0, 0, 0, 2, dls1_.tn, dls1_.h);
    if (dls1_.nhnil < dls1_.mxhnil) goto LABEL_290;
    msg = "DLSODAR-  Above warning has been issued I1 times. ";
    XERRWD(msg, 50, 102, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "     It will not be issued again for this problem.";
    XERRWD(msg, 50, 102, 0, 1, dls1_.mxhnil, 0, 0, 0.0, 0.0);
LABEL_290:
//-----------------------------------------------------------------------
//   DSTODA(NEQ,Y,YH,dls1_.nyh,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,DPRJA,DSOLSY)
//-----------------------------------------------------------------------
    DSTODA(neq, y, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, &ARRAYF(rwork, dls1_.lyh), &ARRAYF(rwork, dls1_.lewt),
        &ARRAYF(rwork, dls1_.lsavf), &ARRAYF(rwork, dls1_.lacor), &ARRAYF(rwork, dls1_.lwm), &ARRAYF(iwork, dls1_.liwm),
        f, jac, &Odepack::DPRJA, &Odepack::DSOLSY, user_data);
    kgo = 1 - dls1_.kflag;
    if (kgo == 1) {
        goto LABEL_300;
    } else if (kgo == 2) {
        goto LABEL_530;
    } else if (kgo == 3) {
        goto LABEL_540;
    }
// ---------------------------------------------------------------------
// Block F.
// The following block handles the case of a successful return from the
// core integrator (KFLAG = 0).
// If a method switch was just made, record TSW, reset MAXORD,
// set JSTART to -1 to signal DSTODA to complete the switch,
// and do extra printing of data if IXPR = 1.
// Then call DRCHEK to check for a root within the last step.
// Then, if no root was found, check for stop conditions.
// -----------------------------------------------------------------------
LABEL_300:
    dls1_.init = 1;
    if (dls1_.meth == dlsa_.mused) goto LABEL_310;
    dlsa_.tsw = dls1_.tn;
    dls1_.maxord = dlsa_.mxordn;
    if (dls1_.meth == 2) dls1_.maxord = dlsa_.mxords;
    if (dls1_.meth == 2) ARRAYF(rwork, dls1_.lwm) = std::sqrt(dls1_.uround);
    dlsa_.insufr = std::min(dlsa_.insufr, 1);
    dlsa_.insufi = std::min(dlsa_.insufi, 1);
    dls1_.jstart = -1;
    if (dlsa_.ixpr == 0) goto LABEL_310;
    if (dls1_.meth == 2) {
        msg = "DLSODAR- A switch to the BDF (stiff) method has occurred    ";
        XERRWD(msg, 60, 105, 0, 0, 0, 0, 0, 0.0, 0.0);
    }
    if (dls1_.meth == 1) {
        msg = "DLSODAR- A switch to the Adams (nonstiff) method occurred   ";
        XERRWD(msg, 60, 106, 0, 0, 0, 0, 0, 0.0, 0.0);
    }
    msg = "     at T = R1,  tentative step size H = R2,  step NST = I1 ";
    XERRWD(msg, 60, 107, 0, 1, dls1_.nst, 0, 2, dls1_.tn, dls1_.h);
LABEL_310:
//
    if (dlsr_.ngc == 0) goto LABEL_315;
    DRCHEK(3, g, neq, y, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh,
        &ARRAYF(rwork, dlsr_.lg0), &ARRAYF(rwork, dlsr_.lg1), &ARRAYF(rwork, dlsr_.lgx), jroot, irt, user_data);
    if (irt != 1) goto LABEL_315;
    dlsr_.irfnd = 1;
    istate      = 3;
    t           = dlsr_.t0;
    goto LABEL_425;
LABEL_315:
//
    if (itask == 1) {
        goto LABEL_320;
    } else if (itask == 2) {
        goto LABEL_400;
    } else if (itask == 3) {
        goto LABEL_330;
    } else if (itask == 4) {
        goto LABEL_340;
    } else if (itask == 5) {
        goto LABEL_350;
    }
// ITASK = 1.  If TOUT has been reached, interpolate. -------------------
LABEL_320:
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_250;
    DINTDY(tout, 0, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, y, iflag);
    t = tout;
    goto LABEL_420;
// ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
LABEL_330:
    if ((dls1_.tn - tout) * dls1_.h >= 0.0) goto LABEL_400;
    goto LABEL_250;
// ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
LABEL_340:
    if ((dls1_.tn - tout) * dls1_.h < 0.0) goto LABEL_345;
    DINTDY(tout, 0, &ARRAYF(rwork, dls1_.lyh), dls1_.nyh, y, iflag);
    t = tout;
    goto LABEL_420;
LABEL_345:
    hmx  = std::abs(dls1_.tn) + std::abs(dls1_.h);
    ihit = std::abs(dls1_.tn - tcrit) <= 100.0 * dls1_.uround * hmx;
    if (ihit) goto LABEL_400;
    tnext   = dls1_.tn + dls1_.h * (1.0 + 4.0 * dls1_.uround);
    if ((tnext - tcrit) * dls1_.h <= 0.0) goto LABEL_250;
    dls1_.h = (tcrit - dls1_.tn) * (1.0 - 4.0 * dls1_.uround);
    if (dls1_.jstart >= 0) dls1_.jstart = -2;
    goto LABEL_250;
// ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
LABEL_350:
    hmx  = std::abs(dls1_.tn) + std::abs(dls1_.h);
    ihit = (std::abs(dls1_.tn - tcrit) <= 100.0 * dls1_.uround * hmx);
// -----------------------------------------------------------------------
//  Block G.
//  The following block handles all successful returns from DLSODAR.
//  If ITASK != 1, Y is loaded from YH and T is set accordingly.
//  ISTATE is set to 2, and the optional outputs are loaded into the
//  work arrays before returning.
// -----------------------------------------------------------------------
LABEL_400:
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = ARRAYF(rwork, i + dls1_.lyh - 1);
    }
    t = dls1_.tn;
    if (itask != 4 && itask != 5) goto LABEL_420;
    if (ihit) t = tcrit;
LABEL_420:
    istate = 2;
LABEL_425:
    ARRAYF(rwork, 11) = dls1_.hu;
    ARRAYF(rwork, 12) = dls1_.h;
    ARRAYF(rwork, 13) = dls1_.tn;
    ARRAYF(rwork, 15) = dlsa_.tsw;
    ARRAYF(iwork, 11) = dls1_.nst;
    ARRAYF(iwork, 12) = dls1_.nfe;
    ARRAYF(iwork, 13) = dls1_.nje;
    ARRAYF(iwork, 14) = dls1_.nqu;
    ARRAYF(iwork, 15) = dls1_.nq;
    ARRAYF(iwork, 19) = dlsa_.mused;
    ARRAYF(iwork, 20) = dls1_.meth;
    ARRAYF(iwork, 10) = dlsr_.nge;
    dlsr_.tlast       = t;
    return;
// -----------------------------------------------------------------------
//  Block H.
//  The following block handles all unsuccessful returns other than
//  those for illegal input.  First the error message routine is called.
//  If there was an error test or convergence test failure, IMXER is set.
//  Then Y is loaded from YH and T is set to dls1_.tn.
//  The optional outputs are loaded into the work arrays before returning.
// -----------------------------------------------------------------------
//  The maximum number of steps was taken before reaching TOUT. ----------
LABEL_500:
    msg = "DLSODAR-  At current T (=R1), MXSTEP (=I1) steps  ";
    XERRWD(msg, 50, 201, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      taken on this call before reaching TOUT     ";
    XERRWD(msg, 50, 201, 0, 1, dls1_.mxstep, 0, 1, dls1_.tn, 0.0);
    istate = -1;
    goto LABEL_580;
// EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
LABEL_510:
    ewti = ARRAYF(rwork, dls1_.lewt + i - 1);
    msg = "DLSODAR-  At t(=r1), ewt(i1) has become r2 <= 0.";
    XERRWD(msg, 50, 202, 0, 1, i, 0, 2, dls1_.tn, ewti);
    istate = -6;
    goto LABEL_580;
// Too much accuracy requested for machine precision. -------------------
LABEL_520:
    msg = "DLSODAR-  At T (=R1), too much accuracy requested ";
    XERRWD(msg, 50, 203, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      for precision of machine..  see TOLSF (=R2) ";
    XERRWD(msg, 50, 203, 0, 0, 0, 0, 2, dls1_.tn, tolsf);
    ARRAYF(rwork, 14) = tolsf;
    istate = -2;
    goto LABEL_580;
// KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
LABEL_530:
    msg = "DLSODAR-  At T(=R1), step size H(=R2), the error  ";
    XERRWD (msg, 50, 204, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      test failed repeatedly or with ABS(H) = HMIN";
    XERRWD (msg, 50, 204, 0, 0, 0, 0, 2, dls1_.tn, dls1_.h);
    istate = -4;
    goto LABEL_560;
// KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
LABEL_540:
    msg = "DLSODAR-  At T (=R1) and step size H (=R2), the   ";
    XERRWD(msg, 50, 205, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      corrector convergence failed repeatedly     ";
    XERRWD(msg, 50, 205, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      or with ABS(H) = HMIN   ";
    XERRWD(msg, 30, 205, 0, 0, 0, 0, 2, dls1_.tn, dls1_.h);
    istate = -5;
    goto LABEL_560;
// rwork length too small to proceed. -----------------------------------
LABEL_550:
    msg = "DLSODAR- At current T(=R1), RWORK length too small";
    XERRWD(msg, 50, 206, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      to proceed.  The integration was otherwise successful.";
    XERRWD(msg, 60, 206, 0, 0, 0, 0, 1, dls1_.tn, 0.0);
    istate = -7;
    goto LABEL_580;
// IWORK length too small to proceed. -----------------------------------
LABEL_555:
    msg = "DLSODAR- At current T(=R1), IWORK length too small";
    XERRWD(msg, 50, 207, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      to proceed.  The integration was otherwise successful.";
    XERRWD(msg, 60, 207, 0, 0, 0, 0, 1, dls1_.tn, 0.0);
    istate = -7;
    goto LABEL_580;
// Compute IMXER if relevant. -------------------------------------------
LABEL_560:
    big   = 0.0;
    imxer = 1;
    for (i = 1; i <= dls1_.n; ++i) {
        size = std::abs(ARRAYF(rwork, i + dls1_.lacor - 1) * ARRAYF(rwork, i + dls1_.lewt - 1));
        if (big >= size) continue;
        big   = size;
        imxer = i;
    }
    ARRAYF(iwork, 16) = imxer;
// Set Y vector, T, and optional outputs. -------------------------------
LABEL_580:
    for (i = 1; i <= dls1_.n; ++i) {
        ARRAYF(y, i) = ARRAYF(rwork, i + dls1_.lyh - 1);
    }
    t = dls1_.tn;
    ARRAYF(rwork, 11) = dls1_.hu;
    ARRAYF(rwork, 12) = dls1_.h;
    ARRAYF(rwork, 13) = dls1_.tn;
    ARRAYF(rwork, 15) = dlsa_.tsw;
    ARRAYF(iwork, 11) = dls1_.nst;
    ARRAYF(iwork, 12) = dls1_.nfe;
    ARRAYF(iwork, 13) = dls1_.nje;
    ARRAYF(iwork, 14) = dls1_.nqu;
    ARRAYF(iwork, 15) = dls1_.nq;
    ARRAYF(iwork, 19) = dlsa_.mused;
    ARRAYF(iwork, 20) = dls1_.meth;
    ARRAYF(iwork, 10) = dlsr_.nge;
    dlsr_.tlast       = t;
    return;
// -----------------------------------------------------------------------
//  Block I.
//  The following block handles all error returns due to illegal input
//  (ISTATE = -3), as detected before calling the core integrator.
//  First the error message routine is called.  If the illegal input
//  is a negative ISTATE, the run is aborted (apparent infinite loop).
// -----------------------------------------------------------------------
LABEL_601:
    msg = "DLSODAR-  ISTATE(=I1) illegal.";
    XERRWD(msg, 30, 1, 0, 1, istate, 0, 0, 0.0, 0.0);
    if (istate < 0) goto LABEL_800;
    goto LABEL_700;
LABEL_602:
    msg = "DLSODAR-  ITASK (=I1) illegal.";
    XERRWD(msg, 30, 2, 0, 1, itask, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_603:
    msg = "DLSODAR-  ISTATE>1 but DLSODAR not initialized.";
    XERRWD (msg, 50, 3, 0, 0, 0, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_604:
    msg = "DLSODAR-  NEQ (=I1) < 1    ";
    XERRWD (msg, 30, 4, 0, 1, neq, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_605:
    msg = "DLSODAR-  ISTATE = 3 and NEQ increased (I1 to I2).";
    XERRWD (msg, 50, 5, 0, 2, dls1_.n, neq, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_606:
    msg = "DLSODAR-  ITOL (=I1) illegal. ";
    XERRWD (msg, 30, 6, 0, 1, itol, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_607:
    msg = "DLSODAR-  IOPT (=I1) illegal. ";
    XERRWD (msg, 30, 7, 0, 1, iopt, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_608:
    msg = "DLSODAR-  JT (=I1) illegal.   ";
    XERRWD (msg, 30, 8, 0, 1, jt, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_609:
    msg = "DLSODAR-  ML (=I1) illegal: <0 or >=NEQ (=I2)";
    XERRWD (msg, 50, 9, 0, 2, ml, neq, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_610:
    msg = "DLSODAR-  MU (=I1) illegal: <0 or >=NEQ (=I2)";
    XERRWD (msg, 50, 10, 0, 2, mu, neq, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_611:
    msg = "DLSODAR-  IXPR (=I1) illegal. ";
    XERRWD (msg, 30, 11, 0, 1, dlsa_.ixpr, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_612:
    msg = "DLSODAR-  MXSTEP (=I1) .lt. 0 ";
    XERRWD (msg, 30, 12, 0, 1, dls1_.mxstep, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_613:
    msg = "DLSODAR-  MXHNIL (=I1) .lt. 0 ";
    XERRWD(msg, 30, 13, 0, 1, dls1_.mxhnil, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_614:
    msg = "DLSODAR-  TOUT (=R1) behind T (=R2)     ";
    XERRWD(msg, 40, 14, 0, 0, 0, 0, 2, tout, t);
    msg = "      Integration direction is given by H0 (=R1)  ";
    XERRWD(msg, 50, 14, 0, 0, 0, 0, 1, h0, 0.0);
    goto LABEL_700;
LABEL_615:
    msg = "DLSODAR-  HMAX (=R1) < 0.0 ";
    XERRWD (msg, 30, 15, 0, 0, 0, 0, 1, hmax, 0.0);
    goto LABEL_700;
LABEL_616:
    msg = "DLSODAR-  HMIN (=R1) < 0.0 ";
    XERRWD(msg, 30, 16, 0, 0, 0, 0, 1, dls1_.hmin, 0.0);
    goto LABEL_700;
LABEL_617:
    msg = "DLSODAR-  RWORK length needed, LENRW(=I1), exceeds LRW(=I2) ";
    XERRWD(msg, 60, 17, 0, 2, lenrw, lrw, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_618:
    msg = "DLSODAR-  IWORK length needed, LENIW(=I1), exceeds LIW(=I2) ";
    XERRWD(msg, 60, 18, 0, 2, leniw, liw, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_619:
    msg = "DLSODAR-  RTOL(I1) is R1 < 0.0       ";
    XERRWD(msg, 40, 19, 0, 1, i, 0, 1, rtoli, 0.0);
    goto LABEL_700;
LABEL_620:
    msg = "DLSODAR-  ATOL(I1) is R1 < 0.0       ";
    XERRWD(msg, 40, 20, 0, 1, i, 0, 1, atoli, 0.0);
    goto LABEL_700;
LABEL_621:
    ewti = ARRAYF(rwork, dls1_.lewt + i - 1);
    msg = "DLSODAR-  EWT(I1) is R1 <= 0.0        ";
    XERRWD(msg, 40, 21, 0, 1, i, 0, 1, ewti, 0.0);
    goto LABEL_700;
LABEL_622:
    msg = "DLSODAR- TOUT(=R1) too close to T(=R2) to start integration.";
    XERRWD(msg, 60, 22, 0, 0, 0, 0, 2, tout, t);
    goto LABEL_700;
LABEL_623:
    msg = "DLSODAR-  ITASK = I1 and TOUT (=R1) behind TCUR - dls1_.hu (= R2) ";
    XERRWD(msg, 60, 23, 0, 1, itask, 0, 2, tout, tp);
    goto LABEL_700;
LABEL_624:
    msg = "DLSODAR-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)  ";
    XERRWD(msg, 60, 24, 0, 0, 0, 0, 2, tcrit, dls1_.tn);
    goto LABEL_700;
LABEL_625:
    msg = "DLSODAR-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)  ";
    XERRWD(msg, 60, 25, 0, 0, 0, 0, 2, tcrit, tout);
    goto LABEL_700;
LABEL_626:
    msg = "DLSODAR-  At start of problem, too much accuracy  ";
    XERRWD(msg, 50, 26, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      requested for precision of machine..  See TOLSF (=R1) ";
    XERRWD(msg, 60, 26, 0, 0, 0, 0, 1, tolsf, 0.0);
    ARRAYF(rwork, 14) = tolsf;
    goto LABEL_700;
LABEL_627:
    msg = "DLSODAR-  Trouble in DINTDY. ITASK = I1, TOUT = R1";
    XERRWD(msg, 50, 27, 0, 1, itask, 0, 1, tout, 0.0);
    goto LABEL_700;
LABEL_628:
    msg = "DLSODAR-  MXORDN (=I1) .lt. 0 ";
    XERRWD(msg, 30, 28, 0, 1, dlsa_.mxordn, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_629:
    msg = "DLSODAR-  MXORDS (=I1) .lt. 0 ";
    XERRWD (msg, 30, 29, 0, 1, dlsa_.mxords, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_630:
     msg = "DLSODAR-  NG (=I1) < 0     ";
    XERRWD(msg, 30, 30, 0, 1, ng, 0, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_631:
    msg = "DLSODAR-  NG changed (from I1 to I2) illegally,   ";
    XERRWD(msg, 50, 31, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      i.e. not immediately after a root was found.";
    XERRWD(msg, 50, 31, 0, 2, dlsr_.ngc, ng, 0, 0.0, 0.0);
    goto LABEL_700;
LABEL_632:
    msg = "DLSODAR-  One or more components of g has a root  ";
    XERRWD(msg, 50, 32, 0, 0, 0, 0, 0, 0.0, 0.0);
    msg = "      too near to the initial point.    ";
    XERRWD(msg, 40, 32, 0, 0, 0, 0, 0, 0.0, 0.0);
//
LABEL_700:
    istate = -3;
    return;
//
LABEL_800:
    msg = "DLSODAR-  Run aborted.. apparent infinite loop.   ";
    XERRWD(msg, 50, 303, 2, 0, 0, 0, 0, 0.0, 0.0);
    return;
}

} // end namespace odepack_cpp