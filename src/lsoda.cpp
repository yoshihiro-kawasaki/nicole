#include "lsoda.hpp"
#include "linear.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

Lsoda_cpp::Lsoda_cpp()
{
    // Initialize arrays.

    ixpr_ = 0;
    mxords_ = 12;
    h_ = 0.0;
    tn_ = 0.0;
    itol_ = 2;

    mord_.resize(2);
    mord_[0] = 12;
    mord_[1] = 5; 

    sm1_.resize(12);
    sm1_[0] = 0.5;
    sm1_[1] = 0.575;
    sm1_[2] = 0.55;
    sm1_[3] = 0.45;
    sm1_[4] = 0.35;
    sm1_[5] = 0.25;
    sm1_[6] = 0.20;
    sm1_[7] = 0.15;
    sm1_[8] = 0.10;
    sm1_[9] = 0.075;
    sm1_[10] = 0.05;
    sm1_[11] = 0.025;

    el_.resize(13, 0.0);
    cm1_.resize(12, 0.0);
    cm2_.resize(5, 0.0);

    elco_.resize(12, std::vector<double>(13, 0.0));
    tesco_.resize(12, std::vector<double>(3, 0.0));

    mxstp0_ = 500;
    mxhnil_ = 10;
}

Lsoda_cpp::~Lsoda_cpp()
{

}

void Lsoda_cpp::TerminateIllegalInput(int *istate)
{
    if (illin_ == 5) {
        std::cerr << "[lsoda] repeated occurrence of illegal input. run aborted.. "
             "apparent infinite loop." << std::endl;
    } else {
        illin_++;
        *istate = -3;
    }
    return;
}


void Lsoda_cpp::TerminateVariousError(std::vector<double> &y, double *t, std::vector<int> &iworks, std::vector<double> &rworks)
{
    for (size_t i = 0; i < n_; ++i) y[i] = yh_[0][i];
    *t = tn_;
    illin_ = 0;
    rworks[10] = hu_;
    rworks[11] = h_;
    rworks[12] = tn_;
    rworks[14] = tsw_;
    iworks[10] = nst_;
    iworks[11] = nfe_;
    iworks[12] = nje_;
    iworks[13] = nqu_;
    iworks[14] = nq_;
    iworks[18] = mused_;
    iworks[19] = meth_;
    return;
}


/*
   The following block handles all successful returns from lsoda.
   If itask != 1, y is loaded from yh_ and t is set accordingly.
   *Istate is set to 2, the illegal input counter is zeroed, and the
   optional outputs are loaded into the work arrays before returning.
*/

void Lsoda_cpp::SuccessReturn(std::vector<double> &y, double *t, int itask, int ihit, double tcrit, int *istate,
    std::vector<int> &iworks, std::vector<double> &rworks)
{
    *t = tn_;

    for (int i = 0; i < n_; ++i) y[i] = yh_[0][i];

    if (itask == 4 || itask == 5) {
        if (ihit) *t = tcrit;
    }

    *istate = 2;
    illin_ = 0;
    rworks[10] = hu_;
    rworks[11] = h_;
    rworks[12] = tn_;
    rworks[14] = tsw_;
    iworks[10] = nst_;
    iworks[11] = nfe_;
    iworks[12] = nje_;
    iworks[13] = nqu_;
    iworks[14] = nq_;
    iworks[18] = mused_;
    iworks[19] = meth_;
}

void Lsoda_cpp::SetTolerance(std::vector<double> atol, std::vector<double> rtol)
{
    int n_atol = atol.size(), n_rtol = rtol.size();
    int i;

    atol_.resize(n_atol);
    rtol_.resize(n_rtol);

    for (i = 0; i < n_atol; ++i) {
        atol_[i] = atol[i];
    }

    for (i = 0; i < n_rtol; ++i) {
        rtol_[i] = rtol[i];
    }

    return;
}

/*
c references..
c 1.  alan c. hindmarsh,  odepack, a systematized collection of ode
c     solvers, in scientific computing, r. s. stepleman et al. (eds.),
c     north-holland, amsterdam, 1983, pp. 55-64.
c 2.  linda r. petzold, automatic selection of methods for solving
c     stiff and nonstiff systems of ordinary differential equations,
c     siam j. sci. stat. comput. 4 (1983), pp. 136-148.
c-----------------------------------------------------------------------
*/

void Lsoda_cpp::Lsoda(LSODA_ODE_SYSTEM_TYPE f, const size_t neq, std::vector<double> &y,
               double *t, double tout, int itol, int itask, int *istate, int iopt, 
               std::vector<int> &iworks, std::vector<double> &rworks, 
               LSODA_ODE_JACOBIAN jac, int jt, void *_data)
{
    assert(tout > *t);

    int mxstp0 = 500, mxhnl0 = 10;

    int iflag = 0, lenyh = 0, ihit = 0;

    double atoli = 0, ayi = 0, big = 0, h0 = 0, hmax = 0, hmx = 0, rh = 0,
           rtoli = 0, tcrit = 0, tdist = 0, tnext = 0, tol = 0, tolsf = 0, tp = 0,
           size = 0, sum = 0, w0 = 0;

    int i;

    /*
       Block a.
       This code block is executed on every call.
       It tests *istate and itask for legality and branches appropriately.
       If *istate > 1 but the flag init shows that initialization has not
       yet been done, an error return occurs.
       If *istate = 1 and tout = t, return immediately.
    */

    if (*istate < 1 || *istate > 3) { 
        std::cerr << "[lsoda] illegal istate = " << *istate << std::endl;
        TerminateIllegalInput(istate);
        return;
    }

    if (itask < 1 || itask > 5) { 
        std::cerr << "[lsoda] illegal itask = " << itask << std::endl;
        TerminateIllegalInput(istate);
        return;
    }

    if (init_ == 0 && (*istate == 2 || *istate == 3)) {
        std::cerr << "[lsoda] istate > 1 but lsoda not initialized" << std::endl;
        TerminateIllegalInput(istate);
        return;
    }

    /*
       Block b.
       The next code block is executed for the initial call ( *istate = 1 ),
       or for a continuation call with parameter changes ( *istate = 3 ).
       It contains checking of all inputs and various initializations.

       First check legality of the non-optional inputs neq, itol, iopt,
       jt, ml, and mu.
    */

    if (*istate == 1 || *istate == 3) {

        ntrep_ = 0;

        if (neq <= 0) {
            std::cerr << "[lsoda] neq = " << neq << " is less than 1." << std::endl;
            TerminateIllegalInput(istate);
            return;
        }

        if (*istate == 3 && neq > n_) { 
            std::cerr << "[lsoda] istate = 3 and neq increased" << std::endl;
            TerminateIllegalInput(istate);
            return;
        }

        n_ = neq; // number of first order ode

        if (itol_ < 1 || itol_ > 4) { 
            std::cerr << "[lsoda] itol = " << itol_ << " illegal" << std::endl;
            TerminateIllegalInput(istate);
            return;
        }

        if (iopt < 0 || iopt > 1) { 
            std::cerr << "[lsoda] iopt = " << iopt << " illegal" << std::endl;
            TerminateIllegalInput(istate);
            return;
        }

        if (jt == 3 || jt < 1 || jt > 5) {
            std::cerr << "[lsoda] jt = " << jt << " illegal" << std::endl;
            TerminateIllegalInput(istate);
            return;
        }

        jtyp_ = jt; // jacobian type indicator

        if (jt > 2) {

            ml_ = iworks[0]; //  default = 0
            mu_ = iworks[1]; //  default = 0

            if (ml_ < 0 || ml_ >= n_) { 
                std::cerr << "[lsoda] ml = " << ml_ << " not between 1 and neq" << std::endl;
                TerminateIllegalInput(istate);
                return;
            }

            if (mu_ < 0 || mu_ >= n_) { 
                std::cerr << "[lsoda] mu = " << mu_ << " not between 1 and neq" <<std:: endl;
                TerminateIllegalInput(istate);
                return;
            }

        }

        /* Next process and check the optional inpus.   */
        /* 
        iopt = an integer flag to specify whether or not any optional
               inputs are being used on this call.  input only.
               the optional inputs are listed separately below.
               iopt = 0 means no optional inputs are being used.
                        default values will be used in all cases.
               iopt = 1 means one or more optional inputs are being used.
        */

        /* Default options. (iopt = 0)*/
        if (iopt == 0) {
        
            ixpr_   = 0;
            mxstep_ = mxstp0;
            mxhnil_ = mxhnl0;
            hmxi_   = 0.0;
            hmin_   = 0.0;
            if (*istate == 1) { 
                h0 = 0.0;
                mxordn_ = mord_[0]; // = 12
                mxords_ = mord_[1]; // = 5
            }

            /* end if ( iopt == 0 )   */ /* Optional inputs.   */

        } else  {  /* if ( iopt = 1 )  */
            
            ixpr_= iworks[4]; 
            if (ixpr_ < 0 || ixpr_ > 1) { 
                std::cerr << "[lsoda] ixpr = " << ixpr_ << " is illegal" << std::endl;
                TerminateIllegalInput(istate);
                return;
            }

            mxstep_ = iworks[5]; 

            if (mxstep_ < 0) { 
                std::cerr << "[lsoda] mxstep = " << mxstep_ << " is illegal (mxstep < 0)" << std::endl;
                TerminateIllegalInput(istate);
                return;
            }

            if (mxstep_ == 0) mxstep_ = mxstp0;

            mxhnil_ = iworks[6]; 

            if (mxhnil_ < 0) { 
                std::cerr << "[lsoda] mxhnil = " << mxhnil_ << " is illegal (mxhnil < 0)" << std::endl;
                TerminateIllegalInput(istate);
                return;
            }

            if (mxhnil_ == 0) mxhnil_ = mxhnl0;

            if (*istate == 1) { 

                h0 = rworks[1];  

                mxordn_ = iworks[7]; 

                if (mxordn_ < 0) { 
                    std::cerr << "[lsoda] mxordn = " << mxordn_ << " is illegal (mxordn < 0)" << std::endl;
                    TerminateIllegalInput(istate);
                    return;
                }

                if (mxordn_ == 0) mxordn_ = 100;

                mxordn_ = std::min(mxordn_, mord_[0]);

                mxords_ = iworks[8];

                if (mxords_ < 0) { 
                    std::cerr << "[lsoda] mxords = " << mxords_ << " is illegal " << std::endl;
                    TerminateIllegalInput(istate);
                    return;
                }

                if (mxords_ == 0) mxords_ = 100; // if mxords is not given, use 100.

                mxords_ = std::min(mxords_, mord_[1]);

                if ((tout - *t) * h0 < 0.0) { 
                    std::cerr << "[lsoda] tout = " << tout << " behind t = " << *t
                         << ". integration direction is given by " << h0 << std::endl;
                    TerminateIllegalInput(istate);
                    return;
                }

            } /* end if ( *istate == 1 )  */

            hmax = rworks[2]; 

            if (hmax < 0.0) { 
                std::cerr << "[lsoda] hmax < 0.0" << std::endl;
                TerminateIllegalInput(istate);
                return;
            }

            hmxi_ = 0.0;

            if (hmax > 0.0) hmxi_ = 1.0 / hmax;

            hmin_ = rworks[3]; 

            if (hmin_ < 0.0) { 
                std::cerr << "[lsoda] hmin < 0.0" << std::endl;
                TerminateIllegalInput(istate);
                return;
            }

        } /* end else   */ /* end iopt = 1   */

    }  /* end if ( *istate == 1 || *istate == 3 )   */

    /*
       If *istate = 1, meth_ is initialized to 1.
       Also allocate memory for yh_, wm_, ewt, savf, acor, ipvt.
    */
   // 60
    if (*istate == 1) {

        /*
           If memory were not freed, *istate = 3 need not reallocate memory.
           Hence this section is not executed by *istate = 3.
        */

        sqrt_eta_ = sqrt(ETA);
        meth_     = 1;
        nyh_      = n_;
        lenyh     = 1 + std::max(mxordn_, mxords_);

        yh_.resize(lenyh, std::vector<double>(nyh_, 0.0));
        wm_.resize(nyh_, std::vector<double>(nyh_, 0.0));
        ewt_.resize(nyh_, 0.0);
        savf_.resize(nyh_, 0.0);
        acor_.resize(nyh_, 0.0);
        ipvt_.resize(nyh_, 0);

    }

    /*
       Check rtol and atol for legality.
    */
    if (*istate == 1 || *istate == 3) {

        rtoli = rtol_[0];
        atoli = atol_[0];

        for (i = 0; i < n_; i++) {
        
            if (itol_ >= 3) rtoli = rtol_[i];
            if (itol_ == 2 || itol_ == 4) atoli = atol_[i];

            if (rtoli < 0.0) {
                std::cerr << "[lsoda] rtol = " << rtoli << " is less than 0.0" << std::endl;
                TerminateIllegalInput(istate);
                return;
            }

            if (atoli < 0.0) { 
                std::cerr << "[lsoda] atol = " << atoli << " is less than 0.0" << std::endl;
                TerminateIllegalInput(istate);
                return;
            }

        } /* end for   */

    }   /* end if ( *istate == 1 || *istate == 3 )   */

    /* If *istate = 3, set flag to signal parameter changes to stoda. */
    if (*istate == 3)  jstart_ = -1;

    /*
       Block c.
       The next block is for the initial call only ( *istate = 1 ).
       It contains all remaining initializations, the initial call to f,
       and the calculation of the initial step size.
       The error weights in ewt are inverted after being loaded.
    */

    if (*istate == 1) {

        tn_     = *t;
        tsw_    = *t;
        maxord_ = mxordn_;

        if (itask == 4 || itask == 5) { 

            tcrit = rworks[0];

            if ((tcrit - tout) * (tout - *t) < 0.0) {
                std::cerr << "[lsoda] itask = 4 or 5 and tcrit behind tout" << std::endl;
                TerminateIllegalInput(istate);
                return;
            }

            if (h0 != 0.0 && (*t + h0 - tcrit) * h0 > 0.0) h0 = tcrit - *t;

        }

        jstart_ = 0;
        nhnil_  = 0;
        nst_    = 0;
        nje_    = 0;
        nslast_ = 0;
        hu_     = 0.;
        nqu_    = 0;
        mused_  = 0;
        miter_  = 0;
        ccmax_  = 0.30;
        maxcor_ = 3;
        msbp_   = 20;
        mxncf_  = 10;

        /* Initial call to f.  */
        assert(yh_.size() == lenyh);
        assert(yh_[0].size() == nyh_);

        (*f)(*t, y, yh_[1], _data);

        nfe_ = 1;

        /* Load the initial value vector in yh_.  */
        for (i = 0; i < n_; i++) yh_[0][i] = y[i];

        /* Load and invert the ewt array.  ( h_ is temporarily set to 1. ) */
        nq_ = 1;
        h_  = 1.0;

        SetErrorWeightVector(y);

        for (i = 0; i < n_; ++i) {

            if (ewt_[i] <= 0.0) {
                std::cerr << "[lsoda] ewt[" << i << "] = " << ewt_[i] << " <= 0.0" << std::endl;
                TerminateVariousError(y, t, iworks, rworks);
                return;
            }

            ewt_[i] = 1.0 / ewt_[i];

        }

        /*
           The coding below computes the step size, h0, to be attempted on the
           first step, unless the user has supplied a value for this.
           First check that tout - *t differs significantly from zero.
           A scalar tolerance quantity tol is computed, as max(rtol[i])
           if this is positive, or max(atol[i]/fabs(y[i])) otherwise, adjusted
           so as to be between 100*ETA and 0.001.
           Then the computed value h0 is given by

              h0^(-2) = 1. / ( tol * w0^2 ) + tol * ( norm(f) )^2

           where   w0     = max( fabs(*t), fabs(tout) ),
                   f      = the initial value of the vector f(t,y), and
                   norm() = the weighted vector norm used throughout, given by
                            the vmnorm function routine, and weighted by the
                            tolerances initially loaded into the ewt array.

           The sign of h0 is inferred from the initial values of tout and *t.
           fabs(h0) is made < fabs(tout-*t) in any case.
        */

        if (h0 == 0.0) { 

            tdist = std::abs(tout - *t);
            w0    = std::max(std::abs(*t), std::abs(tout));

            if (tdist < 2.0 * ETA * w0) { 
                std::cerr << "[lsoda] tout too close to t to start integration" << std::endl;
                TerminateIllegalInput(istate);
                return;
            }

            tol = rtol_[0];

            if (itol_ > 2) { 

                for (i = 1; i < n_; i++) tol = std::max(tol, rtol_[i]);

            }

            if (tol <= 0.0) { 

                atoli = atol_[0];

                for (i = 0; i < n_; ++i) {
                    if (itol_ == 2 || itol_ == 4) atoli = atol_[i];
                    ayi = std::abs(y[i]);
                    if (ayi != 0.0) tol = std::max(tol, atoli / ayi);
                }

            }

            tol = std::max(tol, 100.0 * ETA);
            tol = std::min(tol, 0.001);
            sum = NormOfVectorMax(n_, yh_[1], ewt_);
            sum = 1.0 / (tol * w0 * w0) + tol * sum * sum;
            h0  = 1.0 / std::sqrt(sum);
            h0  = std::min(h0, tdist);
            h0  = h0 * ((tout - *t >= 0.0) ? 1.0 : -1.0);

        } /* end if ( h0 == 0. )   */

        /*
           Adjust h0 if necessary to meet hmax bound.
        */
        rh = std::abs(h0) * hmxi_;
        if (rh > 1.0) h0 /= rh;

        /*
           Load h_ with h0 and scale yh_[2] by h0.
        */
        h_ = h0;
        for (i = 0; i < n_; ++i) yh_[1][i] *= h0;

    } /* if ( *istate == 1 )   */

    /*
       Block d.
       The next code block is for continuation calls only ( *istate = 2 or 3 )
       and is to check stop conditions before taking a step.
    */
    if (*istate == 2 || *istate == 3) {

        nslast_ = nst_;

        switch (itask) {

            case 1: 

                if ((tn_ - tout) * h_ >= 0.0) { 

                    InterpolatedValueOfDerivativeVector(tout, 0, y, &iflag);

                    if (iflag != 0) { 
                        std::cerr << "[lsoda] trouble from intdy, itask = " << itask << ", tout = " << tout << std::endl;
                        TerminateIllegalInput(istate);
                        return;
                    }

                    *t = tout;
                    *istate = 2;
                    illin_ = 0;
                    rworks[10] = hu_;
                    rworks[11] = h_;
                    rworks[12] = tn_;
                    rworks[14] = tsw_;
                    iworks[10] = nst_;
                    iworks[11] = nfe_;
                    iworks[12] = nje_;
                    iworks[13] = nqu_;
                    iworks[14] = nq_;
                    iworks[18] = mused_;
                    iworks[19] = meth_;
                    return;

                }

                break;

        case 2:

            break; 

        case 3:

            tp = tn_ - hu_ * (1.0 + 100.0 * ETA);

            if ((tp - tout) * h_ > 0.0) { 
                std::cerr << "[lsoda] itask = " << itask << " and tout behind tcur - hu" << std::endl;
                TerminateIllegalInput(istate);
                return;
            }

            if ((tn_ - tout) * h_ < 0.0) break;

            SuccessReturn(y, t, itask, ihit, tcrit, istate, iworks, rworks);

            return;

        case 4:

            tcrit = rworks[0]; 

            if ((tn_ - tcrit) * h_ > 0.0) { 
                std::cerr << "[lsoda] itask = 4 or 5 and tcrit behind tcur" << std::endl;
                TerminateIllegalInput(istate);
                return;
            }

            if ((tcrit - tout) * h_ < 0.0) { 
                std::cerr << "[lsoda] itask = 4 or 5 and tcrit behind tout" << std::endl;
                TerminateIllegalInput(istate);
                return;
            }

            if ((tn_ - tout) * h_ >= 0.0) { 

                InterpolatedValueOfDerivativeVector(tout, 0, y, &iflag);

                if (iflag != 0) { 
                    std::cerr << "[lsoda] trouble from intdy, itask = " << itask << ", tout = " << tout << std::endl;
                    TerminateIllegalInput(istate);
                    return;
                }

                *t = tout;

                *istate = 2;
                illin_ = 0;
                rworks[10] = hu_;
                rworks[11] = h_;
                rworks[12] = tn_;
                rworks[14] = tsw_;
                iworks[10] = nst_;
                iworks[11] = nfe_;
                iworks[12] = nje_;
                iworks[13] = nqu_;
                iworks[14] = nq_;
                iworks[18] = mused_;
                iworks[19] = meth_;
                return;

            } else {

                hmx  = std::abs(tn_) + std::abs(h_);
                ihit = std::abs(tn_ - tcrit) <= (100. * ETA * hmx);

                if (ihit) {
                    *t = tcrit;
                    SuccessReturn(y, t, itask, ihit, tcrit, istate, iworks, rworks);
                    return;
                }

                tnext = tn_ + h_ * (1.0 + 4.0 * ETA);
                if ((tnext - tcrit) * h_ <= 0.0) break; 
                h_ = (tcrit - tn_) * (1.0 - 4.0 * ETA);
                if (*istate == 2) jstart_ = -2;
                break;
            }

        case 5:

            tcrit = rworks[0]; 

            if ((tn_ - tcrit) * h_ > 0.0) { 
                std::cerr << "[lsoda] itask = 4 or 5 and tcrit behind tcur" << std::endl;
                TerminateIllegalInput(istate);
                return;
            }

            hmx  = std::abs(tn_) + std::abs(h_);
            ihit = std::abs(tn_ - tcrit) <= (100.0 * ETA * hmx);

            if (ihit) {
                *t = tcrit;
                SuccessReturn(y, t, itask, ihit, tcrit, istate, iworks, rworks); 
                return;
            }

            tnext = tn_ + h_ * (1.0 + 4.0 * ETA);
            if ((tnext - tcrit) * h_ <= 0.0) break; 
            h_ = (tcrit - tn_) * (1.0 - 4.0 * ETA);
            if (*istate == 2) jstart_ = -2;
            break;

        } /* end switch   */


    }   /* end if ( *istate == 2 || *istate == 3 )   */

    /*
       Block e.
       The next block is normally executed for all calls and contains
       the call to the one-step core integrator stoda.
       This is a looping point for the integration steps.
       First check for too many steps being taken, update ewt ( if not at
       start of problem).  Check for too much accuracy being requested, and
       check for h_ below the roundoff level in *t.
    */

    while (1) { 

        if (*istate != 1 || nst_ != 0) { 

            if ((nst_ - nslast_) >= mxstep_) {
                // std::cerr << "[lsoda] " << mxstep_ << " steps taken before reaching tout" << std::endl;
                *istate = -1;
                TerminateVariousError(y, t, iworks, rworks);
                return;
            }

            SetErrorWeightVector(yh_[0]);

            for (i = 0; i < n_; ++i) {

                if (ewt_[i] <= 0.0) {
                    std::cerr << "[lsoda] ewt[" << i << "] = " << ewt_[i] << " <= 0." << std::endl;
                    *istate = -6;
                    TerminateVariousError(y, t, iworks, rworks);
                    return;
                }

                ewt_[i] = 1.0 / ewt_[i];

            }

        } // end if (*istate != 1 || nst_ != 0)

        tolsf = ETA * NormOfVectorMax(n_, yh_[0], ewt_);

        if (tolsf > 0.01) { 

            tolsf = tolsf * 200.0;

            if (nst_ == 0) { 
                std::cerr << "lsoda -- at start of problem, too much accuracy" << std::endl;
                std::cerr << "         requested for precision of machine," << std::endl;
                std::cerr << "         suggested scaling factor = " << tolsf << std::endl;
                TerminateIllegalInput(istate);
                return;
            }

            std::cerr << "lsoda -- at t = " << *t << ", too much accuracy requested" << std::endl;
            std::cerr << "         for precision of machine, suggested" << std::endl;
            std::cerr << "         scaling factor = " << tolsf << std::endl;
            *istate = -2;
            TerminateVariousError(y, t, iworks, rworks);
            return;
        }

        if ((tn_ + h_) == tn_) { 

            nhnil_++;

            if (nhnil_ <= mxhnil_) { 
                // std::cerr << "lsoda -- warning..internal t = " << tn_ << " and h_ = " << h_ << " are" << std::endl;
                // std::cerr << "         such that in the machine, t + h_ = t on the next step" << std::endl;
                // std::cerr << "         solver will continue anyway." << std::endl;
                if (nhnil_ == mxhnil_) {
                    // std::cerr << "lsoda -- above warning has been issued " << nhnil_ << " times, " << std::endl;
                    // std::cerr << "       it will not be issued again for this problem" << std::endl;
                }
            }

        }

        /* Call stoda */
        Stoda(neq, y, f, jac, _data);
        // after call stoda, kflag_ = 0, -1, -2

        if (kflag_ == 0) {

            /*
               Block f.
               The following block handles the case of a successful return from the
               core integrator ( kflag = 0 ).
               If a method switch was just made, record tsw, reset maxord,
               set jstart to -1 to signal stoda to complete the switch,
               and do extra printing of data if ixpr = 1.
               Then, in any case, check for stop conditions.
            */

            init_ = 1;

            if (meth_ != mused_) { 

                //  mused : the method indicator for the last successful step..
                //  1 means adams (nonstiff), 2 means bdf (stiff).

                tsw_    = tn_;
                maxord_ = mxordn_;

                if (meth_ == 2) maxord_ = mxords_;

                jstart_ = -1;

                if (ixpr_) { 

                    if (meth_ == 2) {
                        std::cerr << "[lsoda] a switch to the stiff method has occurred " << std::endl;
                    }

                    if (meth_ == 1) {
                        std::cerr << "[lsoda] a switch to the nonstiff method has occurred" << std::endl;
                    }

                }

            } /* end if ( meth_ != mused )   */

            /*
               itask = 1.
               If tout has been reached, interpolate.
            */
            if (itask == 1) { 

                if ((tn_ - tout) * h_ < 0.0) continue; // while loop continue

                InterpolatedValueOfDerivativeVector(tout, 0, y, &iflag);
                *t = tout;
                *istate = 2;
                illin_ = 0;
                rworks[10] = hu_;
                rworks[11] = h_;
                rworks[12] = tn_;
                rworks[14] = tsw_;
                iworks[10] = nst_;
                iworks[11] = nfe_;
                iworks[12] = nje_;
                iworks[13] = nqu_;
                iworks[14] = nq_;
                iworks[18] = mused_;
                iworks[19] = meth_;

                return;

            }

            /*
               itask = 2.
            */
            if (itask == 2) { 

                SuccessReturn(y, t, itask, ihit, tcrit, istate, iworks, rworks); 
                return;

            }

            /*
               itask = 3.
               Jump to exit if tout was reached.
            */
            if (itask == 3) { 

                if ((tn_ - tout) * h_ >= 0.0) {
                    SuccessReturn(y, t, itask, ihit, tcrit, istate, iworks, rworks);
                    return;
                }

                continue; 
            }

            /*
               itask = 4.
               See if tout or tcrit was reached.  Adjust h_ if necessary.
            */
            if (itask == 4) { 

                if ((tn_ - tout) * h_ >= 0.0) {

                    InterpolatedValueOfDerivativeVector(tout, 0, y, &iflag);
                    *t = tout;
                    *istate = 2;
                    illin_ = 0;
                    rworks[10] = hu_;
                    rworks[11] = h_;
                    rworks[12] = tn_;
                    rworks[14] = tsw_;
                    iworks[10] = nst_;
                    iworks[11] = nfe_;
                    iworks[12] = nje_;
                    iworks[13] = nqu_;
                    iworks[14] = nq_;
                    iworks[18] = mused_;
                    iworks[19] = meth_;
                    return;

                } else { 

                    hmx  = std::abs(tn_) + std::abs(h_);
                    ihit = std::abs(tn_ - tcrit) <= (100.0 * ETA * hmx );

                    if (ihit) {
                        SuccessReturn(y, t, itask, ihit, tcrit, istate, iworks, rworks);
                        return;
                    }

                    tnext = tn_ + h_ * (1.0 + 4.0 * ETA);
                    if ((tnext - tcrit) * h_ <= 0.0) continue;
                    h_ = (tcrit - tn_) * (1.0 - 4.0 * ETA);
                    jstart_ = -2;
                    continue; // while loop 先頭へ

                }
            } /* end if ( itask == 4 )   */

            /*
               itask = 5.
               See if tcrit was reached and jump to exit.
            */
            if (itask == 5) { 

                hmx  = std::abs(tn_) + std::abs(h_);
                ihit = std::abs(tn_ - tcrit) <= (100.0 * ETA * hmx);
                SuccessReturn(y, t, itask, ihit, tcrit, istate, iworks, rworks);
                return;

            }

        } /* end if ( kflag == 0 )   */

        /*
           kflag = -1, error test failed repeatedly or with fabs(h_) = hmin.
           kflag = -2, convergence failed repeatedly or with fabs(h_) = hmin.
        */
        if (kflag_ == -1 || kflag_ == -2) { 

            std::cerr << "lsoda -- at t = " << tn_ << " and step size h_ = " << h_ << " the" << std::endl;

            if (kflag_ == -1) { 
                std::cerr << "         error test failed repeatedly or" << std::endl;
                std::cerr << "         with fabs(h_) = hmin" << std::endl;
                *istate = -4;
            }

            if (kflag_ == -2) {
                std::cerr << "         corrector convergence failed repeatedly or" << std::endl;
                std::cerr << "         with fabs(h_) = hmin" << std::endl;
                *istate = -5;
            }

            big    = 0.0;
            imxer_ = 1;
            for (i = 0; i < n_; ++i) {

                size = std::abs(acor_[i]) * ewt_[i];

                if (big < size) {

                    big = size;
                    imxer_ = i + 1;

                }

            }

            TerminateVariousError(y, t, iworks, rworks);
            return;

        } /* end if ( kflag == -1 || kflag == -2 )   */

    }   /* end while   */

} /* end Lsoda */


void Lsoda_cpp::Stoda(const size_t neq, std::vector<double> &y, LSODA_ODE_SYSTEM_TYPE f,
               LSODA_ODE_JACOBIAN jac, void *_data)
{
    assert(neq == y.size());

    size_t corflag = 0, orderflag = 0;
    size_t m = 0, ncf = 0;
    int i = 0, j = 0, i1 = 0;
    double del = 0.0, delp = 0.0, dsm = 0.0, dup = 0.0, exup = 0.0, r = 0.0,
           rh = 0.0, rh_sub = 0.0, rhup = 0.0, told = 0.0;
    double pdh = 0.0, pnorm = 0.0;

    /*
       stoda performs one step of the integration of an initial value
       problem for a system of ordinary differential equations.
       Note.. stoda is independent of the value of the iteration method
       indicator miter, when this is != 0, and hence is independent
       of the type of chord method used, or the Jacobian structure.
       Communication with stoda is done with the following variables:

       jstart = an integer used for input only, with the following
                values and meanings:

                   0  perform the first step,
                 > 0  take a new step continuing from the last,
                  -1  take the next step with a new value of h_,
                      n, meth_, miter, and/or matrix parameters.
                  -2  take the next step with a new value of h_,
                      but with other inputs unchanged.

       kflag = a completion code with the following meanings:

                 0  the step was successful,
                -1  the requested error could not be achieved,
                -2  corrector convergence could not be achieved,
                -3  fatal error in prja or solsy.

       miter = corrector iteration method:

                 0  functional iteration,
                >0  a chord method corresponding to jacobian type jt.

    */
    kflag_ = 0;
    told   = tn_;
    ncf    = 0;
    ierpj_ = 0;
    iersl_ = 0;
    jcur_  = 0;
    delp   = 0.0;

    /*
       On the first call, the order is set to 1, and other variables are
       initialized.  rmax is the maximum ratio by which h_ can be increased
       in a single step.  It is initially 1.e4 to compensate for the small
       initial h_, but then is normally equal to 10.  If a filure occurs
       (in corrector convergence or error test), rmax is set at 2 for
       the next increase.
       cfode is called to get the needed coefficients for both methods.
    */
    if (jstart_ == 0) {

        lmax_  = maxord_ ;
        nq_    = 1;
        l_     = 2;
        ialth_ = 2;
        rmax_  = 10000.0;
        rc_    = 0.0;
        el0_   = 1.0;
        crate_ = 0.7;
        hold_  = h_;
        nslp_  = 0;
        ipup_  = miter_;

        /*
           Initialize switching parameters.  meth_ = 1 is assumed initially.
        */
        icount_ = 20;
        irflag_ = 0;
        pdest_  = 0.0;
        pdlast_ = 0.0;
        ratio_  = 5.0;

        Cfode(2); // Cfoda(int meth) 

        for (i = 0; i < 5; ++i) {
            cm2_[i] = tesco_[i][1] * elco_[i][i + 1];
        }

        Cfode(1);

        for (i = 0; i < 12; ++i) {
            cm1_[i] = tesco_[i][1] * elco_[i][i + 1];
        }

        Resetcoeff();

    } /* end if ( jstart == 0 )   */

        /*
        The following block handles preliminaries needed when jstart = -1.
        ipup is set to miter to force a matrix update.
        If an order increase is about to be considered ( ialth = 1 ),
        ialth is reset to 2 to postpone consideration one more step.
        If the caller has changed meth_, cfode is called to reset
        the coefficients of the method.
        If h_ is to be changed, yh_ must be rescaled.
        If h_ or meth_ is being changed, ialth is reset to l = nq + 1
        to prevent further changes in h_ for that many steps.
        */
    if (jstart_ == -1) {

        ipup_ = miter_;
        lmax_ = maxord_;

        if (ialth_ == 1) ialth_ = 2;

        if (meth_ != mused_) { 

            Cfode(meth_);
            ialth_ = l_;
            Resetcoeff(); 

        }

        if (h_ != hold_) { 

            rh = h_ / hold_;
            h_ = hold_;
            Scaleh(&rh, &pdh);

        }

    } /* if ( jstart == -1 )   */

    if (jstart_ == -2) {

        if (h_ != hold_) {

            rh = h_ / hold_;
            h_ = hold_;
            Scaleh(&rh, &pdh); 

        }

    } /* if ( jstart == -2 )   */

    /*
       Prediction.
       This section computes the predicted values by effectively
       multiplying the yh_ array by the pascal triangle matrix.
       rc is the ratio of new to old values of the coefficient h_ * el[1].
       When rc differs from 1 by more than ccmax, ipup is set to miter
       to force pjac to be called, if a jacobian is involved.
       In any case, prja is called at least every msbp steps.
    */
    while (1) {

        while (1) {

            if (std::abs(rc_ - 1.0) > ccmax_) ipup_ = miter_; 
            if (nst_ >= nslp_ + msbp_) ipup_ = miter_;
            tn_ += h_;

            for (j = nq_ - 1; j >= 0; j--) { // if nq_= 0, j = 0, i1 = 1
                for (i1 = j; i1 < nq_; i1++) {
                    for (i = 0; i < n_; i++) {

                        yh_[i1][i] += yh_[i1 + 1][i];

                    }
                }
            }

            pnorm = NormOfVectorMax(n_, yh_[0], ewt_);

            Correction(neq, y, f, &corflag, pnorm, &del, &delp, &told, &ncf, &rh, &m, jac, _data); 

            if (corflag == 0) break; // 内側のwhile loopから抜ける

            if (corflag == 1) {

                rh_sub = hmin_ / std::abs(h_);
                rh = std::max(rh, rh_sub);
                Scaleh(&rh, &pdh);
                continue; // 内側のwhile loop先頭へ

            }

            if (corflag == 2) {
                
                kflag_  = -2;
                hold_   = h_;
                jstart_ = 1;
                return; // Stodaから出る
            }

        } /* end inner while ( corrector loop )   */

        /*
           The corrector has converged.  jcur is set to 0
           to signal that the Jacobian involved may need updating later.
           The local error test is done now.
        */ 
        jcur_ = 0;
        if (m == 0) dsm = del / tesco_[nq_ - 1][1];
        if (m > 0) dsm = NormOfVectorMax(n_, acor_, ewt_) / tesco_[nq_ - 1][1];

        if (dsm <= 1.0) {

            /*
               After a successful step, update the yh_ array.
               Decrease icount by 1, and if it is -1, consider switching methods.
               If a method switch is made, reset various parameters,
               rescale the yh_ array, and exit.  If there is no switch,
               consider changing h_ if ialth = 1.  Otherwise decrease ialth by 1.
               If ialth is then 1 and nq < maxord, then acor is saved for
               use in a possible order increase on the next step.
               If a change in h_ is considered, an increase or decrease in order
               by one is considered also.  A change in h_ is made only if it is by
               a factor of at least 1.1.  If not, ialth is set to 3 to prevent
               testing for that many steps.
            */
            kflag_ = 0;
            nst_++;
            hu_    = h_;
            nqu_   = nq_;
            mused_ = meth_;

            for (j = 0; j < l_; ++j) {

                r = el_[j];

                for (i = 0; i < n_; ++i) {

                    yh_[j][i] += r * acor_[i];

                }

            }

            icount_--;

            if (icount_ < 0) {

                MethodSwitch(dsm, pnorm, &pdh, &rh);

                if (meth_ != mused_) { // method swith true

                    rh_sub = hmin_ / std::abs(h_);
                    rh = std::max(rh, rh_sub);
                    Scaleh(&rh, &pdh);
                    rmax_ = 10.0;
                    EndStoda();
                    break;

                }

            }

            /*
               No method switch is being made.  Do the usual step/order selection.
            */
            ialth_--;
            if (ialth_ == 0) {

                rhup = 0.0;
                if (l_!= (lmax_ + 1)) {

                    for (i = 0; i < n_; ++i) savf_[i] = acor_[i] - yh_[lmax_][i];
                    dup = NormOfVectorMax(n_, savf_, ewt_) / tesco_[nq_ - 1][2];
                    exup = 1.0 / (double)(l_ + 1);
                    rhup = 1.0 / (1.4 * std::pow(dup, exup) + 0.0000014);

                }

                OrderSwitch(&rhup, dsm, &pdh, &rh, &orderflag); 

                /*
                   No change in h_ or nq.
                */
                if (orderflag == 0) {

                    EndStoda(); 
                    break; // 外側while loopから抜ける (Stodaから出る)

                }

                /*
                   h_ is changed, but not nq.
                */
                if (orderflag == 1) {

                    rh_sub = hmin_ / std::abs(h_);
                    rh = std::max(rh, rh_sub);
                    Scaleh(&rh, &pdh);
                    rmax_ = 10.0;
                    EndStoda();
                    break; // 外側while loopから抜ける (Stodaから出る)

                }

                /*
                   both nq and h_ are changed.
                */
                if (orderflag == 2) {

                    Resetcoeff(); 
                    rh_sub = hmin_ / std::abs(h_);
                    rh = std::max(rh, rh_sub);
                    Scaleh(&rh, &pdh);
                    rmax_ = 10.0;
                    EndStoda();
                    break; // 外側while loopから抜ける (Stodaから出る)

                }

            } /* end if ( ialth == 0 )   */


            if (ialth_ > 1 || l_ == (lmax_ + 1)) {
        
                EndStoda();
                break;

            }

            for (i = 0; i < n_; i++) yh_[lmax_][i] = acor_[i];
            EndStoda();
            break;

            /* end if ( dsm <= 1. )   */

        } else {
            /*
            The error test failed.  kflag keeps track of multiple failures.
            Restore tn_ and the yh_ array to their previous values, and prepare
            to try the step again.  Compute the optimum step size for this or
            one lower.  After 2 or more failures, h_ is forced to decrease
            by a factor of 0.2 or less.
            */
            kflag_--;
            tn_ = told;

            for (j = nq_ - 1; j >= 0; j--) {
                for (i1 = j; i1 < nq_; ++i1) {
                    for (i = 0; i < n_; ++i) {

                        yh_[i1][i] -= yh_[i1 + 1][i];

                    }
                }
            }

            rmax_ = 2.0;

            if (std::abs(h_) <= hmin_ * 1.00001) {

                kflag_ = -1;
                hold_   = h_;
                jstart_ = 1;
                break; // while loop から抜ける

            }

            if (kflag_ > -3) {

                rhup = 0.0;
                OrderSwitch(&rhup, dsm, &pdh, &rh, &orderflag); 

                if (orderflag == 1 || orderflag == 0) {

                    if (orderflag == 0) rh = std::min(rh, 0.2); 
                    rh_sub = hmin_ * std::abs(h_);
                    rh = std::max(rh, rh_sub);
                    Scaleh(&rh, &pdh);

                }

                if (orderflag == 2) {

                    Resetcoeff();
                    rh_sub = hmin_ * std::abs(h_);
                    rh = std::max(rh, rh_sub);
                    Scaleh(&rh, &pdh);

                }

                continue;

                /* if ( kflag > -3 )   */

            } else {

                /*
                Control reaches this section if 3 or more failures have occurred.
                If 10 failures have occurred, exit with kflag = -1.
                It is assumed that the derivatives that have accumulated in the
                yh_ array have errors of the wrong order.  Hence the first
                derivative is recomputed, and the order is set to 1.  Then
                h_ is reduced by a factor of 10, and the step is retried,
                until it succeeds or h_ reaches hmin.
                */
                
                if (kflag_ == -10) {

                    kflag_ = -1;
                    hold_   = h_;
                    jstart_ = 1;
                    break;

                } else {

                    rh = 0.1;
                    rh_sub = hmin_ / std::abs(h_);
                    rh = std::max(rh, rh_sub);
                    h_ *= rh;
                    for (i = 0; i < n_; i++) {
                        y[i] = yh_[0][i];
                    }
                    (*f)(tn_, y, savf_, _data);
                    nfe_++;
                    for (i = 0; i < n_; ++i) yh_[1][i] = h_ * savf_[i];
                    ipup_ = miter_;
                    ialth_ = 5;
                    if (nq_ == 1) continue; 
                    nq_ = 1;
                    l_  = 2;

                    Resetcoeff();

                    continue;

                }
            } /* end else -- kflag <= -3 */

        }   /* end error failure handling   */

    }     /* end outer while   */

} // end Stoda


void Lsoda_cpp::Prja(const size_t neq, std::vector<double> &y, LSODA_ODE_SYSTEM_TYPE f,
              LSODA_ODE_JACOBIAN jac, void *_data)
{
    size_t i = 0, ier = 0, j = 0;
    double fac = 0.0, hl0 = 0.0, r = 0.0, r0 = 0.0, yj = 0.0;
    double con = 0.0;
    /*
       prja is called by stoda to compute and process the matrix
       P = I - h_ * el[1] * J, where J is an approximation to the Jacobian.
       Here J is computed by finite differencing.
       J, scaled by -h_ * el[1], is stored in wm_.  Then the norm of J ( the
       matrix norm consistent with the weighted max-norm on vectors given
       by vmnorm ) is computed, and J is overwritten by P.  P is then
       subjected to LU decomposition in preparation for later solution
       of linear systems with p as coefficient matrix.  This is done
       by dgefa if miter = 2, and by dgbfa if miter = 5.
    */
    nje_++;
    ierpj_ = 0;
    jcur_  = 1;
    hl0    = h_ * el0_;

    double r_max1;

    /*
       If miter = 1, call jac and multiply by scalar.
    */
   if (miter_ == 1) { 

        if (jac == nullptr) {
            std::cerr << "not jac" << std::endl;
            std::exit(1);
        }

        for (i = 0; i < n_; ++i) {
            for (j = 0; j < n_; ++j) {
                wm_[i][j] = 0.0;
            }
        }

        (*jac)(tn_, y, wm_, _data);

        con = -hl0;
        for (i = 0; i < n_; ++i) {
            for (j = 0; j < n_; ++j) {
                wm_[i][j] *= con;
            }
        } 

        /*
           Compute norm of Jacobian.
        */
        pdnorm_ = Fnorm(n_, wm_, ewt_) / std::abs(hl0);

        /*
           Add identity matrix.
        */
        for (i = 0; i < n_; i++) wm_[i][i] += 1.0;
        /*
           Do LU decomposition on P.
        */
        FactorMatrixByGaussianElimination(wm_, n_, ipvt_, &ier);
        if (ier != 0) ierpj_ = 1;


        return; 

   } else if (miter_ == 2) { // 200

    /*
       If miter = 2, make n calls to f to approximate J.
    */

        fac = NormOfVectorMax(n_, savf_, ewt_);
        r0 = 1000.0 * std::abs(h_) * ETA * ((double)n_) * fac;
        if (r0 == 0.0) r0 = 1.0;

        for (j = 0; j < n_; ++j) {

            yj = y[j];
            r_max1 = sqrt_eta_ * std::abs(yj);
            r = std::max(r_max1, r0 / ewt_[j]);
            y[j] += r;
            fac = -hl0 / r;
            (*f)(tn_, y, acor_, _data);

            for (i = 0; i < n_; ++i) {

                wm_[i][j] = (acor_[i] - savf_[i]) * fac;

            }

            y[j] = yj;
        }

        nfe_ += n_;

        /*
           Compute norm of Jacobian.
        */
        pdnorm_ = Fnorm(n_, wm_, ewt_) / std::abs(hl0);

        /*
           Add identity matrix.
        */
        for (i = 0; i < n_; ++i) wm_[i][i] += 1.0;
        /*
           Do LU decomposition on P.
        */
        FactorMatrixByGaussianElimination(wm_, n_, ipvt_, &ier);
        if (ier != 0) ierpj_ = 1;


        return;

    } else {

        fprintf(stderr, "[prja] miter != 2\n");
        return;

    }

} // end Prja


double Lsoda_cpp::Fnorm(int n, const std::vector<std::vector<double>> &a, const std::vector<double> &w)
{
/*
   This subroutine computes the norm of a full n by n matrix,
   stored in the array a, that is consistent with the weighted max-norm
   on vectors, with weights stored in the array w.

      fnorm = max(i=1,...,n) ( w[i] * sum(j=1,...,n) fabs( a[i][j] ) / w[j] )
*/
    double an = 0, sum = 0;
    int i, j;

    for (i = 0; i < n; ++i) {

        sum = 0.0;

        for (j = 0; j < n; ++j) {
            sum += std::abs(a[i][j]) / w[j];
        }

        an = std::max(an, sum * w[i]);

    }

    return an;
}


void Lsoda_cpp::EndStoda()
{

    double r = 1.0 / tesco_[nqu_ - 1][1];
    for (size_t i = 0; i < n_; ++i) acor_[i] *= r;

    hold_ = h_;
    jstart_ = 1;

    return;
}

void Lsoda_cpp::MethodSwitch(double dsm, double pnorm, double *pdh, double *rh)
{
    int lm1, lm1p1, lm2, lm2p1, nqm1, nqm2;
    double rh1, rh2, rh1it, exm2, dm2, exm1, dm1, alpha, exsm;

    /*
       We are current using an Adams method.  Consider switching to bdf.
       If the current order is greater than 5, assume the problem is
       not stiff, and skip this section.
       If the Lipschitz constant and error estimate are not polluted
       by roundoff, perform the usual test.
       Otherwise, switch to the bdf methods if the last step was
       restricted to insure stability ( irflag = 1 ), and stay with Adams
       method if not.  When switching to bdf with polluted error estimates,
       in the absence of other information, double the step size.

       When the estimates are ok, we make the usual test by computing
       the step size we could have (ideally) used on this step,
       with the current (Adams) method, and also that for the bdf.
       If nq > mxords, we consider changing to order mxords on switching.
       Compare the two step sizes to decide whether to switch.
       The step size advantage must be at least ratio = 5 to switch.
    */
    if (meth_ == 1) { 

        if (nq_ > 5)  return; 

        if (dsm <= (100.0 * pnorm * ETA) || pdest_ == 0.0) {

            if (irflag_ == 0) return; 

            rh2  = 2.0;
            nqm2 = std::min(nq_, mxords_);

        } else {

            exsm  = 1.0 / (double)l_;
            rh1   = 1.0 / (1.2 * std::pow(dsm, exsm) + 0.0000012);
            rh1it = 2.0 * rh1;
            *pdh  = pdlast_ * std::abs(h_);
            if ((*pdh * rh1) > 0.00001) rh1it = sm1_[nq_ - 1] / *pdh;
            rh1   = std::min(rh1, rh1it);

            if (nq_ > mxords_) {

                nqm2  = mxords_;
                lm2   = mxords_ + 1;
                exm2  = 1.0 / (double)lm2;
                lm2p1 = lm2 + 1;
                dm2   = NormOfVectorMax(n_, yh_[lm2p1 - 1], ewt_) / cm2_[mxords_ - 1];
                rh2 = 1.0 / (1.2 * std::pow(dm2, exm2) + 0.0000012);

            } else {

                dm2  = dsm * (cm1_[nq_ - 1] / cm2_[nq_ - 1]);
                rh2  = 1.0 / (1.2 * std::pow(dm2, exsm) + 0.0000012);
                nqm2 = nq_;

            }

            if (rh2 < ratio_ * rh1) return; 
        }

        /*
           The method switch test passed.  Reset relevant quantities for bdf.
        */
        *rh     = rh2;
        icount_ = 20;
        meth_   = 2;
        miter_  = jtyp_;
        pdlast_ = 0.0;
        nq_     = nqm2;
        l_      = nq_ + 1;

        return;

    } /* end if ( meth_ == 1 )   */

    /*
       We are currently using a bdf method, considering switching to Adams.
       Compute the step size we could have (ideally) used on this step,
       with the current (bdf) method, and also that for the Adams.
       If nq > mxordn, we consider changing to order mxordn on switching.
       Compare the two step sizes to decide whether to switch.
       The step size advantage must be at least 5/ratio = 1 to switch.
       If the step size for Adams would be so small as to cause
       roundoff pollution, we stay with bdf.
    */
    // meth_ == 2
    exsm = 1.0 / (double)l_;
    if (mxordn_ < nq_) {

        nqm1 = mxordn_;
        lm1 = mxordn_ + 1;
        exm1 = 1.0 / (double)lm1;
        lm1p1 = lm1 + 1;
        dm1 = NormOfVectorMax(n_, yh_[lm1p1 - 1], ewt_) / cm1_[mxordn_ - 1];
        rh1 = 1.0 / (1.2 * std::pow(dm1, exm1) + 0.0000012);

    } else {

        dm1 = dsm * (cm2_[nq_ - 1] / cm1_[nq_ - 1]);
        rh1 = 1.0 / (1.2 * pow(dm1, exsm) + 0.0000012);
        nqm1 = nq_;
        exm1 = exsm;

    }

    rh1it = 2.0 * rh1;
    *pdh = pdnorm_ * std::abs(h_);
    if ((*pdh * rh1) > 0.00001) rh1it = sm1_[nqm1 - 1] / *pdh;
    rh1 = std::min(rh1, rh1it);
    rh2 = 1.0 / (1.2 * std::pow(dsm, exsm) + 0.0000012);
    if ((rh1 * ratio_) < (5.0 * rh2)) return; 
    alpha = std::max(0.001, rh1);
    dm1 *= std::pow(alpha, exm1);
    if (dm1 <= 1000.0 * ETA * pnorm) return; 

    /*
       The switch test passed.  Reset relevant quantities for Adams.
    */
    *rh     = rh1;
    icount_ = 20;
    meth_   = 1;
    miter_  = 0;
    pdlast_ = 0.;
    nq_     = nqm1;
    l_      = nq_ + 1;

    return;  

} // MethodSwitch

void Lsoda_cpp::Cfode(int meth)
{
    int i, nq, nqm1, nqp1;
    double agamq, fnq, fnqm1, pc[13], pint, ragq, rqfac, rq1fac, tsign, xpin;

    /*
       cfode is called by the integrator routine to set coefficients
       needed there.  The coefficients for the current method, as
       given by the value of meth_, are set for all orders and saved.
       The maximum order assumed here is 12 if meth_ = 1 and 5 if meth_ = 2.
       ( A smaller value of the maximum order is also allowed. )
       cfode is called once at the beginning of the problem, and
       is not called again unless and until meth_ is changed.

       The elco array contains the basic method coefficients.
       The coefficients el[i], 1 < i < nq+1, for the method of
       order nq are stored in elco[nq][i].  They are given by a generating
       polynomial, i.e.,

          l(x) = el[1] + el[2]*x + ... + el[nq+1]*x^nq.

       For the implicit Adams method, l(x) is given by

          dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),   l(-1) = 0.

       For the bdf methods, l(x) is given by

          l(x) = (x+1)*(x+2)*...*(x+nq)/k,

       where   k = factorial(nq)*(1+1/2+...+1/nq).

       The tesco array contains test constants used for the
       local error test and the selection of step size and/or order.
       At order nq, tesco[nq][k] is used for the selection of step
       size at order nq-1 if k = 1, at order nq if k = 2, and at order
       nq+1 if k = 3.
    */

    if (meth == 1) {

        elco_[0][0]   = 1.0;
        elco_[0][1]   = 1.0;
        tesco_[0][0]  = 0.0;
        tesco_[0][1]  = 2.0;
        tesco_[1][0]  = 1.0;
        tesco_[11][2] = 0.0;
        pc[0]         = 1.0;
        rqfac         = 1.0;

        for (nq = 1; nq < 12; ++nq) {

            /*
               The pc array will contain the coefficients of the polynomial
                  p(x) = (x+1)*(x+2)*...*(x+nq-1).
               Initially, p(x) = 1.
            */

            rq1fac = rqfac;
            rqfac  = rqfac / (double)(nq + 1);
            nqm1   = nq - 1;
            fnqm1  = (double)(nqm1 +1);
            nqp1   = nq + 1;

            /*
               Form coefficients of p(x)*(x+nq-1).
            */

            pc[nq] = 0.0;
            for (i = nq; i >= 1; i--) {
                pc[i] = pc[i - 1] + fnqm1 * pc[i];
            }
            pc[0] = fnqm1 * pc[0];

            /*
               Compute integral, -1 to 0, of p(x) and x*p(x).
            */
            pint  = pc[0];
            xpin  = pc[0] / 2.0;
            tsign = 1.0;

            for (i = 1; i <= nq; ++i) {

                tsign = -tsign;
                pint += tsign * pc[i] / (double)(i + 1);
                xpin += tsign * pc[i] / (double)(i + 2);

            }

            /*
               Store coefficients in elco and tesco.
            */
            elco_[nq][0] = pint * rq1fac;
            elco_[nq][1] = 1.0;

            for (i = 1; i <= nq; ++i) {
                elco_[nq][i + 1] = rq1fac * pc[i] / (double)(i + 1);
            }
            agamq = rqfac * xpin;
            ragq  = 1.0 / agamq;
            tesco_[nq][1] = ragq;

            if (nq < 11) tesco_[nqp1][0] = ragq * rqfac / (double)(nqp1 + 1);
            tesco_[nqm1][2] = ragq;

        } /* end for nq   */

        return;

    } /* end if ( meth_ == 1 )   */

    /* meth_ = 2. */
    pc[0]  = 1.0;
    rq1fac = 1.0;

    /*
       The pc array will contain the coefficients of the polynomial
          p(x) = (x+1)*(x+2)*...*(x+nq).
       Initially, p(x) = 1.
    */
    for (nq = 0; nq < 5; nq++) {

        fnq = (double)(nq + 1);
        nqp1 = nq + 1;

        /*
           Form coefficients of p(x)*(x+nq).
        */
        pc[nqp1] = 0.0;

        for (i = nq + 1; i >= 1; --i) {
            pc[i] = pc[i - 1] + fnq * pc[i];
        }
        pc[0] *= fnq;

        /*
           Store coefficients in elco and tesco.
        */
        for (i = 0; i <= nqp1; ++i) {
            elco_[nq][i] = pc[i] / pc[1];
        }
        elco_[nq][1]  = 1.0;
        tesco_[nq][0] = rq1fac;
        tesco_[nq][1] = ((double)(nqp1 + 1)) / elco_[nq][0];
        tesco_[nq][2] = ((double)(nq + 3)) / elco_[nq][0];
        rq1fac /= fnq;

    }

    return;

} // end Cfode


void Lsoda_cpp::OrderSwitch(double *rhup, double dsm, double *pdh, double *rh, size_t *orderflag)
{
    size_t newq = 0;
    double exsm, rhdn, rhsm, ddn, exdn, r, pdh_sub;

    *orderflag = 0;

    exsm = 1.0 / (double)l_;
    rhsm = 1.0 / (1.2 * std::pow(dsm, exsm) + 0.0000012);

    rhdn = 0.0;

    if (nq_ != 1) {

        ddn = NormOfVectorMax(n_, yh_[l_ - 1], ewt_) / tesco_[nq_ - 1][0];
        exdn = 1.0 / (double)nq_;
        rhdn = 1.0 / (1.3 * pow(ddn, exdn) + 0.0000013);

    }

    /*
       If meth_ = 1, limit rh accordinfg to the stability region also.
    */
    if (meth_ == 1) {

        pdh_sub = std::abs(h_) * pdlast_;
        *pdh = std::max(pdh_sub, 0.000001);
        if (l_ < lmax_) *rhup = std::min(*rhup, sm1_[l_ - 1] / *pdh);
        rhsm = std::min(rhsm, sm1_[nq_ - 1] / *pdh);
        if (nq_ > 1) rhdn = std::min(rhdn, sm1_[nq_ - 2] / *pdh);
        pdest_ = 0.0;

    }

    if (rhsm >= *rhup) {

        if (rhsm >= rhdn) {

            newq = nq_;
            *rh = rhsm;

        } else {

            newq = nq_ - 1;
            *rh = rhdn;
            if (kflag_ < 0 && *rh > 1.0) *rh = 1.0;

        }

    } else {

        if (*rhup <= rhdn) {

            newq = nq_ - 1;
            *rh = rhdn;
            if (kflag_ < 0 && *rh > 1.0) *rh = 1.0;

        } else {

            *rh = *rhup;
            if (*rh >= 1.1) {

                r = el_[l_ - 1] / (double)l_;
                nq_ = l_;
                l_ = nq_ + 1;
                for (int i = 0; i < n_; i++) yh_[l_ - 1][i] = acor_[i] * r;

                *orderflag = 2;

                return; 

            } else {

                ialth_ = 3;
                return; 

            }

        }

    }

    /*
       If meth_ = 1 and h_ is restricted by stability, bypass 10 percent test.
    */
    if (meth_ == 1) { 

        if ((*rh * *pdh * 1.00001) < sm1_[newq - 1]) {

            if (kflag_ == 0 && *rh < 1.1) {


                ialth_ = 3;
                return; 

            }

        }

    } else {

        if (kflag_ == 0 && *rh < 1.1) {

            ialth_ = 3;
            return; 

        }

    }

    if (kflag_ <= -2) *rh = std::min(*rh, 0.2);

    /*
       If there is a change of order, reset nq, l, and the coefficients.
       In any case h_ is reset according to rh and the yh_ array is rescaled.
       Then exit or redo the step.
    */
    if (newq == nq_) {

        *orderflag = 1;
        return; 

    }

    nq_ = newq;
    l_ = nq_ + 1;
    *orderflag = 2; 

    return;

} // OrderSwitch

void Lsoda_cpp::Scaleh(double *rh, double *pdh)
{
    double r, a;
    int i;

    /*
       If h_ is being changed, the h_ ratio rh is checked against rmax, hmin,
       and hmxi, and the yh_ array is rescaled.  ialth is set to l = nq + 1
       to prevent a change of h_ for that many steps, unless forced by a
       convergence or error test failure.
    */
    *rh = std::min(*rh, rmax_);
    a   = std::abs(h_) * hmxi_ * (*rh);
    *rh = *rh / std::max(1.0, a);

    /*
       If meth_ = 1, also restrict the new step size by the stability region.
       If this reduces h_, set irflag to 1 so that if there are roundoff
       problems later, we can assume that is the cause of the trouble.
    */
    if (meth_ == 1) { // if meth_ == 2 goto 178 

        irflag_ = 0;

        a    = std::abs(h_) * pdlast_;
        *pdh = std::max(a, 0.000001);

        if ((*rh * *pdh * 1.00001) >= sm1_[nq_ - 1]) { 

            *rh = sm1_[nq_ - 1] / *pdh;
            irflag_ = 1;

        }

    }

    r = 1.0;

    for (size_t j = 1; j < l_; ++j) {

        r *= *rh;
        for (i = 0; i < n_; i++) yh_[j][i] *= r;

    }

    h_    *= (*rh);
    rc_   *= (*rh);
    ialth_ = l_;

    return;

} // Scaleh


void Lsoda_cpp::Correction(const size_t neq, std::vector<double> &y, LSODA_ODE_SYSTEM_TYPE f,
                    size_t *corflag, double pnorm, double *del, double *delp,
                    double *told, size_t *ncf, double *rh, size_t *m,
                    LSODA_ODE_JACOBIAN jac, void *_data)
{
    double rm = 0.0, rate = 0.0, dcon = 0.0;
    int i;

    /*
        *corflag = 0 : corrector converged,
        1 : step size to be reduced, redo prediction,
        2 : corrector cannot converge, failure flag.
    */

    /*
       Up to maxcor corrector iterations are taken.  A convergence test is
       made on the r.m.s. norm of each correction, weighted by the error
       weight vector ewt.  The sum of the corrections is accumulated in the
       vector acor[i].  The yh_ array is not altered in the corrector loop.
    */

    *m       = 0;
    *corflag = 0;
    *del     = 0.0;

    for (i = 0; i < n_; ++i) {
        y[i] = yh_[0][i];
    }

    (*f)(tn_, y, savf_, _data);

    nfe_++;

    /*
       If indicated, the matrix P = I - h_ * el[1] * J is reevaluated and
       preprocessed before starting the corrector iteration.  ipup is set
       to 0 as an indicator that this has been done.
    */
    while (1) {

        if (*m == 0) {  

            if (ipup_ > 0) { 

                Prja(neq, y, f, jac, _data);

                ipup_  = 0;
                rc_    = 1.0;
                nslp_  = nst_;
                crate_ = 0.7;

                if (ierpj_ != 0) { 

                    Corfailure(told, rh, ncf, corflag);
                    // corflag = 1 or 2
                    return; // Correctionから出る

                }

            }

            for (i = 0; i < n_; i++) acor_[i] = 0.0;

        } /* end if ( *m == 0 )   */

        if (miter_ == 0) { 

            /*
               In case of functional iteration, update y directly from
               the result of the last function evaluation.
            */

            for (i = 0; i < n_; ++i) {

                savf_[i] = h_ * savf_[i] - yh_[1][i];
                y[i] = savf_[i] - acor_[i];

            }

            *del = NormOfVectorMax(n_, y, ewt_);

            for (i = 0; i < n_; ++i) {

                y[i] = yh_[0][i] + el_[0] * savf_[i]; 
                acor_[i] = savf_[i];

            }

            /* end functional iteration   */ 

        } else { 

            /*
           In the case of the chord method, compute the corrector error,
           and solve the linear system with that as right-hand side and
           P as coefficient matrix.
             */

            for (i = 0; i < n_; ++i) {

                y[i] = h_ * savf_[i] - (yh_[1][i] + acor_[i]);

            }

            Solsy(y);

            *del = NormOfVectorMax(n_, y, ewt_);

            for (i = 0; i < n_; ++i) {

                acor_[i] += y[i];
                y[i] = yh_[0][i] + el_[0] * acor_[i];

            }

        } /* end chord method   */

        /*
           Test for convergence.  If *m > 0, an estimate of the convergence
           rate constant is stored in crate, and this is used in the test.

           We first check for a change of iterates that is the size of
           roundoff error.  If this occurs, the iteration has converged, and a
           new rate estimate is not formed.
           In all other cases, force at least two iterations to estimate a
           local Lipschitz constant estimate for Adams method.
           On convergence, form pdest = local maximum Lipschitz constant
           estimate.  pdlast is the most recent nonzero estimate.
        */
        if (*del <= 100.0 * pnorm * ETA) break; // while loop 抜ける corflag = 0 Stodaへ戻る

        if (*m != 0 || meth_ != 1) {

            if (*m != 0) {

                rm = 1024.0;
                if (*del <= (1024.0 * *delp)) rm = *del / *delp;
                rate   = std::max(rate, rm);
                crate_ = std::max(0.2 * crate_, rm);

            }

            dcon = *del * std::min(1.0, 1.5 * crate_) / (tesco_[nq_ - 1][1] * conit_);

            if (dcon <= 1.0) {

                double rr = rate / std::abs(h_ * el_[0]);
                pdest_ = std::max(pdest_, rr);
                if (pdest_ != 0.) pdlast_ = pdest_;
                break; // while loop抜ける corflag = 0  Stodaへ戻る

            }

        }

        /*
           The corrector iteration failed to converge.
           If miter != 0 and the Jacobian is out of date, prja is called for
           the next try.   Otherwise the yh_ array is retracted to its values
           before prediction, and h_ is reduced, if possible.  If h_ cannot be
           reduced or mxncf failures have occured, exit with corflag = 2.
        */
        (*m)++;
        if (*m == maxcor_ || (*m >= 2 && *del > 2.0 * *delp)) {

            if (miter_ == 0 || jcur_ == 1) {

                Corfailure(told, rh, ncf, corflag); 
                // corflag = 1 or 2
                return; // Correctionから出る

            }

            ipup_ = miter_;

            /*
               Restart corrector if Jacobian is recomputed. 
            */

            *m   = 0;
            rate = 0.0;
            *del = 0.0;
            for (i = 0; i < n_; i++) y[i] = yh_[0][i];

            (*f)(tn_, y, savf_, _data);

            nfe_++;

            // while loop先頭へ戻る

        } else {
            /*
            Iterate corrector.
            */
            *delp = *del;
            (*f)(tn_, y, savf_, _data);
            nfe_++;

            // while loop先頭へ戻る
        }

    } /* end while   */

} // end Correction

void Lsoda_cpp::Solsy(std::vector<double> &y)
{
/*
   This routine manages the solution of the linear system arising from
   a chord iteration.  It is called if miter != 0.
   If miter is 2, it calls dgesl to accomplish this.
   If miter is 5, it calls dgbsl.

   y = the right-hand side vector on input, and the solution vector
       on output.
*/
    iersl_ = 0;

    if (miter_ == 1 || miter_ == 2) {

        SolveByGaussianElimination(wm_, n_, ipvt_, y, 0);

    } else {

        std::cerr << "solsy -- miter != 1 or 2" << std::endl;
        return;

    }

    return;  

} // End Solsy

void Lsoda_cpp::Corfailure(double *told, double *rh, size_t *ncf, size_t *corflag)
{
    int j, i1, i;

    ncf++;
    rmax_ = 2.0;
    tn_ = *told;

    for (j = nq_ - 1; j >= 0; --j) {
        for (i1 = j; i1 < nq_; ++i1) {
            for (i = 0; i < n_; ++i) {

                yh_[i1][i] -= yh_[i1 + 1][i];

            }
        }
    }

    if (std::abs(h_) <= hmin_ * 1.00001 || *ncf == mxncf_) { 
        
        *corflag = 2; // kflag = -2
        // hold_ = h_;
        // jstart_ = 1;
        return;

    }

    *corflag = 1;
    *rh      = 0.25;
    ipup_    = miter_;

    return;

} // end Corfailure


void Lsoda_cpp::Resetcoeff()
{
/*
   The el vector and related constants are reset
   whenever the order nq is changed, or at the start of the problem.
*/
    int i;

    for (i = 0; i < l_; ++i) el_[i] = elco_[nq_ - 1][i];
    rc_ = rc_ * el_[0] / el0_;
    el0_ = el_[0];
    conit_ = 0.5 / (double)(nq_ + 2);

    return;
}

void Lsoda_cpp::SetErrorWeightVector(const std::vector<double> &ycur)
{
    int i;
    switch (itol_) {

        case 1:

            for (i = 0; i < n_; ++i) {
                ewt_[i] = rtol_[0] * std::abs(ycur[i]) + atol_[0];
            }
            break;

        case 2:

            for (i = 0; i < n_; ++i) {
                ewt_[i] = rtol_[0] * std::abs(ycur[i]) + atol_[i];
            }
            break;

        case 3:

            for (i = 0; i < n_; ++i) {
                ewt_[i] = rtol_[i] * std::abs(ycur[i]) + atol_[0];
            }
            break;

        case 4:

            for (i = 0; i < n_; ++i) {
                ewt_[i] = rtol_[i] * std::abs(ycur[i]) + atol_[i];
            }
            break;

    }

    return;

} // end SetErrorWeightVector

double Lsoda_cpp::NormOfVectorMax(const size_t n, const std::vector<double> &v, const std::vector<double> &w)
{
    double vm = 0.0, a;
    int i;

    for (i = 0; i < n; ++i) {

        a = std::abs(v[i]) * w[i];
        vm = std::max(vm,  a);

    }

    return vm;
}


void Lsoda_cpp::InterpolatedValueOfDerivativeVector(double t, int k, std::vector<double> &dky, int *iflag)
{
    int ic, jp1 = 0;
    double c, r, s, tp;
    int i, j;
    size_t jj;

    *iflag = 0;

    if (k < 0 || k > (int)nq_) { 

        std::cerr << "[intdy] k = " << k << " illegal" << std::endl;
        *iflag = -1;
        return;

    }

    tp = tn_ - hu_ - 100.0 * ETA * (tn_ + hu_);

    if ((t - tp) * (t - tn_) > 0.0) { 

        std::cerr << "InterpolatedValueOfDerivativeVector -- t = " << t << " illegal. t not in interval tcur - hu to tcur" << std::endl;

        *iflag = -2;
        return;

    }

    s  = (t - tn_) / h_;
    ic = 1;

    for (jj = l_ - k; jj <= nq_; ++jj) ic *= jj;

    c = (double)ic;

    for (i = 0; i < n_; ++i) dky[i] = c * yh_[l_ - 1][i];

    for (j = nq_ - 1; j >= k; --j) {

        jp1 = j + 1;
        ic  = 1;

        for (jj = jp1 - k; jj <= j; ++jj) ic *= jj;

        c = (double)ic;

        for (i = 0; i < n_; ++i) {

            dky[i] = c * yh_[jp1 - 1][i] + s * dky[i];

        }

    }

    if (k == 0) return;

    r = pow(h_, (double)(-k));

    for (i = 0; i < n_; ++i) dky[i] *= r;

} // end InterpolatedValueOfDerivativeVector


void Lsoda_cpp::StaticNumber()
{
    std::cout << "nst = " << nst_ << " nfe = " << nfe_ << " nje = " << nje_ << std::endl;
}