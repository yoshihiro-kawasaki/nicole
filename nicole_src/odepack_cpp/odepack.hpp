/**
 * @note
 * 2025/01/21 kawasaki
*/

#ifndef ODEPACK_HPP_
#define ODEPACK_HPP_

#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>

/**
* @namespace odepack_cpp
* @brief
* @details
*/
namespace odepack_cpp
{

/**
* @def ARRAYF
* @brief C++ array index --> Fortran array index
* @details
*/
#ifndef ARRAYF
#define ARRAYF(a, i) (a[(i)-1])
#endif

/**
 * @def MATF
 * @brief Fortran type 2d-array index
 * @details
 *  Column major order
 *  a[nrows][ncols] = a[nrows * ncols]
 *  a[i][j] = a[(i-1) + (j-1) * (rows)]
 * 
*/
#ifndef MATF
#define MATF(a, rows, i, j) (a[(i-1) + (j-1) * (rows)])
#endif

/**
 * @fn SIGN
*/
template <typename T>
constexpr T SIGN(const T x) {
    return (x >= static_cast<T>(0) ? static_cast<T>(1) : static_cast<T>(-1));
}

/**
 * define
*/
using ODEPACK_FUNCTION   = void(*)(const int neq, const double t, const double *y, double *ydot, void *user_data);
using ODEPACK_JACOBIAN1  = void(*)(const int neq, const double t, const double *y, const int ml, const int mu, double *pd, const int nrowpd, void *user_data);
using ODEPACK_JACOBIAN2  = void(*)(const int neq, const double t, const double *y, const int j, int *ian, int *jan, double *pdj, void *user_data);
using ODEPACK_JACOBIAN3  = void(*)(ODEPACK_FUNCTION f, const int neq, const double t, double *y, double *ysv, double *rewt, double *fty, double *v, double hl0, double *wp, int *iwp, int ier, void *user_data);
using ODEPACK_JACOBIAN4  = void(*)(const int neq, const double t, const double *y, double *s, const int ml, const int mu, double *p, const int nrowp, void *user_data);
using ODEPACK_JACOBIAN5  = void(*)(const int neq, const double t, const double *y, double *s, const int mb, const int nb, double *pa, double *pb, double *pc, void *user_data);
using ODEPACK_JACOBIAN6  = void(*)(const int neq, const double t, const double *y, double *s, const int j, int *ian, int *jan, double *pdj, void *user_data);
using ODEPACK_CONSTRAINT = void(*)(const int neq, const double t, const double *y, const int ng, double *gout, void *user_data);
using ODEPACK_PSOL       = void(*)(const int neq, const double t, const double *y, double *fty, double *wk, double hl0, double *wp, int *iwp, double *b, int lr, int irt, void *user_data);
using ODEPACK_RESIDUAL   = void(*)(const int neq, const double t, const double *y, double *s, double *r, int ires, void *user_data);
using ODEPACK_ADDA1      = void(*)(const int neq, const double t, const double *y, const int ml, const int mu, double *p, const int nrowp, void *user_data);
using ODEPACK_ADDA2      = void(*)(const int neq, const double t, const double *y, const int mb, const int nb, double *pa, double *pb, double *pc, void *user_data);
using ODEPACK_ADDA3      = void(*)(const int neq, const double t, const double *y, const int j, int *ian, int *jan, double *p, void *user_data);

/**
 * @fn class Odepack
*/
class Odepack
{
public:
    Odepack() { };
    ~Odepack() { };

    // opkdmain.cpp
    void DLSODE(ODEPACK_FUNCTION f, const int neq, double *y, double &t, const double tout,
        const int itol, double *rtol, double *atol, const int itask, int &istate,
        const int iopt, double *rwork, const int lrw, int *iwork, const int liw,
        ODEPACK_JACOBIAN1 jac, const int mf, void *user_data);

    void DLSODES(ODEPACK_FUNCTION f, const int neq, double *y, double &t, const double tout, 
        const int itol, double *rtol, double *atol, const int itask, int &istate, 
        const int iopt, double *rwork, const int lrw, int *iwork, const int liw, 
        ODEPACK_JACOBIAN2 jac, const int mf, void *user_data);

    void DLSODA(ODEPACK_FUNCTION f, const int neq, double *y, double &t, const double tout,
        const int itol, double *rtol, double *atol, const int itask, int &istate,
        const int iopt, double *rwork, const int lrw, int *iwork, const int liw,
        ODEPACK_JACOBIAN1 jac, const int jt, void *user_data);

    void DLSODAR(ODEPACK_FUNCTION f, const int neq, double *y, double &t, const double tout,
        const int itol, double *rtol, double *atol, const int itask, int &istate,
        const int iopt, double *rwork, const int lrw, int *iwork, const int liw,
        ODEPACK_JACOBIAN1 jac, const int jt, ODEPACK_CONSTRAINT g, const int ng, int *jroot,
        void *user_data);

    void DLSODPK(ODEPACK_FUNCTION f, const int neq, double *y, double &t, const double tout,
        const int itol, double *rtol, double *atol, const int itask, int &istate,
        const int iopt, double *rwork, const int lrw, int *iwork, const int liw,
        ODEPACK_JACOBIAN3 jac, ODEPACK_PSOL psol, const int mf, void *user_data);

    void DLSODKR(ODEPACK_FUNCTION f, const int neq, double *y, double &t, const double tout,
        const int itol, double *rtol, double *atol, const int itask, int &istate,
        const int iopt, double *rwork, const int lrw, int *iwork, const int liw,
        ODEPACK_JACOBIAN3 jac, ODEPACK_PSOL psol, const int mf, ODEPACK_CONSTRAINT g, 
        const int ng, int *jroot, void *user_data);

    void DLSODI(ODEPACK_RESIDUAL res, ODEPACK_ADDA1 adda, ODEPACK_JACOBIAN4 jac, const int neq,
        double *y, double ydoti, double &t, const double tout, const double itol, double *rtol,
        double *atol, const int itask, int &istate, const int iopt, double *Rwork, const int lrw,
        int *iwork, const int liw, const int mf, void *user_data);

    void DLSOIBT(ODEPACK_RESIDUAL res, ODEPACK_ADDA2 adda, ODEPACK_JACOBIAN5 jac, const int neq,
        double *y, double *ydoti, double &T, const double tout, const int itol, const double *rtol,
        const double *atol, const int itask,int &istate, const int iopt, double *rwork, const int lrw,
        int *iwork, const int liw, const int mf, void *user_data);

    void DLSODIS(ODEPACK_RESIDUAL res, ODEPACK_ADDA3 adda, ODEPACK_JACOBIAN6 jac, const int neq,
        double *y, double *ydoti, double &t, const double tout, const int itol, double *rtol,
        double *atol, const int itask, int &istate, const int iopt, double *rwork, const int lrw,
        int *iwork, const int liw, const int mf, void *user_data);

    // opkda2.cpp
    void DGEFA(double *A, const int lda, const int n, int *ipvt, int &info);
    void DGESL(const double *A, const int lda, const int n, int *ipvt, double *b, const int job);
    void DGBFA(double *abd, const int lda, const int n, int ml, int mu, int *ipvt, int &info);
    void DGBSL(double *abd, const int lda, const int n, const int ml, const int mu, int *ipvt, double *b, const int job);
    void DAXPY(const int n, const double da, const double *dx, const int incx, double *dy, const int incy);
    void DCOPY(const int n, const double *dx, const int incx, double *dy, const int incy);
    double DDOT(const int n, const double *dx, const int incx, const double *dy, const int incy);
    double DNRM2(const int n, const double *dx, const int incx);
    void DSCAL(const int n, const double da, double *dx, const int incx);
    int  IDAMAX(const int n, double *dx, const int incx);
    void XERRWD(const std::string msg, const int nmes, const int nerr, const int level, 
        const int ni, const int i1, const int i2, const int nr, const double r1, const double r2);
    int IXSAV(const int ipar, const int ivalue, const bool iset);
    int IUMACH();

private:

    template <typename ODEPACK_JACOBIAN>
    using FUNC_PJAC = void(Odepack::*)(const int neq,  double *y, double *yh, const int nyh, const double *ewt, 
                        double *ftem, double *savf, double *wm, int *iwm, ODEPACK_FUNCTION f, 
                        ODEPACK_JACOBIAN jac, void *user_data);

    using FUNC_SLVS = void(Odepack::*)(double *wm, int *iwm, double *x, double *tem);

    // opkda1.cpp
    double DUMACH();
    void DUMSUM(const double a, const double b, double &c);
    void DCFODE(const int meth, double *elco, double *tesco);
    void DINTDY(const double t, const int k, const double *yh, const int nyh, double *dky, int &iflag);
    void DPREPJ(const int neq,  double *y, double *yh, const int nyh, const double *ewt, 
        double *ftem, double *savf, double *wm, int *iwm, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN1 jac, 
        void *user_data);
    void DSOLSY(double *wm, int *iwm, double *x, double *tem);
    void DSRCOM(double *rsav, int *isav, const int job);

    template <typename ODEPACK_JACOBIAN>
    void DSTODE(const int neq, double *y, double *yh, int nyh, double*yh1, 
        double *ewt, double *savf, double *acor, double *wm, void *iwm_in, 
        ODEPACK_FUNCTION f, ODEPACK_JACOBIAN jac, FUNC_PJAC<ODEPACK_JACOBIAN> pjac, FUNC_SLVS slvs, 
        void *user_data);

    void DEWSET(const int n, const int itol, const double *rtol, const double *atol, const double *ycur,  double *ewt);
    double DVNORM(const int n, const double *v, const double *w);

    void DIPREP(const int neq, double *y, double *rwork, int *ia, int *ja,
        int &ipflag, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN2 jac, void *user_data);
    void DPREP(const int neq, double *y, double *yh, double *savf, double *ewt,
        double *ftem, int *ia, int *ja, double *wk, void *iwk_in, int &ipper, 
        ODEPACK_FUNCTION f, ODEPACK_JACOBIAN2 jac, void *user_data);
    void JGROUP(const int n, int *ia, int *ja, const int maxg, int &ngrp, int *igp, int *jgp, int *incl, int *jdone, int &ier);
    void ADJLR(const int n, const int *isp, int &ldif);
    void CNTNZU(const int n, const int *ia, const int *ja, int &nzsut);
    void DPRJS(const int neq,  double *y, double *yh, const int nyh, const double *ewt, 
        double *ftem, double *savf, double *wk, int *iwk, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN2 jac, 
        void *user_data);
    void DSOLSS(double *wk, int *iwk, double *x, double *tem);
    void DSRCMS(double *rsav, int *isav, const int job);
    void ODRV(const int n, int *ia, int *ja, double *a, int *p, int *ip, int nsp, int *isp, int path, int &flag);
    void MD(const int n, int *ia, int *ja, int max, int *v, int *l, int *head, int *last, int *next, int *mark, int &flag);
    void MDI(const int n, int *ia, int *ja, int max, int *v, int *l, int *head, int *last, int *next, int *mark, int tag, int &flag);
    void MDM(const int vk, int &tail, int *v, int *l, int *last, int *next, int *mark);
    void MDP(int &k, int ek, int &tail, int *v, int *l, int *head, int *last, int *next, int *mark);
    void MDU(int ek, int &dmin, int *v, int *l, int *head, int *last, int *next, int *mark);
    void SRO(const int n, int *ip, int *ia, int *ja, double *a, int *q, int *r, const bool dflag);
    void CDRV(const int n, int *r, int *c, int *ic, int *ia, int *ja, 
        double *a, double *b, double *z, int nsp, int *isp, double *rsp, 
        int &esp, const int path, int &flag);
    void NROC(const int n, int *ic, int *ia, int *ja, double *a, int *jar, 
        double *ar, int *p, int &flag);
    void NSFC(const int n, int *r, const int *ic, const int *ia, const int *ja, 
        const int jlmax, int *il, int *jl, int *ijl, const int jumax, int *iu, int *ju, 
        int *iju, int *q, int *ira, int *jra, int *irac, int *irl, int *jrl, int *iru, 
        int *jru, int &flag);
    void NNFC(const int n, const int *r, const int *c, const int *ic, const int *ia, 
        const int *ja, const double *a, double *z, const double *b, const int lmax, 
        const int *il, const int *jl, const int *ijl, double *l, double *d, const int umax,
        const int *iu, const int *ju, const int *iju, double *u, double *row, double *tmp, 
        int *irl, int *jrl, int &flag);
    void NNSC(const int n, const int *r, const int *c, const int *il, const int *jl, 
        const int *ijl, const double *l, const double *d, const int *iu, const int *ju, 
        const int *iju, const double *u, double *z, const double *b, double *tmp);
    void NNTC(const int n, const int *r, const int *c, const int *il, const int *jl, 
        const int *ijl, const double *l, const double *d, const int *iu, const int *ju, 
        const int *iju, const double *u, double *z, const double *b, double *tmp);

    void DSTODA(const int neq, double *y, double *yh, const int nyh, double *yh1, 
        double *ewt, double *savf, double *acor, double *wm, int *iwm, 
        ODEPACK_FUNCTION f, ODEPACK_JACOBIAN1 jac, FUNC_PJAC<ODEPACK_JACOBIAN1>  pjac, FUNC_SLVS slvs, 
        void *user_data);

    void DPRJA(const int neq,  double *y, double *yh, const int nyh, const double *ewt,  
        double *ftem, double *savf, double *wm, int *iwm, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN1 jac, 
        void *user_data);
    double DMNORM(const int n, const double *v, const double *w);
    double DFNORM(const int n, const double *A, const double *w);
    double DBNORM(const int n, const double *A, const int nra, const int ml, const int mu, const double *w);
    void DSRCMA(double *rsav, int *isav, const int job);
    
    void DRCHEK(const int job, ODEPACK_CONSTRAINT g, const int neq, double *y, double *yh, const int nyh, double *g0, double *g1, double *gx, int *jroot, int &irt, void *user_data);
    void DROOTS(const int ng, const double hmin, int &jflag, double &x0, double &x1, double *g0, double *g1, double *gx, double &x, int *jroot);
    void DSRCAR(double *rsav, int *isav, const int job);

    void DSTODPK(const int neq, double *y,double *yh, const int nyh, double *yh1, double *ewt, double *savf,
        double *savx, double *acor, double *wm, int *iwm, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN3 jac, 
        ODEPACK_PSOL psol, void *user_data);
    void DPKSET(const int neq, double *y, double *ysv, double *ewt, double *ftem, double *savf, double *wm, int *iwm,
        ODEPACK_FUNCTION f, ODEPACK_JACOBIAN3 jac, void *user_data);
    void DSOLPK(const int neq, double *y, double *savf, double *x, double *ewt, double *wm, int *iwm,
        ODEPACK_FUNCTION f, ODEPACK_PSOL psol, void *user_data);
    void DSPIOM(const int neq, double tn, double *y, double *savf, double *b, double *wght, const int n, const int maxl, int kmp,
        double &delta, const double hl0, int jpre, int &mnewt, ODEPACK_FUNCTION f, ODEPACK_PSOL psol, int &npsl, double *x,
        double *v, double *hes, int *ipvt, int &liom, double *wp, int *iwp, double *wk, int &iflag, void *user_data);
    void DATV(const int neq, double *y, double *savf, double *v, double *wght, double *ftem, ODEPACK_FUNCTION f,
        ODEPACK_PSOL psol, double *z, double *vtem, double *wp, int *iwp, double hl0, int &jpre, int &ier, int &npsl,
        void *user_data);
    void DORTHOG(double *vnew, double *v, double *hes, const int n, const int ll, const int ldhes, const int kmp, double &snormw);

    void DSPIGMR(const int neq, double tn, double *y, double *savf, double *b, double *wght, int n, const int maxl,
        int maxlp1, int kmp, double &delta, double hl0, int jpre, const int mnewt, ODEPACK_FUNCTION f, ODEPACK_PSOL psol,
        int &npsl, double *x, double *v, double *hes, double *q, int &lgmr, double *wp, int *iwp, double *wk, double *dl,
        int &iflag, void *user_data);

    void DPCG(const int neq, double tn, double *y, double *savf, double *r, double *wght, int n, const int maxl,
        const double delta, double hl0, const int &jpre, const int mnewt, ODEPACK_FUNCTION f, ODEPACK_PSOL psol,
        int &npsl, double *x, double *p, double *w, double *z, int &lpcg, double *wp,int *iwp, double *wk, int &iflag,
        void *user_data);

    void DPCGS(const int neq, double tn, double *y, double *savf, double *r, double *wght, int n, const int maxl,
        const double delta, double hl0, const int jpre, const int mnewt, ODEPACK_FUNCTION f, ODEPACK_PSOL psol,
        int &npsl, double *x, double *p, double *w, double *z, int &lpcg, double *wp, int *iwp, double *wk, int &iflag,
        void *user_data);

    void DATP(const int neq, double *y, const double *savf, double *p, double *wght, const double hl0, double *wk,
        ODEPACK_FUNCTION f, double *w, void *user_data);

    void DUSOL(const int neq, double tn, double *y, double *savf, double *b, double *wght, int n, const double delta,
        double hl0, const int mnewt, ODEPACK_PSOL psol, int &npsl, double *x, double *wp, int *iwp, double *wk, int &iflag,
        void *user_data);
        
    void DSRCPK(double *rsav, int *isav, const int job);
    void DHEFA(double *a, const int lda, const int n, int *ipvt, int &info, const int job);
    void DHESL(double *a, const int lda, const int n, const int *ipvt, double *b);
    void DHEQR(double *a, const int lda, const int n, double *q, int &info, const int ijob);
    void DHELS(double *A, const int lda, const int n, const double *q, double *b);
    void DLHIN(const int neq, int n, const double t0, double *y0, const double *ydot, ODEPACK_FUNCTION f, const double tout,
        const double *uround, double *ewt, const int itol, const double *atol, double *y, double *temp, double &h0, 
        int &niter, int &ier);

    void DSTOKA(const int neq, double *y, double *yh, const int nyh, double *yh1, double *ewt, double *savf, double *savx,
        double *acor, double *wm, int *iwm, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN1 jac, ODEPACK_PSOL psol, void *user_data);
    void DSETPK(const int neq, double *y, double *ysv, double *ewt, double *ftem, double *savf, int jok, double *wm,int *iwm, 
        ODEPACK_FUNCTION f, ODEPACK_JACOBIAN2 jac, void *user_data);
    void DSRCKR(double *rsav, int*isav, const int job);
    void DAINVG(ODEPACK_RESIDUAL res, ODEPACK_ADDA1 adda, const int neq, double t, double *y, double *ydot, int &miter, int &ml,
        int &mu, double *pw, int *ipvt, int &ier, void *user_data);
    void DSTODI(const int neq, double *y, double *yh, const int nyh, double *yh1, double *ewt, double *savf, double *savr,
        double *acor, double *wm, int *iwm, ODEPACK_RESIDUAL res, ODEPACK_ADDA1 adda, ODEPACK_JACOBIAN1 jac, 
        FUNC_PJAC<ODEPACK_JACOBIAN1> pjac, FUNC_SLVS slvs, void *user_data);
    void DPREPJI(const int neq, double *y, double *yh, const int nyh, double *ewt, double *rtem, double *savr, double *s,
        double *wm, int *iwm, ODEPACK_RESIDUAL res, ODEPACK_JACOBIAN1 jac, ODEPACK_ADDA1 adda, void *user_data);
    void DAIGBT(ODEPACK_RESIDUAL res, ODEPACK_ADDA1 adda, const int neq, double t, double *y, double *ydot, int &mb, int &nb,
        double *pw, int *ipvt, int &ier, void *user_data);
    void DPJIBT(const int neq, double *y, double *yh, const int nyh, double *ewt, double *rtem, double *savr, double *s,
        double *wm, int *iwm, ODEPACK_RESIDUAL res, ODEPACK_JACOBIAN1 jac, ODEPACK_ADDA1 adda, void *user_data);
    void DSLSBT(double *wm, int *iwm, double *x, double *tem);
    void DDECBT(int m, const int n, double *a, double *b, double *c, int *ip, int &ier);
    void DSOLBT(int m, const int n, double *a, double *b, double *c, double *y, int *ip);
    void DIPREPI(const int neq, double *y, double *s, double *rwork, int *ia, int *ja, int *ic, int *jc, int &ipflag,
        ODEPACK_RESIDUAL res, ODEPACK_JACOBIAN1 jac, ODEPACK_ADDA1 adda);
    void DPREPI(const int neq, double *y, double *s, const double *yh, double *savr, const double* ewt, double *rtem,
        const int *ia, const int *ja, const int *ic, const int *jc, double *wk, int *iwk, int &ipper,
        ODEPACK_RESIDUAL res, ODEPACK_JACOBIAN1 jac, ODEPACK_ADDA1 adda, void *user_data);
    void DAINVGS(const int neq, double t, double *y, double *wk, int *iwk, double *tem, double *ydot,int &ier,
        ODEPACK_RESIDUAL res, ODEPACK_ADDA1 adda, void *user_data);
    void DPRJIS(const int neq, double *y, double *yh, const int nyh, const double *ewt, double *rtem, double *savr,
        double *s, double *wk, int *iwk, ODEPACK_RESIDUAL res, ODEPACK_JACOBIAN1 jac, ODEPACK_ADDA1 adda, void *user_data);

    struct DLS001 {
        double coint, crate, el[13], elco[13*12], hold, rmax, tesco[3*12];
        double ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround;
        int    init, mxstep, mxhnil, nhnil, nslast, nyh, ialth, ipup, lmax, meo;
        int    nqnyh, nslp, icf, ierpj, iersl, jcur, jstart, kflag, l;
        int    lyh, lewt, lacor, lsavf, lwm, liwm;
        int    meth, miter, maxord, maxcor, msbp;
        int    mxncf, n, nq, nst, nfe, nje, nqu;
        double rls[218], rlss[6];
        int    ils[37], ilss[34];
    } dls1_;

    struct DLSS01 {
        double con0, conmin, ccmxj, psmall, rbig, seth;
        int    iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp, ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp;
        int    ipa, lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj, nslj, ngp, nlu, nnz, nsp, nzl, nzu;
    } dlss_;

    struct DLSA01 {
        double tsw, cm1[12], cm2[5], pdest, pdlast, ratio, pdnorm;
        int    insufr, insufi, ixpr, icount, irflag, jtyp, mused, mxordn, mxords;
    } dlsa_;

    struct DLSR01 {
        double t0, tlast, toutc;
        int    lg0, lg1, lgx, irfnd, itaskc, ngc, nge;
        double alpha, x2;
        int    imax, last;
    } dlsr_;

    struct DLPK01 {
        double delt, epcon, sqrtn, rsqrtn;
        int    jpre, jacflg, locwp, lociwp, lsavx, kmp, maxl, mnewt, nni, nli, nps, ncfn, ncfl;
    } dlpk_;

};

}

#endif /* ODEPACK_HPP_ */