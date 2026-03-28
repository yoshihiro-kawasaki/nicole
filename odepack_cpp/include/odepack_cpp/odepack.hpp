#ifndef ODEPACK_HPP_
#define ODEPACK_HPP_

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

/**
* @def ARRAY1D
* @brief
*/
#ifndef ARRAY1D
#define ARRAY1D(a, i) (a[(i)-1])
#endif


/**
 * @def ARRAY2D
 * @brief Column major 2D array
*/
#ifndef ARRAY2D
#define ARRAY2D(a, rows, i, j) (a[(i-1) + (j-1) * (rows)])
#endif


/**
 * @def ARRAY3D
 * @brief Column major 3D array
 */
#ifndef ARRAY3D
#define ARRAY3D(a, n1, n2, i, j, k) (a[(i-1) + n1 * ((j-1) + n2 * (k-1))])
#endif

/**
* @namespace odepack_cpp
*/
namespace odepack_cpp {

using odepack_cpp_real = double;

/**
 * @fn SIGN
*/
template <typename T>
constexpr T SIGN(T x) {
    return (x >= static_cast<T>(0) ? static_cast<T>(1) : static_cast<T>(-1));
}


/**
 * ODEPACK_FUNCTION
 * DLSODE, DLSODES, DLSODA, DLSODAR, DLSODPK, DLSODKR, DLSODI, DLSOIBT, DLSODIS
*/
using ODEPACK_FUNCTION = void(*)(
    int neq, odepack_cpp_real t, odepack_cpp_real *y, odepack_cpp_real *ydot, void *user_data
);

/**
 * ODEPACK_JACOBIAN1
 * DLSODE, DLSODA, DLSODAR
 */
using ODEPACK_JACOBIAN1 = void(*)(
    int neq, odepack_cpp_real t, odepack_cpp_real *y, int ml, int mu, odepack_cpp_real *pd, int nrowpd, void *user_data
);

/**
 * ODEPACK_JACOBIAN2
 * DLSODES,
 */
using ODEPACK_JACOBIAN2 = void(*)(
    int neq, odepack_cpp_real t, odepack_cpp_real *y, int j, int *ian, int *jan, odepack_cpp_real *pdj, void *user_data
);

/**
 * ODEPACK_JACOBIAN3
 * DLSODPK, DLSODKR,
 */
using ODEPACK_JACOBIAN3 = void(*)(
    ODEPACK_FUNCTION f, int neq, odepack_cpp_real t, odepack_cpp_real *y, odepack_cpp_real *ysv, odepack_cpp_real *rewt, odepack_cpp_real *fty,
    odepack_cpp_real *v, odepack_cpp_real hl0, odepack_cpp_real *wp, int *iwp, int ier, void *user_data
);

/**
 * ODEPACK_JACOBIAN4
 * DLSODI
 */
using ODEPACK_JACOBIAN4  = void(*)(
    int neq, odepack_cpp_real t, odepack_cpp_real *y, odepack_cpp_real *s, int ml, int mu,
    odepack_cpp_real *p, int nrowp, void *user_data
);

/**
 * ODEPACK_JACOBIAN5
 * DLSOIBT
 */
using ODEPACK_JACOBIAN5  = void(*)(
    int neq, odepack_cpp_real t, odepack_cpp_real *y, odepack_cpp_real *s, int mb, int nb, odepack_cpp_real *pa, odepack_cpp_real *pb, odepack_cpp_real *pc, void *user_data
);

/**
 * ODEPACK_JACOBIAN6
 * DLSODIS
 */
using ODEPACK_JACOBIAN6  = void(*)(
    int neq, odepack_cpp_real t, odepack_cpp_real *y, odepack_cpp_real *s, int j, int *ian, int *jan, odepack_cpp_real *pdj, void *user_data
);

/**
 * ODEPACK_CONSTRAINT
 * DLSODAR, DLSODKR,
 */
using ODEPACK_CONSTRAINT = void(*)(
    int neq, odepack_cpp_real t, odepack_cpp_real *y, int ng, odepack_cpp_real *gout, void *user_data
);

/**
 * ODEPACK_PSOL
 */
using ODEPACK_PSOL = void(*)(
    int neq, odepack_cpp_real t, odepack_cpp_real *y, odepack_cpp_real *fty, odepack_cpp_real *wk, odepack_cpp_real hl0,
    odepack_cpp_real *wp, int *iwp, odepack_cpp_real *b, int lr, int ier, void *user_data
);

/**
 * ODEPACK_RESIDUAL
 * DLSODPK, DLSODI, DLSOIBT, DLSODIS
 */
using ODEPACK_RESIDUAL = void(*)(
    int neq, odepack_cpp_real t, odepack_cpp_real *y, odepack_cpp_real *s, odepack_cpp_real *r, int ires, void *user_data
);

/**
 * ODEPACK_ADDA1
 * DLSODI,
 */
using ODEPACK_ADDA1 = void(*)(
    int neq, odepack_cpp_real t, odepack_cpp_real *y, int ml, int mu, odepack_cpp_real *p, int nrowp, void *user_data
);

/**
 * ODEPACK_ADDA2
 * DLSOIBT
 */
using ODEPACK_ADDA2 = void(*)(
    int neq, odepack_cpp_real t, odepack_cpp_real *y, int mb, int nb, odepack_cpp_real *pa, odepack_cpp_real *pb, odepack_cpp_real *pc, void *user_data
);

/**
 * ODEPACK_ADDA3
 * DLSODIS
 */
using ODEPACK_ADDA3 = void(*)(
    int neq, odepack_cpp_real t, odepack_cpp_real *y, int j, int *ian, int *jan, odepack_cpp_real *p, void *user_data
);


/**
 * @fn class Odepack
*/
class Odepack {
public:
    Odepack() { };
    ~Odepack() { };

    // opkdmain.cpp
    void DLSODE(
        ODEPACK_FUNCTION f, int neq, odepack_cpp_real *y, odepack_cpp_real &t, odepack_cpp_real tout, int itol,
        odepack_cpp_real *rtol, odepack_cpp_real *atol, int itask, int &istate, int iopt, odepack_cpp_real *rwork,
        int lrw, int *iwork, int liw, ODEPACK_JACOBIAN1 jac, int mf, void *user_data
    );

    void DLSODES(
        ODEPACK_FUNCTION f, int neq, odepack_cpp_real *y, odepack_cpp_real &t, odepack_cpp_real tout, int itol,
        odepack_cpp_real *rtol, odepack_cpp_real *atol, int itask, int &istate, int iopt, odepack_cpp_real *rwork,
        int lrw, int *iwork, int liw, ODEPACK_JACOBIAN2 jac, int mf, void *user_data
    );

    void DLSODA(
        ODEPACK_FUNCTION f, int neq, odepack_cpp_real *y, odepack_cpp_real &t, odepack_cpp_real tout, int itol,
        odepack_cpp_real *rtol, odepack_cpp_real *atol, int itask, int &istate, int iopt, odepack_cpp_real *rwork,
        int lrw, int *iwork, int liw, ODEPACK_JACOBIAN1 jac, int jt, void *user_data
    );

    void DLSODAR(
        ODEPACK_FUNCTION f, int neq, odepack_cpp_real *y, odepack_cpp_real &t, odepack_cpp_real tout,
        int itol, odepack_cpp_real *rtol, odepack_cpp_real *atol, int itask, int &istate,
        int iopt, odepack_cpp_real *rwork, int lrw, int *iwork, int liw,
        ODEPACK_JACOBIAN1 jac, int jt, ODEPACK_CONSTRAINT g, int ng, int *jroot,
        void *user_data
    );

    void DLSODPK(
        ODEPACK_FUNCTION f, int neq, odepack_cpp_real *y, odepack_cpp_real &t, odepack_cpp_real tout,
        int itol, odepack_cpp_real *rtol, odepack_cpp_real *atol, int itask, int &istate,
        int iopt, odepack_cpp_real *rwork, int lrw, int *iwork, int liw,
        ODEPACK_JACOBIAN3 jac, ODEPACK_PSOL psol, int mf, void *user_data
    );

    void DLSODKR(
        ODEPACK_FUNCTION f, int neq, odepack_cpp_real *y, odepack_cpp_real &t, odepack_cpp_real tout,
        int itol, odepack_cpp_real *rtol, odepack_cpp_real *atol, int itask, int &istate,
        int iopt, odepack_cpp_real *rwork, int lrw, int *iwork, int liw,
        ODEPACK_JACOBIAN3 jac, ODEPACK_PSOL psol, int mf, ODEPACK_CONSTRAINT g, 
        int ng, int *jroot, void *user_data
    );

    void DLSODI(
        ODEPACK_RESIDUAL res, ODEPACK_ADDA1 adda, ODEPACK_JACOBIAN4 jac, int neq,
        odepack_cpp_real *y, odepack_cpp_real *ydoti, odepack_cpp_real &t, odepack_cpp_real tout, odepack_cpp_real itol, odepack_cpp_real *rtol,
        odepack_cpp_real *atol, int itask, int &istate, int iopt, odepack_cpp_real *Rwork, int lrw,
        int *iwork, int liw, int mf, void *user_data
    );

    void DLSOIBT(
        ODEPACK_RESIDUAL res, ODEPACK_ADDA2 adda, ODEPACK_JACOBIAN5 jac, int neq,
        odepack_cpp_real *y, odepack_cpp_real *ydoti, odepack_cpp_real &T, odepack_cpp_real tout, int itol, odepack_cpp_real *rtol,
        odepack_cpp_real *atol, int itask,int &istate, int iopt, odepack_cpp_real *rwork, int lrw,
        int *iwork, int liw, int mf, void *user_data
    );

    void DLSODIS(
        ODEPACK_RESIDUAL res, ODEPACK_ADDA3 adda, ODEPACK_JACOBIAN6 jac, int neq,
        odepack_cpp_real *y, odepack_cpp_real *ydoti, odepack_cpp_real &t, odepack_cpp_real tout, int itol, odepack_cpp_real *rtol,
        odepack_cpp_real *atol, int itask, int &istate, int iopt, odepack_cpp_real *rwork, int lrw,
        int *iwork, int liw, int mf, void *user_data
    );

private:

    // DSTODE, DSTODA
    template <typename ODEPACK_JACOBIAN>
    using FUNC_PJAC = void(Odepack::*)(
        int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *ewt, odepack_cpp_real *ftem, odepack_cpp_real *savf, odepack_cpp_real *wm, int *iwm, 
        ODEPACK_FUNCTION f, ODEPACK_JACOBIAN jac, void *user_data
    );
    
    // DSTODI
    template <typename ODEPACK_JACOBIAN, typename ODEPACK_ADDA>
    using FUNC_PJAC2 = void(Odepack::*)(
        int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *ewt, 
        odepack_cpp_real *ftem, odepack_cpp_real *savr, odepack_cpp_real *s, odepack_cpp_real *wm, int *iwm, 
        ODEPACK_RESIDUAL res, ODEPACK_JACOBIAN jac, ODEPACK_ADDA adda, void *user_data
    );

    using FUNC_SLVS = void(Odepack::*)(
        odepack_cpp_real *wm, int *iwm, odepack_cpp_real *x, odepack_cpp_real *tem
    );

    // opkda1.cpp
    odepack_cpp_real DUMACH();

    void DUMSUM(odepack_cpp_real a, odepack_cpp_real b, odepack_cpp_real &c);

    void DCFODE(int meth, odepack_cpp_real *elco, odepack_cpp_real *tesco);

    void DINTDY(odepack_cpp_real t, int k, odepack_cpp_real *yh, int nyh, odepack_cpp_real *dky, int &iflag);

    void DPREPJ(int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *ewt,
        odepack_cpp_real *ftem, odepack_cpp_real *savf, odepack_cpp_real *wm, int *iwm, ODEPACK_FUNCTION f,
        ODEPACK_JACOBIAN1 jac, void *user_data
    );

    void DSOLSY(odepack_cpp_real *wm, int *iwm, odepack_cpp_real *x, odepack_cpp_real *tem);

    void DSRCOM(odepack_cpp_real *rsav, int *isav, int job);

    template <typename ODEPACK_JACOBIAN>
    void DSTODE(int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *yh1, odepack_cpp_real *ewt, odepack_cpp_real *savf, odepack_cpp_real *acor,
        odepack_cpp_real *wm, void *iwm_in, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN jac, FUNC_PJAC<ODEPACK_JACOBIAN> pjac,
        FUNC_SLVS slvs, void *user_data
    );

    void DEWSET(
        int n, int itol, odepack_cpp_real *rtol, odepack_cpp_real *atol, odepack_cpp_real *ycur,  odepack_cpp_real *ewt
    );
    
    odepack_cpp_real DVNORM(
        int n, odepack_cpp_real *v, odepack_cpp_real *w
    );

    void DIPREP(
        int neq, odepack_cpp_real *y, odepack_cpp_real *rwork, int *ia, int *ja, int &ipflag, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN2 jac, void *user_data
    );

    void DPREP(
        int neq, odepack_cpp_real *y, odepack_cpp_real *yh, odepack_cpp_real *savf, odepack_cpp_real *ewt, odepack_cpp_real *ftem, int *ia, int *ja,
        odepack_cpp_real *wk, void *iwk_in, int &ipper, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN2 jac, void *user_data
    );

    void JGROUP(
        int n, int *ia, int *ja, int maxg, int &ngrp, int *igp, int *jgp, int *incl, int *jdone, int &ier
    );

    void ADJLR(
        int n, int *isp, int &ldif
    );

    void CNTNZU(
        int n, int *ia, int *ja, int &nzsut
    );

    void DPRJS(
        int neq,  odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *ewt, odepack_cpp_real *ftem, odepack_cpp_real *savf, odepack_cpp_real *wk, int *iwk, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN2 jac, void *user_data
    );

    void DSOLSS(
        odepack_cpp_real *wk, int *iwk, odepack_cpp_real *x, odepack_cpp_real *tem
    );

    void DSRCMS(
        odepack_cpp_real *rsav, int *isav, int job
    );

    void ODRV(
        int n, int *ia, int *ja, odepack_cpp_real *a, int *p, int *ip, int nsp, int *isp, int path, int &flag
    );

    void MD(
        int n, int *ia, int *ja, int max, int *v, int *l, int *head, int *last, int *next, int *mark, int &flag
    );

    void MDI(
        int n, int *ia, int *ja, int max, int *v, int *l, int *head, int *last, int *next, int *mark, int tag, int &flag
    );

    void MDM(
        int vk, int &tail, int *v, int *l, int *last, int *next, int *mark
    );

    void MDP(
        int &k, int ek, int &tail, int *v, int *l, int *head, int *last, int *next, int *mark
    );

    void MDU(
        int ek, int &dmin, int *v, int *l, int *head, int *last, int *next, int *mark
    );

    void SRO(
        int n, int *ip, int *ia, int *ja, odepack_cpp_real *a, int *q, int *r, bool dflag
    );

    void CDRV(
        int n, int *r, int *c, int *ic, int *ia, int *ja, odepack_cpp_real *a, odepack_cpp_real *b, odepack_cpp_real *z, int nsp, int *isp, odepack_cpp_real *rsp, int &esp, int path, int &flag
    );

    void NROC(
        int n, int *ic, int *ia, int *ja, odepack_cpp_real *a, int *jar, odepack_cpp_real *ar, int *p, int &flag
    );

    void NSFC(
        int n, int *r, int *ic, int *ia, int *ja, int jlmax, int *il, int *jl, int *ijl, int jumax, int *iu, int *ju, int *iju, int *q,
        int *ira, int *jra, int *irac, int *irl, int *jrl, int *iru, int *jru, int &flag
    );

    void NNFC(
        int n, int *r, int *c, int *ic, int *ia, int *ja, odepack_cpp_real *a, odepack_cpp_real *z, odepack_cpp_real *b,
        int lmax, int *il, int *jl, int *ijl, odepack_cpp_real *l, odepack_cpp_real *d, int umax, int *iu, int *ju, int *iju, odepack_cpp_real *u,
        odepack_cpp_real *row, odepack_cpp_real *tmp, int *irl, int *jrl, int &flag
    );

    void NNSC(
        int n, int *r, int *c, int *il, int *jl, int *ijl, odepack_cpp_real *l, odepack_cpp_real *d, int *iu, int *ju, int *iju, odepack_cpp_real *u, odepack_cpp_real *z, odepack_cpp_real *b, odepack_cpp_real *tmp
    );

    void NNTC(
        int n, int *r, int *c, int *il, int *jl, 
        int *ijl, odepack_cpp_real *l, odepack_cpp_real *d, int *iu, int *ju, 
        int *iju, odepack_cpp_real *u, odepack_cpp_real *z, odepack_cpp_real *b, odepack_cpp_real *tmp
    );

    void DSTODA(
        int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *yh1, 
        odepack_cpp_real *ewt, odepack_cpp_real *savf, odepack_cpp_real *acor, odepack_cpp_real *wm, int *iwm, 
        ODEPACK_FUNCTION f, ODEPACK_JACOBIAN1 jac, FUNC_PJAC<ODEPACK_JACOBIAN1>  pjac, FUNC_SLVS slvs, 
        void *user_data
    );

    void DPRJA(
        int neq,  odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *ewt,  
        odepack_cpp_real *ftem, odepack_cpp_real *savf, odepack_cpp_real *wm, int *iwm, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN1 jac, 
        void *user_data
    );

    odepack_cpp_real DMNORM(
        int n, odepack_cpp_real *v, odepack_cpp_real *w
    );

    odepack_cpp_real DFNORM(
        int n, odepack_cpp_real *A, odepack_cpp_real *w
    );

    odepack_cpp_real DBNORM(
        int n, odepack_cpp_real *A, int nra, int ml, int mu, odepack_cpp_real *w
    );

    void DSRCMA(
        odepack_cpp_real *rsav, int *isav, int job
    );
    
    void DRCHEK(
        int job, ODEPACK_CONSTRAINT g, int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *g0, odepack_cpp_real *g1, odepack_cpp_real *gx, int *jroot, int &irt, void *user_data
    );

    void DROOTS(
        int ng, odepack_cpp_real hmin, int &jflag, odepack_cpp_real &x0, odepack_cpp_real &x1, odepack_cpp_real *g0, odepack_cpp_real *g1, odepack_cpp_real *gx, odepack_cpp_real &x, int *jroot
    );

    void DSRCAR(
        odepack_cpp_real *rsav, int *isav, int job
    );

    void DSTODPK(
        int neq, odepack_cpp_real *y,odepack_cpp_real *yh, int nyh, odepack_cpp_real *yh1, odepack_cpp_real *ewt, odepack_cpp_real *savf,
        odepack_cpp_real *savx, odepack_cpp_real *acor, odepack_cpp_real *wm, int *iwm, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN3 jac, 
        ODEPACK_PSOL psol, void *user_data
    );

    void DPKSET(
        int neq, odepack_cpp_real *y, odepack_cpp_real *ysv, odepack_cpp_real *ewt, odepack_cpp_real *ftem, odepack_cpp_real *savf, odepack_cpp_real *wm, int *iwm,
        ODEPACK_FUNCTION f, ODEPACK_JACOBIAN3 jac, void *user_data
    );

    void DSOLPK(
        int neq, odepack_cpp_real *y, odepack_cpp_real *savf, odepack_cpp_real *x, odepack_cpp_real *ewt, odepack_cpp_real *wm, int *iwm,
        ODEPACK_FUNCTION f, ODEPACK_PSOL psol, void *user_data
    );

    void DSPIOM(
        int neq, odepack_cpp_real tn, odepack_cpp_real *y, odepack_cpp_real *savf, odepack_cpp_real *b, odepack_cpp_real *wght, int n, int maxl, int kmp,
        odepack_cpp_real &delta, odepack_cpp_real hl0, int jpre, int &mnewt, ODEPACK_FUNCTION f, ODEPACK_PSOL psol, int &npsl, odepack_cpp_real *x,
        odepack_cpp_real *v, odepack_cpp_real *hes, int *ipvt, int &liom, odepack_cpp_real *wp, int *iwp, odepack_cpp_real *wk, int &iflag, void *user_data
    );

    void DATV(
        int neq, odepack_cpp_real *y, odepack_cpp_real *savf, odepack_cpp_real *v, odepack_cpp_real *wght, odepack_cpp_real *ftem, ODEPACK_FUNCTION f,
        ODEPACK_PSOL psol, odepack_cpp_real *z, odepack_cpp_real *vtem, odepack_cpp_real *wp, int *iwp, odepack_cpp_real hl0, int &jpre, int &ier, int &npsl,
        void *user_data
    );

    void DORTHOG(
        odepack_cpp_real *vnew, odepack_cpp_real *v, odepack_cpp_real *hes, int n, int ll, int ldhes, int kmp, odepack_cpp_real &snormw
    );

    void DSPIGMR(
        int neq, odepack_cpp_real tn, odepack_cpp_real *y, odepack_cpp_real *savf, odepack_cpp_real *b, odepack_cpp_real *wght, int n, int maxl,
        int maxlp1, int kmp, odepack_cpp_real &delta, odepack_cpp_real hl0, int jpre, int mnewt, ODEPACK_FUNCTION f, ODEPACK_PSOL psol,
        int &npsl, odepack_cpp_real *x, odepack_cpp_real *v, odepack_cpp_real *hes, odepack_cpp_real *q, int &lgmr, odepack_cpp_real *wp, int *iwp, odepack_cpp_real *wk, odepack_cpp_real *dl,
        int &iflag, void *user_data
    );

    void DPCG(
        int neq, odepack_cpp_real tn, odepack_cpp_real *y, odepack_cpp_real *savf, odepack_cpp_real *r, odepack_cpp_real *wght, int n, int maxl,
        odepack_cpp_real delta, odepack_cpp_real hl0, int &jpre, int mnewt, ODEPACK_FUNCTION f, ODEPACK_PSOL psol,
        int &npsl, odepack_cpp_real *x, odepack_cpp_real *p, odepack_cpp_real *w, odepack_cpp_real *z, int &lpcg, odepack_cpp_real *wp,int *iwp, odepack_cpp_real *wk, int &iflag,
        void *user_data
    );

    void DPCGS(
        int neq, odepack_cpp_real tn, odepack_cpp_real *y, odepack_cpp_real *savf, odepack_cpp_real *r, odepack_cpp_real *wght, int n, int maxl,
        odepack_cpp_real delta, odepack_cpp_real hl0, int jpre, int mnewt, ODEPACK_FUNCTION f, ODEPACK_PSOL psol,
        int &npsl, odepack_cpp_real *x, odepack_cpp_real *p, odepack_cpp_real *w, odepack_cpp_real *z, int &lpcg, odepack_cpp_real *wp, int *iwp, odepack_cpp_real *wk, int &iflag,
        void *user_data
    );

    void DATP(
        int neq, odepack_cpp_real *y, odepack_cpp_real *savf, odepack_cpp_real *p, odepack_cpp_real *wght, odepack_cpp_real hl0, odepack_cpp_real *wk,
        ODEPACK_FUNCTION f, odepack_cpp_real *w, void *user_data
    );

    void DUSOL(
        int neq, odepack_cpp_real tn, odepack_cpp_real *y, odepack_cpp_real *savf, odepack_cpp_real *b, odepack_cpp_real *wght, int n, odepack_cpp_real delta,
        odepack_cpp_real hl0, int mnewt, ODEPACK_PSOL psol, int &npsl, odepack_cpp_real *x, odepack_cpp_real *wp, int *iwp, odepack_cpp_real *wk, int &iflag,
        void *user_data
    );
        
    void DSRCPK(
        odepack_cpp_real *rsav, int *isav, int job
    );

    void DHEFA(
        odepack_cpp_real *a, int lda, int n, int *ipvt, int &info, int job
    );

    void DHESL(
        odepack_cpp_real *a, int lda, int n, int *ipvt, odepack_cpp_real *b
    );

    void DHEQR(
        odepack_cpp_real *a, int lda, int n, odepack_cpp_real *q, int &info, int ijob
    );

    void DHELS(
        odepack_cpp_real *A, int lda, int n, odepack_cpp_real *q, odepack_cpp_real *b
    );

    void DLHIN(
        int neq, int n, odepack_cpp_real t0, odepack_cpp_real *y0, odepack_cpp_real *ydot, ODEPACK_FUNCTION f, odepack_cpp_real tout,
        odepack_cpp_real uround, odepack_cpp_real *ewt, int itol, odepack_cpp_real *atol, odepack_cpp_real *y, odepack_cpp_real *temp, odepack_cpp_real &h0, 
        int &niter, int &ier, void *user_data
    );

    void DSTOKA(
        int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *yh1, odepack_cpp_real *ewt, odepack_cpp_real *savf, odepack_cpp_real *savx,
        odepack_cpp_real *acor, odepack_cpp_real *wm, int *iwm, ODEPACK_FUNCTION f, ODEPACK_JACOBIAN3 jac, ODEPACK_PSOL psol, void *user_data
    );

    void DSETPK(
        int neq, odepack_cpp_real *y, odepack_cpp_real *ysv, odepack_cpp_real *ewt, odepack_cpp_real *ftem, odepack_cpp_real *savf, int jok, odepack_cpp_real *wm,int *iwm, 
        ODEPACK_FUNCTION f, ODEPACK_JACOBIAN3 jac, void *user_data
    );

    void DSRCKR(
        odepack_cpp_real *rsav, int*isav, int job
    );

    void DAINVG(
        ODEPACK_RESIDUAL res, ODEPACK_ADDA1 adda, int neq, odepack_cpp_real t, odepack_cpp_real *y, odepack_cpp_real *ydot, int &miter, int &ml,
        int &mu, odepack_cpp_real *pw, int *ipvt, int &ier, void *user_data
    );
    
    template <typename ODEPACK_ADDA, typename ODEPACK_JACOBIAN>
    void DSTODI(
        int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *yh1, odepack_cpp_real *ewt, odepack_cpp_real *savf, odepack_cpp_real *savr,
        odepack_cpp_real *acor, odepack_cpp_real *wm, void *iwm_in, ODEPACK_RESIDUAL res, ODEPACK_ADDA adda, ODEPACK_JACOBIAN jac, 
        FUNC_PJAC2<ODEPACK_JACOBIAN, ODEPACK_ADDA> pjac, FUNC_SLVS slvs, void *user_data
    );

    void DPREPJI(
        int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *ewt, odepack_cpp_real *rtem, odepack_cpp_real *savr, odepack_cpp_real *s,
        odepack_cpp_real *wm, int *iwm, ODEPACK_RESIDUAL res, ODEPACK_JACOBIAN4 jac, ODEPACK_ADDA1 adda, void *user_data
    );

    void DAIGBT(
        ODEPACK_RESIDUAL res, ODEPACK_ADDA2 adda, int neq, odepack_cpp_real t, odepack_cpp_real *y, odepack_cpp_real *ydot, int &mb, int &nb,
        odepack_cpp_real *pw, int *ipvt, int &ier, void *user_data
    );

    void DPJIBT(
        int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *ewt, odepack_cpp_real *rtem, odepack_cpp_real *savr, odepack_cpp_real *s,
        odepack_cpp_real *wm, int *iwm, ODEPACK_RESIDUAL res, ODEPACK_JACOBIAN5 jac, ODEPACK_ADDA2 adda, void *user_data
    );

    void DSLSBT(
        odepack_cpp_real *wm, int *iwm, odepack_cpp_real *x, odepack_cpp_real *tem
    );

    void DDECBT(
        int m, int n, odepack_cpp_real *a, odepack_cpp_real *b, odepack_cpp_real *c, int *ip, int &ier
    );

    void DSOLBT(
        int m, int n, odepack_cpp_real *a, odepack_cpp_real *b, odepack_cpp_real *c, odepack_cpp_real *y, int *ip
    );
    
    void DIPREPI(
        int neq, odepack_cpp_real *y, odepack_cpp_real *s, odepack_cpp_real *rwork, int *ia, int *ja, int *ic, int *jc, int &ipflag,
        ODEPACK_RESIDUAL res, ODEPACK_JACOBIAN6 jac, ODEPACK_ADDA3 adda, void *user_data
    );

    void DPREPI(
        int neq, odepack_cpp_real *y, odepack_cpp_real *s, odepack_cpp_real *yh, odepack_cpp_real *savr, odepack_cpp_real* ewt, odepack_cpp_real *rtem,
        int *ia, int *ja, int *ic, int *jc, odepack_cpp_real *wk, void *iwk_in, int &ipper,
        ODEPACK_RESIDUAL res, ODEPACK_JACOBIAN6 jac, ODEPACK_ADDA3 adda, void *user_data
    );
    
    void DAINVGS(
        int neq, odepack_cpp_real t, odepack_cpp_real *y, odepack_cpp_real *wk, void *iwk_in, odepack_cpp_real *tem, odepack_cpp_real *ydot,int &ier,
        ODEPACK_RESIDUAL res, ODEPACK_ADDA3 adda, void *user_data
    );
    
    void DPRJIS(
        int neq, odepack_cpp_real *y, odepack_cpp_real *yh, int nyh, odepack_cpp_real *ewt, odepack_cpp_real *rtem, odepack_cpp_real *savr, odepack_cpp_real *s,
        odepack_cpp_real *wk, int *iwk, ODEPACK_RESIDUAL res, ODEPACK_JACOBIAN6 jac, ODEPACK_ADDA3 adda, void *user_data
    );

    struct DLS001 {
        odepack_cpp_real conit, crate, el[13], elco[13*12], hold, rmax, tesco[3*12];
        odepack_cpp_real ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround;
        int    init, mxstep, mxhnil, nhnil, nslast, nyh, ialth, ipup, lmax, meo;
        int    nqnyh, nslp, icf, ierpj, iersl, jcur, jstart, kflag, l;
        int    lyh, lewt, lacor, lsavf, lwm, liwm;
        int    meth, miter, maxord, maxcor, msbp;
        int    mxncf, n, nq, nst, nfe, nje, nqu;
        static const int lenrls = 218, lenils = 37;
        odepack_cpp_real rls[lenrls];
        int    ils[lenils];
    } dls1_;

    struct DLS002 {
        odepack_cpp_real stifr;
        int    newt, nsfi, nslj, njev;
        odepack_cpp_real rls2;
        int    ils2[4];
    } dls2_;
    
    struct DLSS01 {
        odepack_cpp_real con0, conmin, ccmxj, psmall, rbig, seth;
        int    iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp, ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp;
        int    ipa, lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj, nslj, ngp, nlu, nnz, nsp, nzl, nzu;
        odepack_cpp_real rlss[6];
        int ilss[34];
    } dlss_;

    struct DLSA01 {
        odepack_cpp_real tsw, cm1[12], cm2[5], pdest, pdlast, ratio, pdnorm;
        int    insufr, insufi, ixpr, icount, irflag, jtyp, mused, mxordn, mxords;
        odepack_cpp_real rlsa[22];
        int    ilsa[9];
    } dlsa_;

    struct DLSR01 {
        odepack_cpp_real t0, tlast, toutc;
        int    lg0, lg1, lgx, irfnd, itaskc, ngc, nge;
        odepack_cpp_real alpha, x2;
        int    imax, last;
        odepack_cpp_real rlsr[5];
        int    ilsr[9];
    } dlsr_;

    struct DLPK01 {
        odepack_cpp_real delt, epcon, sqrtn, rsqrtn;
        int    jpre, jacflg, locwp, lociwp, lsavx, kmp, maxl, mnewt, nni, nli, nps, ncfn, ncfl;
        odepack_cpp_real rlsp[4];
        int    ilsp[13];
    } dlpk_;
};

// opkda2.cpp
void DGEFA(odepack_cpp_real *A, int lda, int n, int *ipvt, int &info);
void DGESL(odepack_cpp_real *A, int lda, int n, int *ipvt, odepack_cpp_real *b, int job);
void DGBFA(odepack_cpp_real *abd, int lda, int n, int ml, int mu, int *ipvt, int &info);
void DGBSL(odepack_cpp_real *abd, int lda, int n, int ml, int mu, int *ipvt, odepack_cpp_real *b, int job);
void DAXPY(int n, odepack_cpp_real da, odepack_cpp_real *dx, int incx, odepack_cpp_real *dy, int incy);
void DCOPY(int n, odepack_cpp_real *dx, int incx, odepack_cpp_real *dy, int incy);
odepack_cpp_real DDOT(int n, odepack_cpp_real *dx, int incx, odepack_cpp_real *dy, int incy);
odepack_cpp_real DNRM2(int n, odepack_cpp_real *dx, int incx);
void DSCAL(int n, odepack_cpp_real da, odepack_cpp_real *dx, int incx);
int IDAMAX(int n, odepack_cpp_real *dx, int incx);
void XERRWD(std::string msg, int nmes, int nerr, int level, int ni, int i1, int i2, int nr, odepack_cpp_real r1, odepack_cpp_real r2);
int IXSAV(int ipar, int ivalue, bool iset);
int IUMACH();
}

#endif /* ODEPACK_HPP_ */