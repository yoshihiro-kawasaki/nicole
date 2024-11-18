#ifndef LSODE_HPP_
#define LSODE_HPP_

#include "linear.hpp"

#include <array>
#include <cmath>
#include <memory>
#include <vector>

typedef void (*LSODA_ODE_SYSTEM_TYPE)(double t, std::vector<double> &y, std::vector<double> &dydt, void *);

typedef void (*LSODA_ODE_JACOBIAN)(double t, std::vector<double> &y, std::vector<std::vector<double> > &J, void *);

#define ETA  (2.2204460492503131e-16)

class Lsoda_cpp
    :public Linear
{
public:

    Lsoda_cpp();
    ~Lsoda_cpp();

    void Lsoda(LSODA_ODE_SYSTEM_TYPE f, const size_t neq, std::vector<double> &y,
               double *t, double tout, int itol, int itask, int *istate, int iopt, 
               std::vector<int> &iworks, std::vector<double> &rworks, 
               LSODA_ODE_JACOBIAN jac, int jt, void *_data);

    void Correction(const size_t neq, std::vector<double> &y, LSODA_ODE_SYSTEM_TYPE f,
                    size_t *corflag, double pnorm, double *del, double *delp,
                    double *told, size_t *ncf, double *rh, size_t *m,
                    LSODA_ODE_JACOBIAN jac, void *_data);

    void Stoda(const size_t neq, std::vector<double> &y, LSODA_ODE_SYSTEM_TYPE f,
               LSODA_ODE_JACOBIAN jac, void *_data);

    void Prja(const size_t neq, std::vector<double> &y, LSODA_ODE_SYSTEM_TYPE f,
              LSODA_ODE_JACOBIAN jac, void *_data);

    void TerminateIllegalInput(int *istate);
    void TerminateVariousError(std::vector<double> &y, double *t, std::vector<int> &iworks, std::vector<double> &rworks);
    void SuccessReturn(std::vector<double> &y, double *t, int itask, int ihit,
                       double tcrit, int *istate, std::vector<int> &iworks, std::vector<double> &rworks);
    void SetErrorWeightVector(const std::vector<double> &ycur);
    void Resetcoeff();
    void Solsy(std::vector<double> &y);
    void EndStoda();
    void OrderSwitch(double *rhup, double dsm, double *pdh, double *rh, size_t *orderflag);
    void InterpolatedValueOfDerivativeVector(double t, int k, std::vector<double> &dky, int *iflag);
    void Corfailure(double *told, double *rh, size_t *ncf, size_t *corflag);
    void MethodSwitch(double dsm, double pnorm, double *pdh, double *rh);
    void Cfode(int meth);
    void Scaleh(double *rh, double *pdh);
    double Fnorm(int n, const std::vector<std::vector<double>> &a, const std::vector<double> &w);
    double NormOfVectorMax(const size_t n, const std::vector<double> &v, const std::vector<double> &w);

    void StaticNumber();
    void SetTolerance(std::vector<double> atol, std::vector<double> rtol);

private:

    size_t ml_, mu_, imxer_;
    double sqrt_eta_;

    std::vector<size_t> mord_;
    std::vector<double> sm1_;

    std::vector<double> el_;
    std::vector<double> cm1_;
    std::vector<double> cm2_;


    std::vector<std::vector<double> > elco_;
    std::vector<std::vector<double> > tesco_;

    size_t illin_, init_, ierpj_, iersl_, jcur_, l_, miter_, maxord_, maxcor_, msbp_, mxncf_;

    int kflag_, jstart_;

    size_t ixpr_, jtyp_, mused_, mxordn_, mxords_;
    size_t meth_;

    size_t n_, nq_, nst_, nfe_, nje_, nqu_;
    size_t mxstp0_, mxhnl0_;
    size_t mxstep_, mxhnil_;
    size_t nslast_, nhnil_, ntrep_, nyh_;

    double ccmax_, el0_, h_;
    double hmin_, hmxi_, hu_, rc_, tn_;
    double tsw_, pdnorm_;
    double conit_, crate_, hold_, rmax_;

    size_t ialth_, ipup_, lmax_;
    size_t nslp_;
    double pdest_, pdlast_, ratio_;
    int icount_, irflag_;

    std::vector<double> ewt_;
    std::vector<double> savf_;
    std::vector<double> acor_;
    std::vector<std::vector<double>> yh_;
    std::vector<std::vector<double>> wm_;

    std::vector<int> ipvt_;

    int itol_;
    std::vector<double> rtol_;
    std::vector<double> atol_;

public:
    void *param_ = nullptr;

};




#endif /* LSODE_HPP_ */