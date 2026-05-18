#ifndef ITERATED_GSL_INTEGRATOR_H
#define ITERATED_GSL_INTEGRATOR_H
#include "integratorBase.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>

namespace iteratedIntegrals {
/**
 * Integrator using GSL as a backend.
 * @see IntegratorBase
 **/
class GSLIntegrator : public IntegratorBase<> {
public:
  using cmplx = std::complex<double>;
  using StateType = IntegratorBase<>::StateType;
  using FunctionType = IntegratorBase<>::FunctionType;

public:
  /**
   * Constructor passing on all the parameters of IntegratorBase.
   * In addition the integration method can be specified.
   **/
  /**
   * Constructor involving all configuration options. Everything has provided
   * default values.
   * @param epsAbs absolute deviation targeted by the integration
   * @param epsRel relative deviation targeted by the integration
   * @param hStart initial step size
   * @param hMax maximum step size. All values above 1 are equivalent
   * @param startRegulator initial offset to avoid integrable divergences
   * @param switchPoint Cross over point for splitting algorithms
   * @param gslAlgorithm the algorithm used by GSL for integration
   * @param regulator tangent vector for shuffle regularization
   */
  GSLIntegrator(
      double epsAbs = 1e-12, double epsRel = 1e-12, double hStart = 1e-5,
      double hMax = 1, double startRegulator = 0, double switchPoint = 2e-5,
      const gsl_odeiv2_step_type *gslAlgorithm = gsl_odeiv2_step_rk8pd,
      const cmplx &regulator = cmplx{1, 0})
      : IntegratorBase(startRegulator, switchPoint, regulator), epsAbs(epsAbs),
        epsRel(epsRel), hStart(hStart), hMax(hMax), gslAlgorithm(gslAlgorithm) {
  }

private:
  /*
   * Data structures and functions for the differential equations
   */
  struct DglData {
    const std::vector<FunctionType> &eval;
    cmplx zMax;
  };
  struct TreeDglData {
    const std::vector<std::pair<size_t, FunctionType *>> &tree;
    cmplx zMax;
  };

  static int calcDerivatives(double t, const double *y, double *dy,
                             void *params) {
    DglData *data = reinterpret_cast<DglData *>(params);
    const std::vector<FunctionType> &eval = data->eval;
    cmplx zMax = data->zMax;
    cmplx zT = t * zMax;

    cmplx cDer = zMax * eval.front()(zT);
    dy[0] = cDer.real();
    dy[1] = cDer.imag();

    for (size_t i = 1; i < eval.size(); ++i) {
      cmplx cDer =
          eval[i](zT) * zMax * cmplx(y[2 * (i - 1)], y[2 * (i - 1) + 1]);
      dy[2 * i] = cDer.real();
      dy[2 * i + 1] = cDer.imag();
    }
    return GSL_SUCCESS;
  }

  static int calcTreeDerivatives(double t, const double *y, double *dy,
                                 void *params) {
    TreeDglData *data = reinterpret_cast<TreeDglData *>(params);
    const std::vector<std::pair<size_t, FunctionType *>> &tree = data->tree;
    cmplx zMax = data->zMax;
    cmplx zT = t * zMax;

    cmplx cDer = zMax * tree.front().second->operator()(zT);
    dy[0] = cDer.real();
    dy[1] = cDer.imag();

    for (size_t i = 1; i < tree.size(); ++i) {
      cmplx cDer = tree[i].second->operator()(zT) * zMax *
                   cmplx(y[2 * tree[i].first], y[2 * tree[i].first + 1]);
      dy[2 * i] = cDer.real();
      dy[2 * i + 1] = cDer.imag();
    }

    return GSL_SUCCESS;
  }

protected:
  /*
   * Overloads of the virtual members
   */

  void
  callTreeIntegrator(const cmplx &zMax,
                     const std::vector<std::pair<size_t, FunctionType *>> &tree,
                     StateType &state, const double &start) override {
    TreeDglData dat{tree, zMax};
    gsl_odeiv2_system sys = {calcTreeDerivatives, nullptr, 2 * tree.size(),
                             &dat};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_standard_new(
        &sys, gslAlgorithm, hStart, epsAbs, epsRel, 1, 1);
    gsl_odeiv2_driver_set_hmax(d, hMax);

    double t = start;
    // the standard requires the reinterpret_cast to work
    int status = gsl_odeiv2_driver_apply(
        d, &t, 1, reinterpret_cast<double *>(state.data()));

    if (status != GSL_SUCCESS)
      std::cerr << "GSL error, return code = " << status << " reached t=" << t
                << "\n";

    gsl_odeiv2_driver_free(d);
  }

  void callBasicIntegrator(const cmplx &zMax,
                           const std::vector<FunctionType> &integral,
                           StateType &state, const double &start) override {
    DglData dat{integral, zMax};
    gsl_odeiv2_system sys = {calcDerivatives, nullptr, 2 * integral.size(),
                             &dat};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_standard_new(
        &sys, gslAlgorithm, hStart, epsAbs, epsRel, 1, 1);
    gsl_odeiv2_driver_set_hmax(d, hMax);

    double t = start;
    // the standard requires the reinterpret_cast to work
    int status = gsl_odeiv2_driver_apply(
        d, &t, 1, reinterpret_cast<double *>(state.data()));

    if (status != GSL_SUCCESS)
      std::cerr << "GSL error, return code = " << status << " reached t=" << t
                << "\n";

    gsl_odeiv2_driver_free(d);
  }

protected:
  double epsAbs;
  double epsRel;
  double hStart;
  double hMax;
  const gsl_odeiv2_step_type *gslAlgorithm;
};
} // namespace iteratedIntegrals

#endif
