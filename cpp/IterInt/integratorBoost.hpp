#ifndef ITERATED_BOOST_INTEGRATOR_H
#define ITERATED_BOOST_INTEGRATOR_H
#include "integratorBase.hpp"
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/stepper/generation.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>

namespace iteratedIntegrals {
namespace {
namespace odeint = boost::numeric::odeint;
}

/**
 * Integrator using Boost as a backend.
 * @see IntegratorBase
 * @tparam stepper The stepper to use for integration
 **/
template <typename nr = double, typename cmplx = std::complex<nr>,
          typename stepper =
              odeint::runge_kutta_fehlberg78<std::vector<cmplx>, nr>>
class BoostIntegrator : public IntegratorBase<nr, cmplx> {
public:
  using StateType = std::vector<cmplx>;
  using FunctionType = IntegratorBase<nr, cmplx>::FunctionType;

public:
  /**
   * Constructor involving all configuration options. Everything has provided
   * default values.
   * @param epsAbs absolute deviation targeted by the integration
   * @param epsRel relative deviation targeted by the integration
   * @param hStart initial step size
   * @param hMax maximum step size. All values above 1 are equivalent
   * @param startRegulator initial offset to avoid integrable divergences
   * @param switchPoint Cross over point for splitting algorithms
   * @param regulator tangent vector for shuffle regularization
   */
  BoostIntegrator(std::type_identity_t<nr> epsAbs = 1e-12,
                  std::type_identity_t<nr> epsRel = 1e-12,
                  std::type_identity_t<nr> hStart = 1e-5,
                  std::type_identity_t<nr> hMax = 1,
                  std::type_identity_t<nr> startRegulator = 0,
                  std::type_identity_t<nr> switchPoint = 2e-5,
                  const std::type_identity_t<cmplx> &regulator = cmplx{1, 0})
      : IntegratorBase<nr, cmplx>(startRegulator, switchPoint, regulator),epsAbs(epsAbs),epsRel(epsRel), hStart(hStart), hMax(hMax) {}

private:
  /*
   * structs for the different kinds of differential equation
   */

  struct System {
    void operator()(const StateType &y, StateType &dy, const nr t) {
      cmplx zT = t * zMax;
      dy[0] = zMax * eval.front()(zT);

      for (size_t i = 1; i < eval.size(); ++i) {
        dy[i] = eval[i](zT) * zMax * y[i - 1];
      }
    }

    const std::vector<FunctionType> &eval;
    const cmplx &zMax;
  };

  struct TreeSystem {
    void operator()(const StateType &y, StateType &dy, const nr t) {
      cmplx zT = t * zMax;
      dy[0] = zMax * tree.front().second->operator()(zT);

      for (size_t i = 1; i < tree.size(); ++i) {
        dy[i] = tree[i].second->operator()(zT) * zMax * y[tree[i].first];
      }
    }
    const std::vector<std::pair<size_t, FunctionType *>> &tree;
    const cmplx &zMax;
  };

protected:
  /*
   * Implementations of the virtual members.
   */

  void
  callTreeIntegrator(const cmplx &zMax,
                     const std::vector<std::pair<size_t, FunctionType *>> &tree,
                     StateType &state, const nr &start) override {
    auto basicStepper = stepper();
    odeint::integrate_adaptive(
        odeint::make_controlled(epsAbs, epsRel, hMax, basicStepper),
        TreeSystem{tree, zMax}, state, start, nr{1}, hStart);
  }

  void callBasicIntegrator(const cmplx &zMax,
                           const std::vector<FunctionType> &integral,
                           StateType &state, const nr &start) override {
    auto basicStepper = stepper();
    odeint::integrate_adaptive(
        odeint::make_controlled(epsAbs, epsRel, hMax, basicStepper),
        System{integral, zMax}, state, start, nr{1}, hStart);
  }

protected:
  nr epsAbs;
  nr epsRel;
  nr hStart;
  nr hMax;
};
} // namespace iteratedIntegrals
#endif
