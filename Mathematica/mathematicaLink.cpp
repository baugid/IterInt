#include "integratorBase.hpp"
#include "integratorGSL.hpp"
#include <cstdint>

extern "C" {
struct cmplxType {
  double re, im;
};

cmplxType invFunc(std::complex<double> x) {
  if (x == std::complex<double>{})
    return {0, 0};
  auto res = 1. / x;
  return {res.real(), res.imag()};
}

cmplxType shift1(std::complex<double> x) {
  auto res = 1. / (x - 1.);
  return {res.real(), res.imag()};
}

cmplxType shift3(std::complex<double> x) {
  auto res = 1. / (x - 3.);
  return {res.real(), res.imag()};
}
cmplxType shift5(std::complex<double> x) {
  auto res = 1. / (x - 5.);
  return {res.real(), res.imag()};
}

cmplxType
plainIntegrate(std::complex<double> zMax,
                 std::complex<double> (**functions)(std::complex<double>),
                 std::complex<double> *residues, uint64_t count, double epsAbs,
                 double epsRel, double hStart, double hMax,
                 double startRegulator, std::complex<double> regulator) {
  iteratedIntegrals::GSLIntegrator integrator(epsAbs, epsRel, hStart, hMax,
                                              startRegulator, 2e-5,
                                              gsl_odeiv2_step_rk8pd, regulator);

  std::vector<iteratedIntegrals::Integrand<>> integrands;
  integrands.reserve(count);
  for (uint64_t i = 0; i < count; ++i) {
    integrands.emplace_back(residues[i], functions[i]);
  }

  std::complex<double> res = integrator.plainIntegrate(zMax, integrands);

  return *reinterpret_cast<cmplxType *>(&res);
}

cmplxType splittingPlainIntegrate(
    std::complex<double> zMax,
    std::complex<double> (**functions)(std::complex<double>),
    std::complex<double> *residues, uint64_t count, double epsAbs,
    double epsRel, double hStart, double hMax, double startRegulator,
    double baseRange, std::complex<double> regulator) {
  iteratedIntegrals::GSLIntegrator integrator(epsAbs, epsRel, hStart, hMax,
                                              startRegulator, baseRange,
                                              gsl_odeiv2_step_rk8pd, regulator);

  std::vector<iteratedIntegrals::Integrand<>> integrands;
  integrands.reserve(count);
  for (uint64_t i = 0; i < count; ++i) {
    integrands.emplace_back(residues[i], functions[i]);
  }

  std::complex<double> res =
      integrator.splittingPlainIntegrate(zMax, integrands);

  return *reinterpret_cast<cmplxType *>(&res);
}

void deleteCmplxArray(cmplxType *arr) { delete[] arr; }

cmplxType *
treeIntegrate(std::complex<double> zMax,
              std::complex<double> (**functions)(std::complex<double>),
              std::complex<double> *residues, uint64_t count,
              uint64_t **integrals, uint64_t integralCount, double epsAbs,
              double epsRel, double hStart, double hMax, double startRegulator,
              std::complex<double> regulator) {
  iteratedIntegrals::GSLIntegrator integrator(epsAbs, epsRel, hStart, hMax,
                                              startRegulator, 2e-5,
                                              gsl_odeiv2_step_rk8pd, regulator);

  std::vector<iteratedIntegrals::Integrand<>> integrands;
  integrands.reserve(count);
  for (uint64_t i = 0; i < count; ++i) {
    integrands.emplace_back(residues[i], functions[i]);
  }

  std::vector<std::vector<iteratedIntegrals::Integrand<> *>> integs;

  for (uint64_t i = 0; i < integralCount; ++i) {
    auto &newInt = integs.emplace_back();
    for (uint64_t j = 0; j < integrals[i][0]; ++j) {
      newInt.push_back(&integrands[integrals[i][j + 1]]);
    }
  }

  std::vector<std::complex<double>> res =
      integrator.treeIntegrate(zMax, integs);

  cmplxType *output = new cmplxType[res.size()];
  for (size_t i = 0; i < res.size(); ++i) {
    output[i] = *reinterpret_cast<cmplxType *>(&res[i]);
  }
  return output;
}

cmplxType *
splittingTreeIntegrate(std::complex<double> zMax,
                       std::complex<double> (**functions)(std::complex<double>),
                       std::complex<double> *residues, uint64_t count,
                       uint64_t **integrals, uint64_t integralCount,
                       double epsAbs, double epsRel, double hStart, double hMax,
                       double startRegulator, double baseRange,
                       std::complex<double> regulator) {
  iteratedIntegrals::GSLIntegrator integrator(epsAbs, epsRel, hStart, hMax,
                                              startRegulator, baseRange,
                                              gsl_odeiv2_step_rk8pd, regulator);

  std::vector<iteratedIntegrals::Integrand<>> integrands;
  integrands.reserve(count);
  for (uint64_t i = 0; i < count; ++i) {
    integrands.emplace_back(residues[i], functions[i]);
  }

  std::vector<std::vector<iteratedIntegrals::Integrand<> *>> integs;

  for (uint64_t i = 0; i < integralCount; ++i) {
    auto &newInt = integs.emplace_back();
    for (uint64_t j = 0; j < integrals[i][0]; ++j) {
      newInt.push_back(&integrands[integrals[i][j + 1]]);
    }
  }

  std::vector<std::complex<double>> res =
      integrator.splittingTreeIntegrate(zMax, integs);

  cmplxType *output = new cmplxType[res.size()];
  for (size_t i = 0; i < res.size(); ++i) {
    output[i] = *reinterpret_cast<cmplxType *>(&res[i]);
  }
  return output;
}
}
