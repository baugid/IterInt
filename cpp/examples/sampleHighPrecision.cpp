#include <iostream>
#include <iomanip>

//Headers for arbitrary precision
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>

//Load the header for BoostIntegrator
#include "../IterInt/integratorBoost.hpp"


namespace multiprec = boost::multiprecision;

int main() {
  using nr = multiprec::mpfr_float;
  using cmplx = multiprec::mpc_complex;

  //Set up the arbitrary precision types for 50 digits
  nr::default_precision(50);
  cmplx::default_precision(50);
  std::cout << std::setprecision(40);

  //Construct a BoostIntegrator using them
  iteratedIntegrals::BoostIntegrator<nr, cmplx> integrator(1e-40, 1e-30);

  //The remainder is mostly equivalent to sampleBoost.cpp and sampleGSL.cpp
  //Slight differences are there to ensure successful type deduction
  iteratedIntegrals::Integrand<cmplx> cot{1, [](const cmplx& x){if(x == 0.) return cmplx{}; else return cmplx{1}/tan(x);}};
  iteratedIntegrals::Integrand<cmplx> sq{0, [](const cmplx& x){return x * x;}};
  iteratedIntegrals::Integrand<cmplx> inv{1, [](const cmplx& x){if(x == 0.) return cmplx{}; else return cmplx{1}/x;}};

  cmplx res = integrator.plainIntegrate(
    cmplx{3/2., 0},
    {cot, sq, inv}
  );
  std::cout << res << std::endl;

  auto resList = integrator.treeIntegrate(
    cmplx{3/2., 0},{{&cot,&sq,&inv},{&cot,&inv,&sq}});
  std::cout << resList[0]<< "\t" << resList[1] << std::endl;

  resList = integrator.treeIntegrate(
      cmplx{3/2., 0},
      {{&cot, &sq, &inv}, {&cot, &inv, &sq}},
      [&inv](auto integ){
        if(integ[0].isRegulated() && integ[0].fct == &inv)
          return std::optional(cmplx{});
        return std::optional<cmplx>{};
      });
  std::cout << resList[0]<< "\t" << resList[1] << std::endl;
}
