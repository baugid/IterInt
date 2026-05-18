#include <iostream>
#include <iomanip>
#include <cmath>

//Load the header for BoostIntegrator
#include "../IterInt/integratorBoost.hpp"


int main() {
  using cmplx = std::complex<double>;
  std::cout << std::setprecision(16);

  //Construct a BoostIntegrator
  iteratedIntegrals::BoostIntegrator integrator(1e-12, 0);

  //The remainder of this file is identical to sampleGSL.cpp

  iteratedIntegrals::Integrand<> cot{1, [](const cmplx& x){if(x == 0.) return cmplx{}; else return 1./tan(x);}};
  iteratedIntegrals::Integrand<> sq{0, [](const cmplx& x){return x * x;}};
  iteratedIntegrals::Integrand<> inv{1, [](const cmplx& x){if(x == 0.) return cmplx{}; else return 1./x;}};

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
