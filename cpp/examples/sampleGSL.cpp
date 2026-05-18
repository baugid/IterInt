#include <iostream>
#include <iomanip>

//Load the header for GSLIntegrator
#include "../IterInt/integratorGSL.hpp"


int main() {
  using cmplx = std::complex<double>;
  std::cout << std::setprecision(16);

  //Construct a GSLIntegrator with 12 digits of absolute precision
  iteratedIntegrals::GSLIntegrator integrator(1e-12, 0);

  //Define three integration kernels
  iteratedIntegrals::Integrand<> cot{1, [](const cmplx& x){if(x == 0.) return cmplx{}; else return 1./tan(x);}};
  iteratedIntegrals::Integrand<> sq{[](const cmplx& x){return x * x;}};
  iteratedIntegrals::Integrand<> inv{1, [](const cmplx& x){if(x == 0.) return cmplx{}; else return 1./x;}};

  //Evaluate the first integral
  cmplx res = integrator.plainIntegrate(
    cmplx{3/2., 0},
    {cot, sq, inv}
  );
  std::cout << res << std::endl;

  //Now we compute two integrals simultaneously. Here the integrals are given by lists of pointers
  auto resList = integrator.treeIntegrate(
    cmplx{3/2., 0},
    {{&cot, &sq, &inv}, {&cot, &inv, &sq}}
  );
  std::cout << resList[0]<< "\t" << resList[1] << std::endl;

  //An optional parameter can be given, which skips the calculation of certain integrals, when they are trivial.
  //The type hidden behind auto is std::span<iteratedIntegrals::TaggedFunction>
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
