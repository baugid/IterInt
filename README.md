# IterInt

A small library for C++ or Mathematica to quickly numerically evaluate iterated Integrals.
Detailed information can be found in the accompaning paper: [2606.02744](http://arxiv.org/abs/2606.02744)

## Getting Started

To obtain the code please just clone/download this repository.
Both the C++ and Mathematica versions provide a similar feature set and are generally independent from each other.
We also provide an interface to call the C++ code from within Mathematica.
In this repository both versions as well as accompaning example codes are provided.
For a detailed description of the available feature set we refer again to the [paper](http://arxiv.org/abs/2606.02744).

### Mathematica

To use the library no further dependencies are necessary. The main parts of the library should work with Mathematica 10.0 or newer.
It is validated to work with Mathematica 12. In case the interface to the C++ code is to be used, Mathematica 13.3 or newer is required.

The example notebook can be easily evaluated, if `IterInt.m` is located in the same directory.
In general the library can directly be loaded using
```Mathematica
Get["IterInt`"]
```
if it can be found via the current `$PATH`.

### C++
The C++ code depends on [GSL](https://www.gnu.org/software/gsl/) or [Boost](https://www.boost.org/).
For some examples also [GNU MPC](https://www.multiprecision.org/) and [GNU MPFR](https://www.mpfr.org/) are required.
Furthermore, features from C+++20 are used.

In the following we describe how to compile all provided source code.
We assume that these commands are invoked in the respective directory and `g++` is installed and working.

The `sampleGSL.cpp` file needs to be linked against GSL:
```sh
g++ --std=c++20 -lgsl sampleGSL.cpp
```

For `sampleBoost.cpp` no linker flags are necessary:
```sh
g++ --std=c++20 sampleBoost.cpp
```

The arbitrary precision example `sampleHighPrecision.cpp` requires linking against the corresponding libraries:
```sh
g++ --std=c++20 -lmpfr -lmpc sampleHighPrecision.cpp
```

The library to use the C++ code from within Mathematica is provided both as source code and as a compiled version.
To compile this yourself, if necessary due to a different system architecture or operating system, run
```sh
g++ --std=c++20 -lgsl -fPIC -shared -o libIntegrate.so mathematicaLink.cpp
```
in the respective directory.

## License

This library is made available under the GPL-3.0 license. For further information see the LICENSE file.
