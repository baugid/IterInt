#pragma once
#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <ranges>
#include <span>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace iteratedIntegrals {
// private namespace for internal helpers
namespace internals {
// simple hash function for containers of objects of type T with corresponding
// hash function memberHash
template <typename T, typename memberHash = std::hash<std::remove_cv_t<T>>>
struct ContainerHash {
  template <typename U> size_t operator()(const U &vec) const {
    size_t seed = vec.size();
    for (const T &i : vec) {
      seed ^= memberHash{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

// deep equals for spans. This is necessary for std::unordered_set
template <typename T> struct SpanEquals {
  bool operator()(const std::span<T> &a, const std::span<T> &b) const {
    if (a.size() != b.size())
      return false;
    for (size_t i = 0; i < a.size(); ++i)
      if (a[i] != b[i])
        return false;
    return true;
  }
};

// implementation of the shuffle product for short lengths using bit magic
template <std::forward_iterator iter1, std::forward_iterator iter2, typename F>
void fastShuffleProdImpl(iter1 aBegin, size_t aLen, iter2 bBegin, size_t bLen,
                         F callback) {
  using elemType = std::remove_cv_t<std::remove_reference_t<decltype(*bBegin)>>;
  using ull = unsigned long long;

  size_t totalLength = aLen + bLen;
  std::vector<elemType> entry(totalLength);

  ull stop = 1 << totalLength;

  ull current = (1 << aLen) - 1; // one bits at the back

  while (current < stop) {
    iter1 aRunning = aBegin;
    iter2 bRunning = bBegin;

    ull shifted = current;
    for (size_t i = 0; i < totalLength; i++) {
      if (shifted & 1)
        entry[i] = *(aRunning++);
      else
        entry[i] = *(bRunning++);
      shifted >>= 1;
    }
    callback(entry);

    // next permutation see
    // https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation

    ull t =
        current | (current - 1); // t gets v's least significant 0 bits set to 1
    // Next set to 1 the most significant bit to change,
    // set to 0 the least significant ones, and add the necessary 1 bits.
    current = (t + 1) | (((~t & -~t) - 1) >> (std::countr_zero(current) + 1));
  }
}

// implementation of the shuffle product for arbitrary lengths
template <std::forward_iterator iter1, std::forward_iterator iter2, typename F>
void slowShuffleProdImpl(iter1 aBegin, size_t aLen, iter2 bBegin, size_t bLen,
                         F callback) {
  using elemType = std::remove_cv_t<std::remove_reference_t<decltype(*bBegin)>>;
  std::vector<uint8_t> currentPerm(aLen + bLen, 0);
  std::fill(currentPerm.begin() + aLen, currentPerm.end(), 1);

  std::vector<elemType> entry(currentPerm.size());
  do {
    iter1 aRunning = aBegin;
    iter2 bRunning = bBegin;

    for (size_t i = 0; i < currentPerm.size(); i++) {
      if (currentPerm[i] == 0)
        entry[i] = *(aRunning++);
      else
        entry[i] = *(bRunning++);
    }
    callback(entry);
  } while (std::next_permutation(currentPerm.begin(), currentPerm.end()));
}

/**
 * shuffle product implementation, which calls a callback for each term.
 **/
template <std::forward_iterator iter1, std::forward_iterator iter2,
          std::forward_iterator iter3, std::forward_iterator iter4, typename F>
void shuffleProduct(iter1 aBegin, iter2 aEnd, iter3 bBegin, iter4 bEnd,
                    F callback) {
  using elemType = std::remove_cv_t<std::remove_reference_t<decltype(*bBegin)>>;
  if (aBegin == aEnd && bBegin == bEnd) {
    std::vector<elemType> out;
    callback(out);
    return;
  }

  if (aBegin == aEnd) {
    std::vector<elemType> out;
    out.insert(out.end(), bBegin, bEnd);
    callback(out);
    return;
  }

  if (bBegin == bEnd) {
    std::vector<elemType> out;
    out.insert(out.end(), aBegin, aEnd);
    callback(out);
    return;
  }

  size_t aLen = std::distance(aBegin, aEnd);
  size_t bLen = std::distance(bBegin, bEnd);
  if (aLen + bLen < 8 * sizeof(unsigned long long)) {
    fastShuffleProdImpl(aBegin, aLen, bBegin, bLen, callback);
  } else {
    slowShuffleProdImpl(aBegin, aLen, bBegin, bLen, callback);
  }
}
} // namespace internals

/**
 * @brief class modelling an integration kernel, which consists out of a
 *function and its residue in zero.
 * @tparam cmplx the data type for complex numbers to use. It should work
 *similar to std::complex
 **/
template <typename cmplx = std::complex<double>> class Integrand {
public:
  using FunctionType = std::function<cmplx(const cmplx &)>;

public:
  /**
   * Default constructor for an empty object.
   **/
  Integrand() {}

  /**
   * Constructor that takes both the residue and the function as a parameter.
   * The latter is given as a std::function object
   **/
  Integrand(cmplx residue, const FunctionType &function)
      : residue(residue), fun(function) {}

  /**
   * Constructor that takes only the function as a parameter. The residue is
   * defaulted to zero.
   **/
  Integrand(const FunctionType &function) : Integrand(cmplx{}, function) {}

  /**
   * Constructor overload for arbitrary callable types. The parameter is
   * forwarded to the constructor of std::function.
   **/
  template <typename F>
    requires std::constructible_from<FunctionType, F>
  Integrand(F func)
      : Integrand(static_cast<FunctionType>(std::function{func})) {}

  /**
   * Constructor overload for arbitrary callable types. The second parameter is
   * forwarded to the constructor of std::function.
   **/
  template <typename F>
  Integrand(cmplx residue, F func)
    requires std::constructible_from<FunctionType, F>
      : Integrand(residue, static_cast<FunctionType>(std::function{func})) {}

  /**
   * Returns a function that corresponds to the stored function minus its pole.
   * Using the returned object after this object is destroyed leads to undefined
   * behaviour.
   **/
  FunctionType regulated() const {
    if (residue == cmplx{})
      return fun;

    return [*this](const cmplx &x) -> cmplx {
      if (x == cmplx{})
        return cmplx{};
      return fun(x) - residue / x;
    };
  }

  /**
   * Returns a function representing only the pole term.
   **/
  FunctionType poleTerm() const {
    return [res = residue](const cmplx &x) -> cmplx {
      if (x == cmplx{})
        return cmplx{};

      return res / x;
    };
  }

  /**
   * Decays this object into its stored std::function.
   **/
  operator FunctionType() const { return fun; }

  /**
   * Returns true, iff the residue of the stored function is zero.
   **/
  bool isFinite() const { return residue == cmplx{}; }

  /**
   * Returns the residue of this function.
   **/
  cmplx getResidue() const { return residue; }

private:
  cmplx residue;
  FunctionType fun;
};

/**
 * Struct representing a reference to a function through its pointer
 * and some flags storing additional information.
 **/
template <typename cmplx = std::complex<double>> struct TaggedFunction {
  const Integrand<cmplx> *fct;
  unsigned char flags;

  /**
   * enum representing all present flags. A function can be regulated, the pole
   * term, or finite. Generally every combination, even redundant ones, are
   * possible
   */
  enum : unsigned char { REGULATED = 1, POLE = 2, FINITE = 4 };

  /*
   * Helper functions to check the flags
   */

  bool isRegulated() const { return flags & REGULATED; }
  bool isPole() const { return flags & POLE; }
  bool isNormal() const { return (flags & (REGULATED | POLE)) == 0; }
  bool isFinite() const { return flags & FINITE; }

  TaggedFunction toRegulated() const {
    return TaggedFunction(fct, flags | REGULATED);
  }

  TaggedFunction toPole() const { return TaggedFunction(fct, flags | POLE); }

  /**
   * Obtain a callable object derived from this TaggedFunction
   **/
  Integrand<cmplx>::FunctionType toCallable() const {
    if (isNormal())
      return *fct;
    else if (isRegulated())
      return fct->regulated();
    else // must be pole
      return fct->poleTerm();
  }

  bool operator==(const TaggedFunction &) const = default;
};
} // namespace iteratedIntegrals

/**
 * std::hash specialization for iteratedIntegrals::TaggedFunction objects.
 **/
template <typename cmplx>
struct std::hash<iteratedIntegrals::TaggedFunction<cmplx>> {
  size_t operator()(const iteratedIntegrals::TaggedFunction<cmplx> &f) const {
    size_t hash1 =
        std::hash<const iteratedIntegrals::Integrand<cmplx> *>{}(f.fct);
    size_t hash2 = std::hash<unsigned char>{}(f.flags);
    return hash1 ^ (hash2 + 0x9e3779b9 + (hash1 << 6) + (hash1 >> 2));
  }
};

namespace iteratedIntegrals {

/**
 * Class that performs the regularization of iterated integrals.
 * @tparam nr Floating point data type
 * @tparam cmplx data type for complex numbers, which must be compatible with nr
 */
template <typename nr = double, typename cmplx = std::complex<nr>>
class Regulator {
public:
  using FunctionType = Integrand<cmplx>::FunctionType;
  /**
   * Data type for a single regulated integral.
   * The outer vector represents the sum over multiple integrals, represented by
   * the second entry of the pair, with a prefactor, which is given by the first
   * entry.
   **/
  using regulatedIntegral =
      std::vector<std::pair<cmplx, std::vector<TaggedFunction<cmplx>>>>;

public:
  /**
   * constructs the Regulator from the tangent vector in zero.
   **/
  Regulator(const std::type_identity_t<cmplx> &tangent) : tangent(tangent) {}

  /**
   * regulates the iterated integral given by an upper integration bound and a
   * vector of integration kernels using the residues provided.
   * @param zMax upper integration bound
   * @param taggedIntegr integration kernels as TaggedFunction objects
   **/
  std::pair<cmplx, regulatedIntegral>
  regulate(const cmplx &zMax,
           const std::vector<TaggedFunction<cmplx>> &taggedIntegr) {
    cmplx factor = cmplx{1, 0};
    using std::log;
    cmplx logV = log(zMax / tangent);

    regulatedIntegral res;

    for (size_t i = 0; i < taggedIntegr.size(); i++) {
      if (factor == cmplx{})
        break;

      rFunction(taggedIntegr, i, factor, res);

      factor *=
          taggedIntegr[i].fct->getResidue() * logV / static_cast<nr>(i + 1);
    }

    return {factor, res};
  }

private:
  /**
   * implements the R operator for regularization.
   * @param offset denotes the number of kernels being skipped.
   * @param factor the prefactor infront of R
   * @param res output parameter
   **/
  void rFunction(const std::vector<TaggedFunction<cmplx>> &integrals,
                 size_t offset, const cmplx factor, regulatedIntegral &res) {
    if (offset >= integrals.size())
      return;

    std::vector<TaggedFunction<cmplx>> integral;
    integral.resize(integrals.size() - offset);

    cmplx prefac = -factor;
    for (size_t i = offset; i < integrals.size(); i++) {
      prefac = -prefac;

      auto poles = integrals | std::views::take(i) | std::views::drop(offset) |
                   std::views::reverse |
                   std::views::transform(&TaggedFunction<cmplx>::toPole);

      integral[0] = integrals[i].toRegulated();

      internals::shuffleProduct(
          integrals.begin() + i + 1, integrals.end(), poles.begin(),
          poles.end(),
          [&integral, &res, prefac](std::vector<TaggedFunction<cmplx>> &shuf) {
            std::move(shuf.begin(), shuf.end(), integral.begin() + 1);
            res.emplace_back(prefac, integral);
          });

      if (integrals[i].isFinite())
        break;
    }
  }

private:
  cmplx tangent;
};

namespace internals {
/**
 * Simple wrapper caching the last value of its member function.
 **/
template <typename cmplx> struct CachingFunction {
  cmplx operator()(const cmplx &z) {
    if (z == lastZ)
      return lastVal;

    lastZ = z;
    return lastVal = fun(z);
  }

  std::function<cmplx(const cmplx &)> fun;
  cmplx lastZ{-1};
  cmplx lastVal{};
};
} // namespace internals

/**
 * Abstract base class for the integrators.
 * @tparam nr floating point type
 * @tparam cmplx complex number type
 **/
template <typename nr = double, typename cmplx = std::complex<nr>>
class IntegratorBase {
public:
  using StateType = std::vector<cmplx>;
  using FunctionType = std::function<cmplx(const cmplx &)>;

public:
  /**
   * Constructor involving configuration options. Everywhere sensible defaults are provided
   * @param startRegulator initial offset to avoid integrable divergences
   * @param switchPoint Cross over point for splitting algorithms
   * @param regulator tangent vector for shuffle regularization
   */
  IntegratorBase(std::type_identity_t<nr> startRegulator = 0,
                 std::type_identity_t<nr> switchPoint = 2e-5,
                 const std::type_identity_t<cmplx> &regulator = cmplx{1, 0})
      : switchPoint(switchPoint), reg(regulator), startRegulator(startRegulator) {}

  /*
   * Overloads to make the check function optional.
   */

  cmplx plainIntegrate(cmplx zMax,
                         const std::vector<Integrand<cmplx>> &integral) {
    return plainIntegrate(
        zMax, integral,
        [](const std::span<TaggedFunction<cmplx>> &) -> std::optional<cmplx> {
          return std::optional<cmplx>();
        });
  }

  cmplx
  splittingPlainIntegrate(cmplx zMax,
                            const std::vector<Integrand<cmplx>> &integral) {
    return splittingPlainIntegrate(
        zMax, integral,
        [](const std::span<TaggedFunction<cmplx>> &) -> std::optional<cmplx> {
          return std::optional<cmplx>();
        });
  }

  std::vector<cmplx>
  treeIntegrate(cmplx zMax,
                const std::vector<std::vector<Integrand<cmplx> *>> &integrals) {
    return treeIntegrate(
        zMax, integrals,
        [](std::span<TaggedFunction<cmplx>>) -> std::optional<cmplx> {
          return std::optional<cmplx>();
        });
  }

  std::vector<cmplx> splittingTreeIntegrate(
      cmplx zMax,
      const std::vector<std::vector<Integrand<cmplx> *>> &integrals) {
    return splittingTreeIntegrate(
        zMax, integrals,
        [](std::span<TaggedFunction<cmplx>>) -> std::optional<cmplx> {
          return std::optional<cmplx>();
        });
  }

  /**
   * evaluates a single integral using the most basic approach.
   * @param zMax the upper integration bound
   * @param integral the integration kernels in order
   * @param checkIntegral a function that returns an optional when given a span
   * of TaggedFunctions. Its value, if present, replaces the passed in integral
   **/
  cmplx plainIntegrate(cmplx zMax,
                         const std::vector<Integrand<cmplx>> &integral,
                         auto checkIntegral)
    requires std::is_convertible_v<
        decltype(checkIntegral(std::span<TaggedFunction<cmplx>>{})),
        std::optional<cmplx>>
  {
    std::vector<FunctionType> integrExpr;
    std::vector<TaggedFunction<cmplx>> taggedIntegr;
    integrExpr.reserve(integral.size());
    taggedIntegr.reserve(integral.size());

    for (size_t i = 0; i < integral.size(); i++) {
      taggedIntegr.emplace_back(&integral[i], TaggedFunction<cmplx>::FINITE *
                                                  integral[i].isFinite());
    }

    // just integrate everything completely
    auto [res, regulated] = reg.regulate(zMax, taggedIntegr);

    StateType state;
    for (auto &[fac, integ] : regulated) {
      std::optional<cmplx> checkRes = checkIntegral(std::span{integ});
      if (checkRes) {
        res += fac * (*checkRes);
        continue;
      }

      state.resize(integ.size());

      for (size_t i = 0; i < integ.size(); i++) {
        integrExpr.push_back(integ[i].toCallable());
      }

      res += fac * evalIntegral(zMax, integrExpr, state, startRegulator);
      integrExpr.clear();
      state.clear();
    }

    return res;
  }

  /**
   * evaluates a single integral by splitting the integration path and
   * performing the regulation only for the first section of them.
   * @param zMax the upper integration bound
   * @param integral the integration kernels in order
   * @param checkIntegral a function that returns an optional when given a span
   * of TaggedFunctions. Its value, if present, replaces the passed in integral
   **/
  cmplx splittingPlainIntegrate(cmplx zMax,
                                  const std::vector<Integrand<cmplx>> &integral,
                                  auto checkIntegral)
    requires std::is_convertible_v<
        decltype(checkIntegral(std::span<TaggedFunction<cmplx>>{})),
        std::optional<cmplx>>
  {
    // init data on functions
    std::vector<FunctionType> integrExpr;
    std::vector<TaggedFunction<cmplx>> taggedIntegr;

    integrExpr.reserve(integral.size());
    taggedIntegr.reserve(integral.size());

    for (size_t i = 0; i < integral.size(); i++) {
      taggedIntegr.emplace_back(&integral[i], TaggedFunction<cmplx>::FINITE *
                                                  integral[i].isFinite());
    }

    // vector for all regulated partial integrals
    std::vector<decltype(reg.regulate(zMax, taggedIntegr))> regulatedExpr;
    regulatedExpr.reserve(integral.size());

    // list of all integrals to compute
    std::vector<std::span<TaggedFunction<cmplx>>> toCompute;

    // regulate everything
    {
      // Helper for regularization calculation
      std::vector<TaggedFunction<cmplx>> toRegulate;
      toRegulate.reserve(integral.size());

      for (size_t i = 0; i < integral.size(); i++) {
        toRegulate.push_back(taggedIntegr[i]);
        // regulate for switch over point
        auto regRes = reg.regulate(zMax * switchPoint, toRegulate);
        regulatedExpr.emplace_back(std::move(regRes));

        for (auto &integral : regulatedExpr.back().second)
          toCompute.emplace_back(integral.second);
      }
    }

    // longest integrals have to be computed first. This ensures nothing is
    // done redundantly
    std::sort(toCompute.begin(), toCompute.end(),
              [](const auto &first, const auto &second) {
                return second.size() < first.size();
              });

    // map to prevent recompution of any iterated integral
    std::unordered_map<std::span<TaggedFunction<cmplx>>, cmplx,
                       internals::ContainerHash<TaggedFunction<cmplx>>,
                       internals::SpanEquals<TaggedFunction<cmplx>>>
        valueMap;
    valueMap.reserve(toCompute.size());

    for (const auto &integ : toCompute) {
      auto mapIter = valueMap.try_emplace(integ, cmplx{});
      if (!mapIter.second) // already calculated
        continue;

      std::optional<cmplx> integCheck = checkIntegral(integ);
      if (integCheck) { // user calculated
        mapIter.first->second = *integCheck;
      } else { // we need to work
        // construct integrand for baseEval
        for (size_t i = 0; i < integ.size(); i++) {
          integrExpr.push_back(integ[i].toCallable());
        }

        auto resVec = baseEval(zMax * switchPoint, integrExpr);

        // fill results into map
        mapIter.first->second = resVec.back();
        for (size_t i = 1; i < integ.size(); ++i) {
          valueMap.try_emplace(integ.first(i), resVec[i - 1]);
        }

        // we clear the integral for next use without reallocation
        integrExpr.clear();
      }
    }

    // regulated initial values
    std::vector<cmplx> regulatedVals(integral.size());

    for (size_t i = 0; i < regulatedExpr.size(); i++) {
      cmplx res = regulatedExpr[i].first;
      auto &regIntegrals = regulatedExpr[i].second;
      for (auto &[fac, integ] : regIntegrals) {
        res += fac * valueMap[integ];
      }

      regulatedVals[i] = std::move(res);
    }

    return evalIntegral(zMax, {integral.begin(), integral.end()}, regulatedVals,
                        switchPoint);
  }

  /**
   * computes a list of iterated integrals at once by combining them into trees.
   * @param zMax the upper integration bound
   * @param integrals a vector of the integrals represented by a vector of
   * pointers to integrands
   * @param checkIntegral a function that returns an optional when given a span
   * of TaggedFunctions. Its value, if present, replaces the passed in integral
   **/
  std::vector<cmplx>
  treeIntegrate(cmplx zMax,
                const std::vector<std::vector<Integrand<cmplx> *>> &integrals,
                auto checkIntegral)
    requires std::is_convertible_v<
        decltype(checkIntegral(std::span<TaggedFunction<cmplx>>{})),
        std::optional<cmplx>>
  {
    using spanSet =
        std::unordered_set<std::span<TaggedFunction<cmplx>>,
                           internals::ContainerHash<TaggedFunction<cmplx>>,
                           internals::SpanEquals<TaggedFunction<cmplx>>>;
    // map the integrands to TaggedFunctions
    std::unordered_map<TaggedFunction<cmplx>, FunctionType> fctMap;

    for (size_t i = 0; i < integrals.size(); i++) {
      for (auto integrand : integrals[i]) {
        fctMap.try_emplace(
            TaggedFunction<cmplx>(integrand, TaggedFunction<cmplx>::FINITE *
                                                 integrand->isFinite()),
            internals::CachingFunction<cmplx>(*integrand));
      }
    }

    // regulate everything and keep track of the appearing integrals
    // Here also all integrals are owned
    std::vector<decltype(reg.regulate(zMax, {}))> regulatedIntegrals;
    regulatedIntegrals.reserve(integrals.size());

    std::unordered_map<std::span<TaggedFunction<cmplx>>, cmplx,
                       internals::ContainerHash<TaggedFunction<cmplx>>,
                       internals::SpanEquals<TaggedFunction<cmplx>>>
        integralValueMap;

    std::vector<std::vector<std::pair<size_t, FunctionType *>>> groups;
    std::unordered_map<TaggedFunction<cmplx>, size_t> headToIndex;
    std::unordered_map<std::span<TaggedFunction<cmplx>>, size_t,
                       internals::ContainerHash<TaggedFunction<cmplx>>,
                       internals::SpanEquals<TaggedFunction<cmplx>>>
        integToIndex;

    {
      spanSet allHeadIntegrals;
      for (std::vector<TaggedFunction<cmplx>> toRegulate;
           const auto &integral : integrals) {
        toRegulate.clear();
        for (size_t i = 0; i < integral.size(); i++)
          toRegulate.emplace_back(integral[i], TaggedFunction<cmplx>::FINITE *
                                                   integral[i]->isFinite());

        regulatedIntegrals.push_back(reg.regulate(zMax, toRegulate));

        // add all integrals to the relevant sets
        for (auto &regPair : regulatedIntegrals.back().second) {
          std::span<TaggedFunction<cmplx>> integ{regPair.second.data(),
                                                 regPair.second.size()};
          if (allHeadIntegrals.contains(integ)) // already handled
            continue;
          allHeadIntegrals.emplace(integ);

          std::optional<cmplx> checkRes = checkIntegral(integ);
          if (checkRes) { // value externally provided
            integralValueMap[{integ.data(), integ.size()}] = *checkRes;
            continue;
          }

          // create and add to trees if necessary
          auto headIter = headToIndex.try_emplace(integ.front(), groups.size());
          if (headIter.second) { // new head
            auto fun = fctMap.try_emplace(
                integ.front(),
                internals::CachingFunction<cmplx>(integ.front().toCallable()));
            groups.emplace_back(std::initializer_list{
                std::pair{size_t{0}, &(fun.first->second)}});
            integToIndex.emplace(integ.first(1), 0);
          }

          auto &grp = groups[headIter.first->second];

          size_t prevIdx = 0;
          for (size_t i = 2; i <= integ.size(); ++i) {
            std::span subInteg{integ.data(), i};
            auto elemIter = integToIndex.try_emplace(subInteg, grp.size());
            if (elemIter.second) { // new subintegral
              auto fun = fctMap.try_emplace(subInteg.back(),
                                            internals::CachingFunction<cmplx>(
                                                subInteg.back().toCallable()));
              grp.emplace_back(prevIdx, &(fun.first->second));
            }
            prevIdx = elemIter.first->second;
          }
        }
      }
    }

    std::vector<std::vector<cmplx>> integralVals;
    // execute integration
    for (const auto &grp : groups) {
      integralVals.push_back(evalTree(zMax, grp));
    }

    std::vector<cmplx> res;
    res.reserve(integrals.size());
    // construct result from integralValueMap and integration
    // drop const for span creation
    for (auto &regInt : regulatedIntegrals) {
      cmplx val = regInt.first;
      for (auto &integs : regInt.second) {
        auto iter = integralValueMap.find(integs.second);
        cmplx integVal;
        if (iter != integralValueMap.end()) {
          integVal = iter->second;
        } else {
          integVal =
              integralVals[headToIndex[integs.second.front()]][integToIndex[{
                  integs.second.data(), integs.second.size()}]];
        }

        val += integs.first * integVal;
      }
      res.push_back(val);
    }

    return res;
  }

  /**
   * computes a list of iterated integrals by combining them into trees
   * and only regulating them along a short segment.
   * @param zMax the upper integration bound
   * @param integrals a vector of the integrals represented by a vector of
   * pointers to Integrand<cmplx> objects.
   * @param checkIntegral a function that returns an optional when given a span
   * of TaggedFunctions. Its value, if present, replaces the passed in integral
   **/
  std::vector<cmplx> splittingTreeIntegrate(
      cmplx zMax, const std::vector<std::vector<Integrand<cmplx> *>> &integrals,
      auto checkIntegral)
    requires std::is_convertible_v<
        decltype(checkIntegral(std::span<TaggedFunction<cmplx>>{})),
        std::optional<cmplx>>
  {

    std::vector<std::vector<Integrand<cmplx> *>> subIntegralList;
    std::unordered_map<Integrand<cmplx> *, FunctionType> fctMap;
    std::vector<std::vector<std::pair<size_t, FunctionType *>>> groups;
    std::unordered_map<Integrand<cmplx> *, size_t> headToIndex;
    std::unordered_map<std::span<Integrand<cmplx> *const>, size_t,
                       internals::ContainerHash<Integrand<cmplx> *const>,
                       internals::SpanEquals<Integrand<cmplx> *const>>
        integToIndex;

    for (const auto &integ : integrals)
      for (const auto &integr : integ)
        fctMap.try_emplace(integr, internals::CachingFunction<cmplx>(*integr));

    // determine all necessary Subintegrals
    for (auto &integ : integrals) {
      // create new tree
      auto curTree = headToIndex.try_emplace(integ.front(), groups.size());
      if (curTree.second) {
        groups.emplace_back(std::initializer_list{
            std::pair<size_t, FunctionType *>{0, &fctMap[integ.front()]}});
        std::span<Integrand<cmplx> *const> span = std::span(integ.data(), 1);
        integToIndex.emplace(span, 0);
        subIntegralList.emplace_back(std::initializer_list{integ.front()});
      }

      size_t prevIdx = 0;
      for (size_t i = 2; i <= integ.size(); ++i) {
        auto &grp = groups[(curTree.first->second)];
        std::span subInteg{integ.data(), i};
        auto subIter = integToIndex.try_emplace(subInteg, grp.size());
        if (subIter.second) { // add new subintegral
          grp.emplace_back(prevIdx, &fctMap[integ[i - 1]]);
          subIntegralList.emplace_back(subInteg.begin(), subInteg.end());
        }
        prevIdx = subIter.first->second;
      }
    }

    // evaluate them using treeIntegrate

    std::vector<cmplx> initVals =
        treeIntegrate(switchPoint * zMax, subIntegralList, checkIntegral);

    std::vector<std::vector<cmplx>> integralState;
    integralState.reserve(groups.size());
    // initialize state with default values
    for (size_t i = 0; i < groups.size(); ++i) {
      integralState.emplace_back(groups[i].size());
    }

    // now we fill in the data from initVals
    for (size_t i = 0; i < initVals.size(); ++i) {
      integralState[headToIndex[subIntegralList[i].front()]][integToIndex[{
          subIntegralList[i].cbegin(), subIntegralList[i].cend()}]] =
          initVals[i];
    }

    // directly call the integration routine to update integralState
    for (size_t i = 0; i < groups.size(); ++i) {
      callTreeIntegrator(zMax, groups[i], integralState[i], switchPoint);
    }

    // finally extract the result from the calculated data
    std::vector<cmplx> res;
    res.reserve(integrals.size());

    for (const auto &integ : integrals) {
      res.push_back(
          integralState[headToIndex[integ.front()]][integToIndex[integ]]);
    }

    return res;
  }

  /*Virtual destructor*/
  virtual ~IntegratorBase() = default;

protected:
  /**
   * evaluates a given integral starting in zMax*startRegulator with zero
   * initial values. This method can be overloaded, but in general this should
   * not be necessary.
   * @param zMax the upper integration bound
   * @param integral the integral to directly evaluate
   * @returns the vector of all partial integrals evaluated at zMax
   **/
  virtual StateType baseEval(const cmplx &zMax,
                             const std::vector<FunctionType> &integral) {
    StateType state(integral.size());
    callBasicIntegrator(zMax, integral, state, startRegulator);
    return state;
  }

  /**
   * evaluates a single integral at zMax with given starting point and initial
   * values. It only returns the value of the full integral. state can be
   * modified. This method can be overloaded, but this is rarely necessary.
   * @param zMax the upper integration bound
   * @param integral the integral to be evaluated
   * @param state the initial state, which might be modified during execution
   * @param start the value of t between zero and one where the integration
   * should be started
   * @returns the value of the full integral
   */
  virtual cmplx evalIntegral(const cmplx &zMax,
                             const std::vector<FunctionType> &integral,
                             StateType &state, const nr &start) {
    callBasicIntegrator(zMax, integral, state, start);
    return state.back();
  }

  /**
   * evaluates the family of integrals passed in as a tree starting in
   * startRegulator*zMax with zero initial conditions.
   * @param zMax the upper integration bound
   * @param tree the tree of integrals to evaluate. The zeroth entry has no
   * parent. For the others the index of the parent is given by the first
   * element of the pair.
   * @returns a vector containing the value of the i-th integral in tree at
   * index i
   **/
  virtual std::vector<cmplx>
  evalTree(const cmplx &zMax,
           const std::vector<std::pair<size_t, FunctionType *>> &tree) {
    StateType state(tree.size());
    callTreeIntegrator(zMax, tree, state, startRegulator);
    return state;
  }

  /**
   * purely virtual function, which should compute the iterated integrals
   * corresponding to the given tree at zMax with given initial conditions and
   * starting point
   * @param zMax the upper integration bound
   * @param tree the tree of integrals to evaluate. The zeroth entry has no
   * parent. For the others the index of the parent is given by the first
   * element of the pair
   * @param state the initial state to use. Its indices match the parameter
   * tree. At the end it should contain the result
   * @param start the starting point between zero and one
   **/
  virtual void
  callTreeIntegrator(const cmplx &zMax,
                     const std::vector<std::pair<size_t, FunctionType *>> &tree,
                     StateType &state, const nr &start) = 0;

  /**
   * purely virtual function, which should compute the iterated integral
   * corresponding to the given vector at zMax with given initial conditions and
   * starting point
   * @param zMax the upper integration bound
   * @param integral the vector of integration kernels to use
   * @param state the initial state to use. Its indices match the parameter
   * tree. At the end it should contain the result
   * @param start the starting point between zero and one
   **/
  virtual void callBasicIntegrator(const cmplx &zMax,
                                   const std::vector<FunctionType> &integral,
                                   StateType &state, const nr &start) = 0;

protected:
  nr switchPoint;
  Regulator<nr, cmplx> reg;
  nr startRegulator;
};

} // namespace iteratedIntegrals
