#ifndef MKE2_INCLUDE_UTIL_H_
#define MKE2_INCLUDE_UTIL_H_

#include <array>

static constexpr std::array<double, 3> tol = {1e-1, 1e-5, 1e-9};

enum class Tolerance { SINGLE = 1, DOUBLE = 2 };

bool isnear(double lhs, double rhs, Tolerance t = Tolerance::SINGLE);
double sqr(double x);

#endif // MKE2_INCLUDE_UTIL_H_