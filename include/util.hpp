#ifndef MKE2_INCLUDE_UTIL_HPP_
#define MKE2_INCLUDE_UTIL_HPP_

#include <array>

enum class Tolerance { ZERO = 0, SINGLE = 1, DOUBLE = 2, TRIPLE = 3 };

bool isnear(double lhs, double rhs, Tolerance t = Tolerance::DOUBLE);
double sqr(double x);

#endif // MKE2_INCLUDE_UTIL_HPP_