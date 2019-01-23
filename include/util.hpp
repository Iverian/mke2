#ifndef MKE2_INCLUDE_UTIL_HPP_
#define MKE2_INCLUDE_UTIL_HPP_

#include "common_types.hpp"


enum class Tolerance { ZERO = 0, SINGLE = 1, DOUBLE = 2, TRIPLE = 3 };

bool isnear(Value lhs, Value rhs, Tolerance t = Tolerance::DOUBLE);
bool iszero(Value x, Tolerance t = Tolerance::DOUBLE);
Value sqr(Value x);

#endif // MKE2_INCLUDE_UTIL_HPP_