#ifndef MKE2_INCLUDE_COMMON_TYPES_HPP_
#define MKE2_INCLUDE_COMMON_TYPES_HPP_

#include <cstdlib>
#include <ostream>
#include <utility>


using Value = double;
using Index = size_t;
using Index2d = std::pair<Index, Index>;

std::ostream& operator<<(std::ostream& os, const Index2d& obj);

#endif // MKE2_INCLUDE_COMMON_TYPES_HPP_