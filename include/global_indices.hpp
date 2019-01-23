#ifndef MKE2_SRC_MKE2_GLOBAL_INDICES_HPP_
#define MKE2_SRC_MKE2_GLOBAL_INDICES_HPP_

#include "common_types.hpp"
#include "debug.hpp"


enum class Coord : Index { X, Y, Z };

template <class CoordIndex>
Index _g(const Index& node, const CoordIndex& coord, const Index& m)
{
    check_if(node < m && Index(coord) < Triangulation::DIM,
             "Index out of range");

#ifndef MKE2_GLOBAL_INDICES_DOF_FIRST
    return node * Triangulation::DIM + Index(coord);
#else  // MKE2_GLOBAL_INDICES_DOF_FIRST
    return Index(coord) * m + node;
#endif // MKE2_GLOBAL_INDICES_DOF_FIRST
}

template <class CoordIndex>
Index _v(const Index& node, const CoordIndex& coord)
{
    return _g(node, coord, Triangulation::N);
}

template <class CoordIndex>
Index _s(const Index& node, const CoordIndex& coord)
{
    return _g(node, coord, Triangulation::SN);
}

#endif // MKE2_SRC_MKE2_GLOBAL_INDICES_HPP_