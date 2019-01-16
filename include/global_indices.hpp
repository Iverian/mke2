#ifndef MKE2_SRC_MKE2_GLOBAL_INDICES_HPP_
#define MKE2_SRC_MKE2_GLOBAL_INDICES_HPP_

#include "debug.hpp"
#include "triangulation.hpp"

enum class Coord : Triangulation::Index { X, Y, Z };

template <class CoordIndex>
Triangulation::Index _g(const Triangulation::Index& node,
                        const CoordIndex& coord, const Triangulation::Index& m)
{
    check_if(node < m && Triangulation::Index(coord) < Triangulation::DIM,
             "Index out of range");

#ifndef MKE2_GLOBAL_INDICES_DOF_FIRST
    return node * Triangulation::DIM + Triangulation::Index(coord);
#else  // MKE2_GLOBAL_INDICES_DOF_FIRST
    return Triangulation::Index(coord) * m + node;
#endif // MKE2_GLOBAL_INDICES_DOF_FIRST
}

template <class CoordIndex>
Triangulation::Index _v(const Triangulation::Index& node,
                        const CoordIndex& coord)
{
    return _g(node, coord, Triangulation::N);
}

template <class CoordIndex>
Triangulation::Index _s(const Triangulation::Index& node,
                        const CoordIndex& coord)
{
    return _g(node, coord, Triangulation::SN);
}

#endif // MKE2_SRC_MKE2_GLOBAL_INDICES_HPP_