#ifndef MKE2_SRC_MKE2_GLOBAL_INDICES_HPP_
#define MKE2_SRC_MKE2_GLOBAL_INDICES_HPP_

#include <triangulation.hpp>

Triangulation::Index _g(const Triangulation::Index& node,
                        const Triangulation::Index& coord,
                        const Triangulation::Index& m);
Triangulation::Index _v(const Triangulation::Index& node,
                        const Triangulation::Index& coord);
Triangulation::Index _s(const Triangulation::Index& node,
                        const Triangulation::Index& coord);

enum Coord : Triangulation::Index { X, Y, Z };

#endif // MKE2_SRC_MKE2_GLOBAL_INDICES_HPP_