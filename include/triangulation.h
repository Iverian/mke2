#ifndef MKE2_INCLUDE_TRIANGULATION_H_
#define MKE2_INCLUDE_TRIANGULATION_H_

#include <array>
#include <memory>
#include <vector>

#include "tetrahedron.h"

class Triangulation : public std::vector<Tetrahedron> {
    using Super = std::vector<Tetrahedron>;

public:
};

Triangulation triangulate_cuboid(const std::array<double, 3> dim,
                                 size_t scale);

#endif // MKE2_INCLUDE_TRIANGULATION_H_