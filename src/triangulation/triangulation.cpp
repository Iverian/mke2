#include <point3d.h>
#include <tetrahedron.h>
#include <triangulation.h>
#include <util.h>

#include <cmath>

using namespace std;

Triangulation triangulate_cuboid(const array<double, 3> dim, size_t scale)
{
    Triangulation result;

    array<size_t, 3> m;
    array<double, 3> s;
    auto d = *min_element(begin(dim), end(dim));

    for (size_t i = 0; i < 3; ++i) {
        m[i] = max(2ull, scale * size_t(dim[i] / d));
        s[i] = dim[i] / (m[i] - 1);
    }

    result.reserve(6 * m[0] * m[1] * m[2]);
    for (size_t i = 0; i < m[0]; ++i) {
        auto x0 = i * s[0];
        auto x1 = (i + 1) * s[0];
        for (size_t j = 0; j < m[1]; ++j) {
            auto y0 = i * s[1];
            auto y1 = (i + 1) * s[1];
            for (size_t k = 0; k < m[2]; ++k) {
                auto z0 = i * s[2];
                auto z1 = (i + 1) * s[2];

                auto p = combine({x0, x1}, {y0, y1}, {z0, z1});
                result.push_back({p[0], p[1], p[2], p[6]});
                result.push_back({p[0], p[1], p[4], p[6]});
                result.push_back({p[1], p[3], p[6], p[7]});
                result.push_back({p[1], p[5], p[6], p[7]});
                result.push_back({p[1], p[2], p[3], p[6]});
                result.push_back({p[1], p[4], p[5], p[6]});
            }
        }
    }
    result.shrink_to_fit();

    return result;
}
