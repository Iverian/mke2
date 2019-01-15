#include <export.hpp>

#include "global_indices.hpp"

#include <tuple>

using namespace std;

void mv2_export(ostream& os, const Triangulation& t, const Vec& values)
{
    auto m = t.nodes().size();
    auto n = t.triangles().size();

    vector<Point3d> nodes(m);
    for (auto& n : t.nodes()) {
        auto [p, i] = n;
        nodes[i] = p;
    }

    os << m << " 3 3 Ux Uy Uz " << endl;
    for (Triangulation::Index i = 0; i < m; ++i) {
        auto& p = nodes[i];
        auto [x, y, z]
            = tie(values[_g(i, Coord::X, m)], values[_g(i, Coord::Y, m)],
                  values[_g(i, Coord::Z, m)]);

        os << (i + 1) << " " << p[0] << " " << p[1] << " " << p[2] << " ";
        os << x << " " << y << " " << z << endl;
    }

    os << n << " 3 3 BC_id mat_id mat_id_Out" << endl;
    auto& tr = t.triangles();
    for (Triangulation::Index i = 0; i < n; ++i) {
        auto& e = tr[i];
        os << (i + 1) << " " << (e[0].index() + 1) << " " << (e[1].index() + 1)
           << " " << (e[2].index() + 1) << " 1 1 0" << endl;
    }
}
