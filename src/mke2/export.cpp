#include <export.hpp>

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
    for (size_t i = 0; i < m; ++i) {
        auto& p = nodes[i];
        auto [x, y, z] = tie(values[i], values[m + i], values[2 * m + i]);
        os << (i + 1) << " " << p[0] << " " << p[1] << " " << p[2] << " ";
        os << x << " " << y << " " << z << endl;
    }

    os << n << " 3 3 BC_id mat_id mat_id_Out" << endl;
    auto& tr = t.triangles();
    for (size_t i = 0; i < n; ++i) {
        auto& e = tr[i];
        os << (i + 1) << " " << (e[0]->second + 1) << " " << (e[1]->second + 1)
           << " " << (e[2]->second + 1) << " 1 1 0" << endl;
    }
}

void csv_export(std::ostream& os, const Triangulation& t, const Vec& values)
{
    static constexpr char delim = ',';

    auto m = t.nodes().size();

    os << "x" << delim << "y" << delim << "z" << delim << "Ux" << delim << "Uy"
       << delim << "Uz" << endl;
    for (auto& n : t.nodes()) {
        auto [p, i] = n;
        auto [x, y, z] = tie(values[i], values[m + i], values[2 * m + i]);

        os << p[0] << delim << p[1] << delim << p[2] << delim << x << delim
           << y << delim << z << endl;
    }
}
