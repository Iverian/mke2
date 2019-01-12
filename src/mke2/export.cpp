#include <export.hpp>

#include <tuple>

using namespace std;

void mv2_export(ostream& os, const Triangulation& t, const Vec& values)
{
    auto m = t.nodes().size();
    auto n = t.elems().size();
    auto q = t.third().size();

    vector<Point3d> nodes(m);
    for (auto& n : t.nodes()) {
        auto [p, i] = n;
        nodes[i] = p;
    }

    os << m << " " << 3 << " " << 3 << " Ux Uy Uz" << endl;
    for (size_t i = 0; i < m; ++i) {
        auto& p = nodes[i];
        auto [x, y, z] = tie(values[i], values[m + i], values[2 * m + i]);
        os << (i + 1) << " " << p[0] << " " << p[1] << " " << p[2] << " " << x
           << " " << y << " " << z << endl;
    }

    os << n << " " << 4 << " " << 4 << " A B C D" << endl;
    auto& elems = t.elems();
    for (size_t i = 0; i < elems.size(); ++i) {
        auto& e = elems[i];
        os << (i + 1) << " " << (e[0]->second + 1) << " " << (e[1]->second + 1)
           << " " << (e[2]->second + 2) << " " << (e[3]->second + 3) << " ";
        os << "0 0 0 0" << endl;
    }

    os << q << " " << 3 << endl;
    for (auto& e : t.third()) {
        os << (e[0]->second + 1) << " " << (e[1]->second + 1) << " "
           << (e[2]->second + 1) << endl;
    }
}

void csv_export(std::ostream& os, const Triangulation& t, const Vec& values)
{
    auto m = t.nodes().size();

    os << "x; y; z; Ux; Uy; Uz" << endl;
    for (auto& n : t.nodes()) {
        auto [p, i] = n;
        auto [x, y, z] = tie(values[i], values[m + i], values[2 * m + i]);

        os << p[0] << "; " << p[1] << "; " << p[2] << "; " << x << "; " << y
           << "; " << z << endl;
    }
}
