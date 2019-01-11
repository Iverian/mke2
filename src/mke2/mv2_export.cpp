#include <mv2_export.hpp>

using namespace std;

void mv2_export(ostream& os, const Triangulation& t, const Vec& values)
{
    auto m = t.nodes().size();
    os << m << 3 << " " << 3 << " Ux Uy Uz" << endl;

    for (auto& n : t.nodes()) {
        auto [p, i] = n;
        auto& x = values[i];
        auto& y = values[m + i];
        auto& z = values[2 * m + i];
        os << i << " " << p[0] << " " << p[1] << " " << p[2] << " " << x << " "
           << y << " " << z << endl;
    }
}