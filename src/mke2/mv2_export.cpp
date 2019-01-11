#include <mv2_export.hpp>

#include <fmt/ostream.h>

using namespace std;

void mv2_export(ostream& os, const Triangulation& t, const Vec& values)
{
    auto m = t.nodes().size();
    fmt::print(os, "{} 3 3 Ux Uy Uz\n", m);

    for (auto& n : t.nodes()) {
        auto [p, i] = n;
        auto& x = values[i];
        auto& y = values[m + i];
        auto& z = values[2 * m + i];
        fmt::print(os, "{:3d} {:.6f} {:.6f} {:.6f} {:.9f} {:.9f} {:.9f}", i,
                   p[0], p[1], p[2], x, y, z);
    }
}