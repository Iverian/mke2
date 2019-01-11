#include <export.hpp>
#include <global_eq_builder.hpp>
#include <local_eq_v17.hpp>
#include <sparce_matrix.hpp>
#include <triangulation.hpp>

#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace std;

static constexpr auto xdim = 40.;
static constexpr auto ydim = 40.;
static constexpr auto zdim = 200.;

int main(int argc, char const* argv[])
{
    auto t = Triangulation::cuboid({xdim, ydim, zdim}, 4);
    auto gen = make_shared<LocalEqV17>();

    GlobalEqBuilder build(t, gen);
    build.get();
    auto res
        = solve_cg(build.mat(), build.vec(), Vec(3 * t.nodes().size(), 0.));

    if (argc > 1) {
        ofstream of(argv[1]);
        mv2_export(of, t, res);
    } else {
        mv2_export(cout, t, res);
    }
    return 0;
}
