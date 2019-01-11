#include <iostream>

#include <global_eq_builder.hpp>
#include <local_eq_v17.hpp>
#include <mv2_export.hpp>
#include <sparce_matrix.hpp>
#include <triangulation.hpp>

#include <fstream>
#include <memory>

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

    ofstream of("out/v17_dbg.mv2");
    mv2_export(of, t, res);

    return 0;
}
