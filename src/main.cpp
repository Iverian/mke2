#include <debug.hpp>
#include <export.hpp>
#include <global_eq_builder.hpp>
#include <local_eq_v17.hpp>
#include <sparce_matrix.hpp>
#include <triangulation.hpp>

#include <chrono>
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
    static constexpr double tol = 1e-2;

    try {
        auto t = Triangulation::cuboid({xdim, ydim, zdim}, 2);
        auto gen = make_shared<LocalEqV17>();

        GlobalEqBuilder build(t, gen);
        build.get();

        auto start = chrono::high_resolution_clock::now();
        auto x0 = Vec(3 * t.nodes().size(), 1000.);
        auto res = solve_cg(build.mat(), build.vec(), x0, tol);
        auto finish = chrono::high_resolution_clock::now();
        cdbg << chrono::duration_cast<chrono::milliseconds>(finish - start)
                    .count()
             << endl;

        if (argc > 1) {
            ofstream of(argv[1]);
            mv2_export(of, t, res);
        } else {
            mv2_export(cout, t, res);
        }
    } catch (const exception& err) {
        cerr << "CRITICAL: " << err.what() << endl;
        return 1;
    }
    return 0;
}
