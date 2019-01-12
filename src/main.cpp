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

static constexpr auto xdim = 200.;
static constexpr auto ydim = 40.;
static constexpr auto zdim = 40.;

int main(int argc, char const* argv[])
{
    try {
        auto t = Triangulation::cuboid({xdim, ydim, zdim}, 4);
        auto gen = make_shared<LocalEqV17>();

        GlobalEqBuilder build(t, gen);
        build.get();

        auto start = chrono::high_resolution_clock::now();

        auto x0 = Vec(3 * t.nodes().size(), 1000.);
        auto res = solve_bcg(build.mat(), build.vec(), x0);

        auto finish = chrono::high_resolution_clock::now();
        cerrd << "Iteration finished in "
              << chrono::duration_cast<chrono::milliseconds>(finish - start)
                     .count()
              << "ms" << endl;

        t.extract_triangles();
        if (argc > 1) {
            ofstream of(argv[1]);
            mv2_export(of, t, res);
        } else {
            mv2_export(cout, t, res);
        }

        if (argc > 2) {
            ofstream of(argv[2]);
            csv_export(of, t, res);
        }
    } catch (const std::exception& err) {
        cerr << "CRITICAL " << err.what() << endl;
        std::exit(1);
    }

    return 0;
}
