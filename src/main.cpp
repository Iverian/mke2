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

// #define THIRD_SEP

using namespace std;

static constexpr auto xdim = 200.;
static constexpr auto ydim = 40.;
static constexpr auto zdim = 40.;
static constexpr size_t scale = 4;

int main(int argc, char const* argv[])
{
    try {
        double init = 0;
        if (argc > 1) {
            init = atof(argv[1]);
        }

        auto t = Triangulation::cuboid({xdim, ydim, zdim}, scale);

#ifdef THIRD_SEP
        LocalEqV17 gen;
        auto glob = build_global_system(t, gen);
#else
        auto glob = build_global_system(t, LocalEqGen(v17));
#endif
        auto start = chrono::high_resolution_clock::now();

        auto x0 = Vec(glob.second.size(), init);
        auto res = solve(glob.first, glob.second, move(x0));

        auto finish = chrono::high_resolution_clock::now();
        cout << "Iteration finished in "
             << chrono::duration_cast<chrono::milliseconds>(finish - start)
                    .count()
             << "ms" << endl;

        t.extract_triangles();
        if (argc > 2) {
            ofstream of(argv[2]);
            mv2_export(of, t, res);
        } else {
            mv2_export(cout, t, res);
        }

    } catch (const std::exception& err) {
        cerr << "CRITICAL " << err.what() << endl;
        std::exit(1);
    }

    return 0;
}
