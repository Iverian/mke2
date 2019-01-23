#include <csr_matrix.hpp>
#include <debug.hpp>
#include <export.hpp>
#include <global_eq_builder.hpp>
#include <local_eq_gen.hpp>
#include <triangulation.hpp>
#include <lup_factor.hpp>

#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

// #define DENSE_SOLVE

using namespace std;

static constexpr auto xdim = 200.;
static constexpr auto ydim = 40.;
static constexpr auto zdim = 40.;
static constexpr size_t scale = 10;

int main(int argc, char const* argv[])
{
    try {
        double init = 0;
        if (argc > 2) {
            init = atof(argv[2]);
        }
        auto t = Triangulation::cuboid({xdim, ydim, zdim}, scale);

        auto start = chrono::high_resolution_clock::now();

#ifndef DENSE_SOLVE
        auto [mat, vec] = build_global_system(t, LocalEqGen(gen_local));
        auto x0 = Vec(vec.size(), init);
        auto res = psolve(mat, vec, move(x0));
#else  // DENSE_SOLVE
        auto [mat, vec] = build_global_system_dense(t, LocalEqGen(gen_local));
        auto res = LupFactor(mat).factor().solve(vec);
#endif // DENSE_SOLVE

        auto finish = chrono::high_resolution_clock::now();
        cout << "Iteration finished in "
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

    } catch (const std::exception& err) {
        cerr << "CRITICAL " << err.what() << endl;
        std::exit(1);
    }

    return 0;
}
