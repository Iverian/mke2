#include <gtest/gtest.h>

#include <csr_matrix.hpp>
#include <debug.hpp>
#include <dense_matrix.hpp>

#include <chrono>
#include <random>

using namespace std;

static constexpr Index N = 2000;

TEST(TestSparse, test_modify)
{
    DokMatrix a(N);
    DenseMatrix b({N, N}, 0);

    default_random_engine eng;
    uniform_real_distribution<Value> dr;
    uniform_int_distribution<Index> di(0, N - 1);

    auto seed = default_random_engine::result_type(
        chrono::system_clock::now().time_since_epoch().count() % 1000000);
    eng.seed(seed);

    for (auto count = 0; count < 30; ++count) {
        auto value = dr(eng);
        auto i = di(eng);
        auto j = di(eng);

        if (eng() % 2) {
            a.add(i, j, value);
            b(i, j) += value;
        } else {
            a.sub(i, j, value);
            b(i, j) -= value;
        }
    }

    for (Index i = 0; i < N; ++i) {
        for (Index j = 0; j < N; ++j) {
            ASSERT_DOUBLE_EQ(a(i, j), b(i, j));
        }
    }

    SUCCEED();
}
TEST(TestSparse, test_solve_rand)
{
    DokMatrix da(N);
    Vec x(N, 0.);

    default_random_engine eng;
    uniform_real_distribution<Value> dr(-1., 1.);
    uniform_int_distribution<Index> di(0, N - 1);

    eng.seed(default_random_engine::result_type(
        chrono::system_clock::now().time_since_epoch().count()));
    for (Index i = 0; i < N; ++i) {
        da.set(i, i, dr(eng));
        x[i] = dr(eng);
    }

    for (auto count = 0; count < 2 * N; ++count) {
        auto value = dr(eng);
        auto i = di(eng);
        auto j = di(eng);

        da.add(i, j, value);
        da.add(j, i, value);
    }

    CsrMatrix a(da);
    auto b = a * x;
    auto y = solve(a, b, Vec(N, 0.));

    ASSERT_EQ(x, y);
}

TEST(TestSparse, test_solve)
{
    CsrMatrix m(5, {1, 0, 0, 2, 0, 0, 4, 0, 1, 0, 0, 0, 3,
                    0, 0, 2, 1, 0, 3, 1, 0, 0, 0, 1, 1});
    Vec x {1000, 2000, 100, 4000, 1};
    Vec x0 {0, 0, 0, 0, 0};
    auto b = m * x;
    auto x1 = solve(m, b, x0);

    ASSERT_EQ(x, x1);
}