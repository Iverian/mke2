#include <gtest/gtest.h>

#include <debug.hpp>
#include <dense_matrix.hpp>
#include <sparce_matrix.hpp>

#include <chrono>
#include <random>

using namespace std;

static constexpr AbstractMatrix::Index N = 100;

TEST(TestSparce, test_modify)
{
    SparceMatrix a(N);
    DenseMatrix b({N, N}, 0);

    default_random_engine eng;
    uniform_real_distribution<double> dr;
    uniform_int_distribution<AbstractMatrix::Index> di(0, N - 1);

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
    a.remove_zeroes();

    for (AbstractMatrix::Index i = 0; i < N; ++i) {
        for (AbstractMatrix::Index j = 0; j < N; ++j) {
            if (!isnear(a(i, j), b(i, j))) {
                coutd << "seed=" << seed << ", i=" << i << ", j=" << j << endl;
            }
            ASSERT_DOUBLE_EQ(a(i, j), b(i, j));
        }
    }

    SUCCEED();
}

TEST(TestSparce, test_clean_up)
{
    SparceMatrix a(N);

    default_random_engine eng;
    uniform_real_distribution<double> dr;
    uniform_int_distribution<AbstractMatrix::Index> di(0, N - 1);

    eng.seed(default_random_engine::result_type(
        chrono::system_clock::now().time_since_epoch().count()));

    for (auto count = 0; count < 30; ++count) {
        auto value = dr(eng);
        auto i = di(eng);
        auto j = di(eng);

        a.add(i, j, value);
        a.sub(i, j, value);
    }
    a.remove_zeroes();

    ASSERT_EQ(a.non_zero(), 0);
}

TEST(TestSparce, test_solve_rand)
{
    SparceMatrix a(N);
    Vec x(N, 0.);

    default_random_engine eng;
    uniform_real_distribution<double> dr;
    uniform_int_distribution<AbstractMatrix::Index> di(0, N - 1);

    eng.seed(default_random_engine::result_type(
        chrono::system_clock::now().time_since_epoch().count()));
    for (AbstractMatrix::Index i = 0; i < N; ++i) {
        a.set(i, i, dr(eng));
        x[i] = dr(eng);
    }

    for (auto count = 0; count < 30; ++count) {
        auto value = dr(eng);
        auto i = di(eng);
        auto j = di(eng);

        a.add(i, j, value);
        a.add(j, i, value);
    }
    a.remove_zeroes();

    auto b = a * x;
    auto y = solve(a, b, Vec(N, 0.));

    ASSERT_EQ(x, y);
}

TEST(TestSparce, test_solve)
{
    SparceMatrix m(5, {1, 0, 0, 2, 0, 0, 4, 0, 1, 0, 0, 0, 3,
                       0, 0, 2, 1, 0, 3, 1, 0, 0, 0, 1, 1});
    Vec x {1000, 2000, 100, 4000, 1};
    Vec x0 {0, 0, 0, 0, 0};
    auto b = m * x;
    auto x1 = solve(m, b, x0);

    ASSERT_EQ(x, x1);
}