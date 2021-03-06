#include <dense_matrix.hpp>
#include <gtest/gtest.h>
#include <lup_factor.hpp>
#include <vec.hpp>


#include <random>

using namespace std;

static constexpr Index N = 525 * 3;

TEST(TestLup, test_solve)
{
    default_random_engine eng;
    uniform_real_distribution<Value> dr;

    DenseMatrix mat({N, N}, 0.);
    Vec x(N, 0.);

    for (Index i = 0; i < N; ++i) {
        x[i] = dr(eng);
        for (Index j = 0; j < N; ++j) {
            mat(i, j) = dr(eng);
        }
    }

    auto b = mat * x;
    auto x1 = LupFactor(mat).factor().solve(b);

    ASSERT_EQ(x, x1);
}

TEST(TestLup, test_inverse)
{
    default_random_engine eng;
    uniform_real_distribution<Value> dr;

    DenseMatrix mat({N, N}, 0.);
    auto e3 = DenseMatrix::eye(N);

    for (Index i = 0; i < N; ++i) {
        for (Index j = 0; j < N; ++j) {
            mat(i, j) = dr(eng);
        }
    }

    auto inv = LupFactor(mat).factor().inverse();

    ASSERT_EQ(mat * inv, e3);
    ASSERT_EQ(inv * mat, e3);
}

TEST(TestLup, test_det)
{
    default_random_engine eng;
    uniform_real_distribution<Value> dr;

    DenseMatrix mat({N, N}, 0.);

    Value det = 1;
    for (Index i = 0; i < N; ++i) {
        for (Index j = i; j < N; ++j) {
            mat(i, j) = dr(eng);
        }
        det *= mat(i, i);
    }

    auto det1 = LupFactor(mat).factor().det();
    ASSERT_DOUBLE_EQ(det, det1);
}