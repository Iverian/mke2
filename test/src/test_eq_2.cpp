#include <gtest/gtest.h>

#include <global_eq_builder.hpp>
#include <global_indices.hpp>
#include <local_eq_v17.hpp>
#include <triangulation.hpp>
#include <util.hpp>

#include <memory>
#include <random>

using namespace std;

struct TestEq2 : ::testing::Test {
    struct Data {
        Data()
            : t(Triangulation::cuboid({200, 40, 40}, 4))
            , glob(build_global_system(t, LocalEqGen(v17)))
        {
        }

        Triangulation t;
        pair<SparseMatrix, Vec> glob;
    };

protected:
    void SetUp()
    {
        d = make_unique<Data>();
    }

    void TearDown()
    {
        d.reset();
    }

    unique_ptr<Data> d;
};

TEST_F(TestEq2, test_symmetric)
{
    auto& mat = d->glob.first;

    for (SparseMatrix::Index i = 0; i < mat.shape().m; ++i) {
        auto k = mat.indptr()[i];
        auto last = mat.indptr()[i + 1];

        for (; k < last; ++k) {
            auto j = mat.indices()[k];
            auto v = mat.data()[k];

            ASSERT_TRUE(isfinite(v));
            ASSERT_TRUE(isnear(v, mat(j, i)));
        }
    }

    SUCCEED();
}

TEST_F(TestEq2, test_zero_rows)
{
    auto& mat = d->glob.first;

    for (SparseMatrix::Index i = 0; i < mat.shape().m; ++i) {
        auto first = mat.indptr()[i];
        auto last = mat.indptr()[i + 1];

        auto cond = (first < last) || !isnear(mat.diag()[i], 0);
        if (!cond) {
            coutd << "i=" << i;
        }
        ASSERT_TRUE(cond);
    }

    SUCCEED();
}

TEST_F(TestEq2, test_zero_cols)
{
    auto& mat = d->glob.first;
    auto m = mat.shape().m;
    vector<bool> nonzero(m, false);
    size_t count = 0;

    for (AbstractMatrix::Index j = 0; j < m; ++j) {
        if (!isnear(mat.diag()[j], 0)) {
            nonzero[j] = true;
            ++count;
        }
    }

    for (auto& j : mat.indices()) {
        if (count >= m) {
            break;
        }

        if (!nonzero[j]) {
            nonzero[j] = true;
            ++count;
        }
    }

    ASSERT_EQ(count, m);
}

TEST_F(TestEq2, test_positive)
{
    auto& mat = d->glob.first;
    auto m = mat.shape().m;

    default_random_engine eng;
    uniform_real_distribution<double> dr(-1., 1.);

    Vec x(m, 0.), y(m, 0.);
    for (size_t count = 0; count < 100; ++count) {
        for (auto& xi : x) {
            xi = dr(eng);
        }
        dot(y, mat, x);

        auto res = dot(x, y);
        if (res < 0) {
            coutd << "count=" << count << ", res=" << res << endl;
        }
        ASSERT_GT(res, 0);
    }
}

TEST_F(TestEq2, strict_first_boundary_condition)
{
    auto& t = d->t;
    auto& [mat, vec] = d->glob;
    const auto m = mat.shape().m;

    for (auto& n : t.first()[0]) {
        auto gn = n.index();

        for (auto p : {Coord::X, Coord::Z}) {
            auto gnp = _g(gn, p, m);

            auto a = mat.indptr()[gnp];
            auto b = mat.indptr()[gnp + 1];
            ASSERT_EQ(a, b);
            ASSERT_DOUBLE_EQ(mat.diag()[gnp], 1);

            for (SparseMatrix::Index j = 0; j < m; ++j) {
                if (gnp == j) {
                    continue;
                }
                ASSERT_DOUBLE_EQ(mat(gnp, j), 0.);
            }
        }
    }

    for (auto& n : t.first()[1]) {
        auto gn = n.index();

        for (auto p : {Coord::Y}) {
            auto gnp = _g(gn, p, m);

            auto a = mat.indptr()[gnp];
            auto b = mat.indptr()[gnp + 1];
            ASSERT_EQ(a, b);
            ASSERT_DOUBLE_EQ(mat.diag()[gnp], 1);

            for (SparseMatrix::Index j = 0; j < m; ++j) {
                if (gnp == j) {
                    continue;
                }
                ASSERT_DOUBLE_EQ(mat(gnp, j), 0.);
            }
        }
    }
}