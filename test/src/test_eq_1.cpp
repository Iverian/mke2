#include <gtest/gtest.h>

#include <global_eq_builder.hpp>
#include <local_eq_v17.hpp>
#include <triangulation.hpp>
#include <util.hpp>

#include <memory>
#include <random>

using namespace std;

struct TestEq1 : ::testing::Test {
    struct Data {
        Data()
            : t(Triangulation::cuboid({200, 40, 40}, 4))
            , gen()
            , glob(build_global_system(t, gen))
        {
        }

        Triangulation t;
        LocalEqV17 gen;
        pair<SparceMatrix, Vec> glob;
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

TEST_F(TestEq1, test_symmetric)
{
    auto& mat = d->glob.first;

    for (SparceMatrix::Index i = 0; i < mat.shape().m; ++i) {
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

TEST_F(TestEq1, test_zero_rows)
{
    auto& mat = d->glob.first;

    for (SparceMatrix::Index i = 0; i < mat.shape().m; ++i) {
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

TEST_F(TestEq1, test_zero_cols)
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

TEST_F(TestEq1, test_positive)
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