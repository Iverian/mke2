#include <gtest/gtest.h>

#include <global_eq_builder.hpp>
#include <global_indices.hpp>
#include <local_eq_gen.hpp>
#include <triangulation.hpp>
#include <util.hpp>

#include <memory>
#include <random>

using namespace std;

struct TestGlobalEq : ::testing::Test {
    struct Data {
        Data()
            : t(Triangulation::cuboid({200, 40, 40}, 4))
            , glob(build_global_system(t, LocalEqGen(gen_local)))
        {
        }

        Triangulation t;
        pair<CsrMatrix, Vec> glob;
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
TEST_F(TestGlobalEq, test_symmetric)
{
    auto& mat = d->glob.first;

    for (Index i = 0; i < mat.shape().first; ++i) {
        auto k = mat.indptr()[i];
        auto last = mat.indptr()[i + 1];

        for (; k < last; ++k) {
            auto j = mat.indices()[k];
            auto v = mat.data()[k];

            ASSERT_TRUE(isfinite(v));
            ASSERT_DOUBLE_EQ(v, mat(j, i));
        }
    }

    SUCCEED();
}

TEST_F(TestGlobalEq, test_zero_rows)
{
    auto& mat = d->glob.first;

    for (Index i = 0; i < mat.shape().first; ++i) {
        auto first = mat.indptr()[i];
        auto last = mat.indptr()[i + 1];

        ASSERT_TRUE(!iszero(mat.diag()[i]) || first < last);
    }

    SUCCEED();
}

TEST_F(TestGlobalEq, test_zero_cols)
{
    auto& mat = d->glob.first;
    auto m = mat.shape().first;
    vector<bool> nonzero(m, false);
    size_t count = 0;

    for (Index i = 0; i < m; ++i) {
        if (!iszero(mat.diag()[i])) {
            nonzero[i] = true;
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

TEST_F(TestGlobalEq, test_zero_diag)
{
    auto& mat = d->glob.first;
    auto m = mat.shape().first;

    for (Index i = 0; i < m; ++i) {
        ASSERT_TRUE(!iszero(mat.diag()[i]));
    }
    SUCCEED();
}

TEST_F(TestGlobalEq, test_positive)
{
    auto& mat = d->glob.first;
    auto m = mat.shape().first;

    default_random_engine eng;
    uniform_real_distribution<Value> dr(-1., 1.);

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
    SUCCEED();
}

TEST_F(TestGlobalEq, strict_first_boundary_condition)
{
    auto& t = d->t;
    auto& [mat, vec] = d->glob;
    const auto m = mat.shape().first;

    for (auto& fc : t.first()) {
        for (auto& n : fc) {
            auto gn = n.index();
            auto cn = t.on_first(n);

            for (Index p = 0; p < Triangulation::DIM; ++p) {
                auto cnp = t.on_first(cn, p);
                if (cnp) {
                    auto gnp = _g(gn, p, m);
                    auto a = mat.indptr()[gnp];
                    auto b = mat.indptr()[gnp + 1];

                    ASSERT_EQ(a, b);
                    ASSERT_DOUBLE_EQ(vec[gnp], 0.);
                    ASSERT_DOUBLE_EQ(mat.diag()[gnp], 1.);

                    for (Index j = 0; j < m; ++j) {
                        if (gnp == j) {
                            continue;
                        }
                        ASSERT_DOUBLE_EQ(mat(gnp, j), 0.);
                    }
                }
            }
        }
    }

    SUCCEED();
}