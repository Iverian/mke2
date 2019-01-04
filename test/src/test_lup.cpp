#include <gtest/gtest.h>

#include <lup_factor.h>

#include <array>
#include <memory>

using namespace std;

struct Data {
    Data()
        : mat({3, 3}, {1, 3, 6, 5, 3, 8, 5, 7, 1})
        , lup(mat)
    {
    }
    DenseMatrix mat;
    LupFactor lup;
};

struct TestLup : ::testing::Test {
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

TEST_F(TestLup, test_factor)
{
    DenseMatrix r({3, 3}, {5, 3, 8, 1, 4, -7, 0.2, 0.6, 8.6});
    ASSERT_EQ(d->lup.mat(), r);
    ASSERT_EQ(d->lup.pivot(), LupFactor::PivotType({1, 2, 0, 5}));
}

TEST_F(TestLup, test_solve)
{
    array<Vec, 4> x = {Vec{1, 0, 0}, Vec{0, 1, 0}, Vec{0, 0, 1}, Vec{1, 1, 1}};
    for (auto& i : x) {
        ASSERT_EQ(d->lup.solve(d->mat * i), i);
    }
}

TEST_F(TestLup, test_inverse)
{
    DenseMatrix inv({3, 3},
                    {-0.30813953, 0.22674419, 0.03488372, 0.20348837,
                     -0.16860465, 0.12790698, 0.11627907, 0.04651163,
                     -0.06976744});
    ASSERT_EQ(d->lup.invert(), inv);
}

TEST_F(TestLup, test_det)
{
    ASSERT_DOUBLE_EQ(d->lup.det(), 172.);
}
