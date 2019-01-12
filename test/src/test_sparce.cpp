#include <gtest/gtest.h>

#include <debug.hpp>
#include <sparce_matrix.hpp>

using namespace std;

TEST(TestSparce, test_modify)
{
    SparceMatrix m({3, 3});

    m.fetch_add(0, 0, 1);
    m.fetch_add(0, 0, 2);
    m.fetch_sub(0, 0, 3);

    m.fetch_add(1, 0, 4);
    m.fetch_add(1, 1, 2);
    m.fetch_add(1, 2, 5);

    m.fetch_set(2, 0, 0);
    m.fetch_add(2, 2, 3);
    m.clean_up();

    coutd << m << endl;

    ASSERT_DOUBLE_EQ(m(0, 0), 0);
    ASSERT_DOUBLE_EQ(m(1, 0), 4);
    ASSERT_DOUBLE_EQ(m(1, 1), 2);
    ASSERT_DOUBLE_EQ(m(1, 2), 5);
    ASSERT_DOUBLE_EQ(m(2, 0), 0);
    ASSERT_DOUBLE_EQ(m(2, 2), 3);
}

TEST(TestSparce, test_solve_1)
{
    SparceMatrix m({5, 5}, {1, 0, 0, 0, 0, 0, 4, 0, 1, 3, 0, 0, 3,
                            0, 0, 2, 1, 0, 3, 1, 1, 0, 0, 1, 0});
    Vec x{1000, 2000, 100, 4000, 1};
    Vec x0{0, 0, 0, 0, 0};
    auto b = m * x;
    auto x1 = solve_cg(m, b, x0, 100);

    ASSERT_EQ(x, x1);
}

TEST(TestSparce, test_solve_2)
{
    SparceMatrix m({5, 5}, {1, 0, 0, 0, 0, 0, 4, 0, 1, 3, 0, 0, 3,
                            0, 0, 2, 1, 0, 3, 1, 1, 0, 0, 1, 0});
    Vec x{1000, 2000, 100, 4000, 1};
    Vec x0{0, 0, 0, 0, 0};
    auto b = m * x;
    auto x1 = solve_bcg(m, b, x0, 100);

    ASSERT_EQ(x, x1);
}
