#include <gtest/gtest.h>

#include <debug.hpp>
#include <sparce_matrix.hpp>

using namespace std;

TEST(TestSparce, test_find)
{
    SparceMatrix m({3, 3});
    m.fetch_add(0, 0, 1);
    m.fetch_add(1, 1, 2);
    m.fetch_add(2, 2, 3);
    m.fetch_add(1, 0, 4);
    coutd << m << endl;

    ASSERT_DOUBLE_EQ(m(0, 0), 1);
    ASSERT_DOUBLE_EQ(m(1, 1), 2);
    ASSERT_DOUBLE_EQ(m(2, 2), 3);
    ASSERT_DOUBLE_EQ(m(1, 0), 4);
}