#include <gtest/gtest.h>

#include <debug.hpp>
#include <global_eq_builder.hpp>
#include <local_eq_v17.hpp>
#include <triangulation.hpp>

using namespace std;

TEST(TestTriangBasic, test_tetra)
{
    Triangulation t({1, 1, 1});
    auto p = t.append_nodes({Point3d {0, 0, 0}, Point3d {0, 0, 1},
                             Point3d {0, 1, 0}, Point3d {1, 0, 0}});
    t.append_elem(p[0], p[1], p[2], p[3]);

    auto [mat, vec] = build_global_system(t, LocalEqGen(v17));

    for (CsrMatrix::Index i = 0; i < mat.shape().m; ++i) {
        auto a = mat.indptr()[i];
        auto b = mat.indptr()[i + 1];
        ASSERT_EQ(b - a, 1);
        ASSERT_DOUBLE_EQ(mat.data()[a], 1);
        ASSERT_DOUBLE_EQ(vec[i], 0.);
    }

    SUCCEED();
}

TEST(TestTriangBasic, test_msh)
{
    auto t = Triangulation::from_msh("brick.msh", {200, 40, 40});
    SUCCEED();
}