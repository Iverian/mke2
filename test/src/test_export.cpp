#include <gtest/gtest.h>

#include <debug.hpp>
#include <export.hpp>
#include <triangulation.hpp>
#include <vec.hpp>

#include <fstream>

using namespace std;

TEST(TestExport, test_cube)
{
    auto t = Triangulation::cuboid({1, 1, 1}, 1);
    t.extract_triangles();
    auto v = Vec(Triangulation::DIM * t.nodes().size(), 0);
    ofstream of("cube.mv2");
    mv2_export(of, t, v);
}

TEST(TestExport, test_target)
{
    auto t = Triangulation::cuboid({200, 40, 40}, 4);
    t.extract_triangles();
    auto v = Vec(Triangulation::DIM * t.nodes().size(), 0);
    ofstream of("target.mv2");
    mv2_export(of, t, v);
}