#include <iostream>

#include <triangulation.hpp>

using namespace std;

static constexpr auto xdim = 40.;
static constexpr auto ydim = 40.;
static constexpr auto zdim = 200.;

int main(int argc, char const* argv[])
{
    auto mesh = triangulate_cuboid({xdim, ydim, zdim}, 4);
    return 0;
}
