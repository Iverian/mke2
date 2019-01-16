#include <debug.hpp>
#include <global_indices.hpp>

#define NODE_FIRST

using Index = Triangulation::Index;

#ifdef NODE_FIRST

Index _g(const Index& node, const Index& coord, const Index& m)
{
    check_if(node < m && coord < Triangulation::DIM, "Index out of range");

    return node * Triangulation::DIM + coord;
}

#else

Index _g(const Index& node, const Index& coord, const Index& m)
{
    check_if(node < m && coord < Triangulation::DIM, "Index out of range");

    return coord * m + node;
}

#endif // NODE_FIRST

Index _v(const Index& node, const Index& coord)
{
    return _g(node, coord, Triangulation::N);
}

Index _s(const Index& node, const Index& coord)
{
    return _g(node, coord, Triangulation::SN);
}