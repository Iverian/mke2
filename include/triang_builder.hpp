#ifndef MKE2_INCLUDE_TRIANG_BUILDER_H_
#define MKE2_INCLUDE_TRIANG_BUILDER_H_

#include <unordered_map>
#include <vector>

#include "dense_matrix.hpp"
#include "point3d.hpp"
#include "triangulation.hpp"

class TriangBuilder {
public:
    static constexpr size_t N = 4;

    using Index = size_t;
    using NodeContainer = std::unordered_map<Point3d, Index>;
    using NodePtr = NodeContainer::const_pointer;
    using FiniteElement = std::array<NodePtr, N>;
    using SurfaceElement = std::array<NodePtr, N - 1>;
    using FiniteElementData = std::array<Point3d, N>;
    using SurfaceElementData = std::array<Point3d, N - 1>;

    explicit TriangBuilder(const std::array<double, 3>& dim);

    const NodeContainer& nodes() const;
    const std::vector<FiniteElement>& elems() const;
    const std::array<std::vector<NodePtr>, 2>& first() const;
    const std::vector<SurfaceElement>& third() const;

    NodePtr append_node(const Point3d& p);
    std::vector<NodePtr> append_nodes(const std::vector<Point3d>& vp);
    void append_elem(const FiniteElement& e);

    FiniteElementData data(const FiniteElement& e) const;
    SurfaceElementData data(const SurfaceElement& e) const;
    SurfaceElement face(const FiniteElement& e, Index i) const;

    static TriangBuilder cuboid(const std::array<double, 3>& dim,
                                size_t scale);

protected:
    bool on_third(const SurfaceElement& e) const;

private:
    std::array<double, 3> dim_;

    NodeContainer nodes_;
    std::vector<FiniteElement> elems_;
    std::array<std::vector<NodePtr>, 2> first_;
    std::vector<SurfaceElement> third_;
};

#endif // MKE2_INCLUDE_TRIANG_BUILDER_H_