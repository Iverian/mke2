#ifndef MKE2_INCLUDE_TRIANG_BUILDER_H_
#define MKE2_INCLUDE_TRIANG_BUILDER_H_

#include <unordered_map>
#include <vector>

#include "dense_matrix.hpp"
#include "point3d.hpp"

class Triangulation {
public:
    static constexpr AbstractMatrix::Index N = 4;
    static constexpr AbstractMatrix::Index SN = 3;

    struct OnFirst {
        bool on_sigma_1;
        bool on_sigma_2;

        operator bool()
        {
            return on_sigma_1 || on_sigma_2;
        }
    };

    using Index = AbstractMatrix::Index;
    using NodeContainer = std::unordered_map<Point3d, Index>;
    using NodePtr = NodeContainer::const_pointer;
    using FiniteElement = std::array<NodePtr, N>;
    using SurfaceElement = std::array<NodePtr, SN>;
    using FiniteElementData = std::array<Point3d, N>;
    using SurfaceElementData = std::array<Point3d, SN>;

    explicit Triangulation(const std::array<double, 3>& dim);

    const std::array<double, 3>& dim() const;
    const NodeContainer& nodes() const;
    const std::vector<SurfaceElement>& triangles() const;
    const std::vector<FiniteElement>& elems() const;
    const std::array<std::vector<NodePtr>, 2>& first() const;
    const std::vector<SurfaceElement>& third() const;

    NodePtr append_node(const Point3d& p);
    std::vector<NodePtr> append_nodes(const std::vector<Point3d>& vp);
    void append_elem(const FiniteElement& e);

    FiniteElementData data(const FiniteElement& e) const;
    SurfaceElementData data(const SurfaceElement& e) const;
    SurfaceElement face(const FiniteElement& e, Index i) const;

    bool is_boundary(const SurfaceElement& e) const;
    void extract_triangles();

    bool on_third(const SurfaceElement& e) const;
    OnFirst on_first(const NodePtr& n) const;

    static Triangulation cuboid(std::array<double, 3> dim, size_t scale);

private:
    std::array<double, 3> dim_;

    Index size_;
    NodeContainer nodes_;
    std::vector<FiniteElement> elems_;
    std::vector<SurfaceElement> triangles_;
    std::array<std::vector<NodePtr>, 2> first_;
    std::vector<SurfaceElement> third_;
};

namespace std {
template <size_t N>
struct hash<array<Triangulation::NodePtr, N>> {
    static constexpr size_t seed = 16651;
    using argument_type = array<Triangulation::NodePtr, N>;

    size_t operator()(const argument_type& key) const
    {
        size_t result = seed;
        for (auto i = 0; i < N; ++i) {
            result = (result << 1) ^ h(key[i]);
        }
        return result;
    }

private:
    hash<Triangulation::NodePtr> h;
};
}

#endif // MKE2_INCLUDE_TRIANG_BUILDER_H_