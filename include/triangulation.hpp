#ifndef MKE2_INCLUDE_TRIANG_BUILDER_H_
#define MKE2_INCLUDE_TRIANG_BUILDER_H_

#include <array>
#include <unordered_map>
#include <utility>
#include <vector>

#include "dense_matrix.hpp"
#include "point3d.hpp"

class Triangulation {
public:
    static constexpr AbstractMatrix::Index DIM = 3;
    static constexpr AbstractMatrix::Index N = DIM + 1;
    static constexpr AbstractMatrix::Index SN = DIM;

    struct OnFirst {
        bool on_sigma_1;
        bool on_sigma_2;

        operator bool() const noexcept
        {
            return on_sigma_1 || on_sigma_2;
        }
    };

    using Index = AbstractMatrix::Index;
    using NodeContainer = std::unordered_map<Point3d, Index>;

    class NodePtr {
    public:
        Index index() const;
        const Point3d& point() const;

        operator NodeContainer::const_pointer() const
        {
            return ptr;
        }

        NodeContainer::const_pointer ptr;
    };

    struct SurfaceElement : std::array<NodePtr, SN> {
        using Super = std::array<NodePtr, SN>;
        using Data = std::array<Point3d, SN>;

        SurfaceElement(const Super& data = Super());

        Data data() const;
    };

    struct FiniteElement : std::array<NodePtr, N> {
        using Super = std::array<NodePtr, N>;
        using Data = std::array<Point3d, N>;

        FiniteElement(const Super& data = Super());

        Data data() const;
        SurfaceElement face(Index i) const;
    };

    explicit Triangulation(const std::array<double, DIM>& dim);

    const std::array<double, DIM>& dim() const;
    const NodeContainer& nodes() const;
    const std::vector<SurfaceElement>& triangles() const;
    const std::vector<FiniteElement>& elems() const;
    const std::array<std::vector<NodePtr>, 2>& first() const;
    const std::vector<SurfaceElement>& third() const;

    NodePtr append_node(const Point3d& p);
    std::vector<NodePtr> append_nodes(const std::vector<Point3d>& vp);

    template <class... Args>
    void append_elem(Args... args);

    bool is_boundary(const SurfaceElement& e) const;
    void extract_triangles();

    bool on_third(const SurfaceElement& e) const;
    OnFirst on_first(const NodePtr& n) const;
    OnFirst coord_on_first(OnFirst node, Index coord) const;

    static Triangulation cuboid(std::array<double, DIM> dim, size_t scale);

private:
    std::array<double, DIM> dim_;

    Index size_;
    NodeContainer nodes_;
    std::vector<FiniteElement> elems_;
    std::vector<SurfaceElement> triangles_;
    std::array<std::vector<NodePtr>, 2> first_;
    std::vector<SurfaceElement> third_;
};

template <class... Args>
void Triangulation::append_elem(Args... args)
{
    elems_.emplace_back(FiniteElement({std::move(args)...}));
    auto& e = elems_.back();

    for (Index i = 0; i < N; ++i) {
        auto f = e.face(i);
        if (on_third(f)) {
            third_.emplace_back(move(f));
        }
    }
}

// struct Triangulation::NodePtr : Triangulation::NodeContainer::const_pointer
// {
//     Index index() const;
// };

// struct Triangulation::FiniteElement
//     : std::array<Triangulation::NodePtr, Triangulation::N> {
//     using Data = std::array<Point3d, N>;

//     Data data() const;
//     SurfaceElement face() const;
// };

// struct Triangulation::SurfaceElement
//     : std::array<Triangulation::NodePtr, Triangulation::SN> {
//     using Data = std::array<Point3d, N>;

//     Data data() const;
// };

namespace std {

template <size_t N>
struct hash<array<Triangulation::NodePtr, N>> {
    static constexpr size_t seed = 16651;
    using argument_type = array<Triangulation::NodePtr, N>;

    size_t operator()(const argument_type& key) const
    {
        size_t result = seed;
        for (size_t i = 0; i < N; ++i) {
            result = (result << 1) ^ h(key[i].ptr);
        }
        return result;
    }

private:
    hash<Triangulation::NodeContainer::const_pointer> h;
};

template <>
struct hash<Triangulation::FiniteElement> {
    using argument_type = Triangulation::FiniteElement;

    size_t operator()(const argument_type& key) const
    {
        return h(key);
    }

private:
    hash<array<Triangulation::NodePtr, Triangulation::N>> h;
};

template <>
struct hash<Triangulation::SurfaceElement> {
    using argument_type = Triangulation::SurfaceElement;

    size_t operator()(const argument_type& key) const
    {
        return h(key);
    }

private:
    hash<array<Triangulation::NodePtr, Triangulation::SN>> h;
};

}

#endif // MKE2_INCLUDE_TRIANG_BUILDER_H_