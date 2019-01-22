#ifndef MKE2_INCLUDE_TRIANG_BUILDER_H_
#define MKE2_INCLUDE_TRIANG_BUILDER_H_

#include <array>
#include <iostream>
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

    class OnFirst : public std::array<bool, 2> {
    public:
        operator bool() const noexcept
        {
            return (*this)[0] || (*this)[1];
        }
    };

    using Index = AbstractMatrix::Index;
    using NodeContainer = std::unordered_map<Point3d, Index>;

    class NodePtr {
    public:
        Index index() const;
        const Point3d& point() const;
        NodeContainer::const_pointer ptr;
    };

    struct SurfaceElement : std::array<NodePtr, SN> {
        using Super = std::array<NodePtr, SN>;
        using Data = std::array<Point3d, SN>;

        SurfaceElement(const Super& data = Super());

        Data data() const;
        double area() const;

        friend std::ostream& operator<<(std::ostream& os,
                                        const SurfaceElement& obj);
    };

    struct FiniteElement : std::array<NodePtr, N> {
        using Super = std::array<NodePtr, N>;
        using Data = std::array<Point3d, N>;

        FiniteElement(const Super& data = Super());

        Data data() const;
        double volume() const;
        SurfaceElement face(Index i) const;

        friend std::ostream& operator<<(std::ostream& os,
                                        const FiniteElement& obj);
    };

    Triangulation();
    explicit Triangulation(const std::array<double, DIM>& dim);

    const std::array<double, DIM>& dim() const;
    const NodeContainer& nodes() const;
    const std::vector<std::pair<Index, SurfaceElement>>& triangles() const;
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
    OnFirst on_first(OnFirst cond, Index coord) const;

    friend std::ostream& operator<<(std::ostream& os,
                                    const Triangulation& obj);

    static Triangulation cuboid(std::array<double, DIM> dim,
                                std::array<Index, DIM> size);
    static Triangulation cuboid(std::array<double, DIM> dim, Index scale);
    static Triangulation from_msh(const char* filename,
                                  std::array<double, DIM> dim);

private:
    std::array<double, DIM> dim_;

    Index size_;
    NodeContainer nodes_;
    std::vector<FiniteElement> elems_;
    std::vector<std::pair<Index, SurfaceElement>> triangles_;
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