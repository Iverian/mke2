#ifndef MKE2_INCLUDE_TRIANGULATION_HPP_
#define MKE2_INCLUDE_TRIANGULATION_HPP_

#include <array>
#include <memory>
#include <unordered_map>
#include <vector>

#include "dense_matrix.hpp"
#include "point3d.hpp"

class Triang3d {
    struct NodeImpl;

public:
    struct Node;

    using NodePtr = std::shared_ptr<Node>;
    using NodeContainer = std::vector<NodePtr>;
    using Index = NodeContainer::size_type;

    using IndexSet = std::vector<Index>;
    using PointContainer = std::unordered_map<Point3d, IndexSet>;
    using PointPtr = PointContainer::pointer;

    static constexpr Index N = 4;
    static constexpr Index npos = Index(-1);

    void append(const std::array<Point3d, N> vertices);
    void set_adjacent();

    const NodeContainer& nodes() const;
    void reserve(Index node_count);
    void shrink_to_fit();

    friend std::ostream& operator<<(std::ostream& os, const Triang3d& obj);

protected:
    Index find_adjacent(const IndexSet& a, const IndexSet& b,
                        const IndexSet& c, Index cur);

private:
    PointContainer vertices_;
    NodeContainer nodes_;
};

struct Triang3d::Node {
    friend class Triang3d;
    explicit Node(const Triang3d& parent);

    const Point3d& v(Index i) const;
    NodePtr adj(Index i) const;
    Index iadj(Index i) const;

    std::vector<Index> boundary() const;
    std::array<Point3d, N - 1> face(Index i) const;
    Point3d normal(Index i) const;
    DenseMatrix vertex_view() const;

private:
    const Triang3d* parent_;
    std::array<PointPtr, N> vindices_;
    std::array<Index, N> adjacent_;
};

Triang3d triangulate_cuboid(std::array<double, 3> dim, size_t scale);

#endif // MKE2_INCLUDE_TRIANGULATION_HPP_