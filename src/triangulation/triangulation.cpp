#include <debug.hpp>
#include <triangulation.hpp>
#include <util.hpp>

#include <algorithm>
#include <iterator>
#include <unordered_set>

#include <mke2/global_indices.hpp>

using namespace std;

Triangulation::Triangulation(const array<double, 3>& dim)
    : dim_(dim)
    , size_(0)
    , nodes_()
    , elems_()
    , first_()
    , third_()
{
}

const array<double, 3>& Triangulation::dim() const
{
    return dim_;
}

const Triangulation::NodeContainer& Triangulation::nodes() const
{
    return nodes_;
}

const vector<Triangulation::SurfaceElement>& Triangulation::triangles() const
{
    return triangles_;
}

const vector<Triangulation::FiniteElement>& Triangulation::elems() const
{
    return elems_;
}

const array<vector<Triangulation::NodePtr>, 2>& Triangulation::first() const
{
    return first_;
}

const vector<Triangulation::SurfaceElement>& Triangulation::third() const
{
    return third_;
}

Triangulation::NodePtr Triangulation::append_node(const Point3d& p)
{
    auto [it, flag] = nodes_.insert({p, size_});
    auto result = NodePtr {&(*it)};
    if (flag) {
        ++size_;
        auto first_flag = on_first(result);
        if (first_flag.on_sigma_1) {
            first_[0].push_back(result);
        }
        if (first_flag.on_sigma_2) {
            first_[1].push_back(result);
        }
    }

    return result;
}

vector<Triangulation::NodePtr>
Triangulation::append_nodes(const vector<Point3d>& vp)
{
    vector<NodePtr> result;
    result.reserve(vp.size());
    transform(begin(vp), end(vp), back_inserter(result),
              [this](auto& p) { return append_node(p); });
    return result;
}

// void Triangulation::append_elem(const FiniteElement& e)
// {
//     elems_.push_back(e);
//     for (Index i = 0; i < N; ++i) {
//         auto f = e.face(i);
//         if (on_third(f)) {
//             third_.emplace_back(move(f));
//         }
//     }
// }

Triangulation::Index Triangulation::NodePtr::index() const
{
    return ptr->second;
}

const Point3d& Triangulation::NodePtr::point() const
{
    return ptr->first;
}

Triangulation::SurfaceElement::SurfaceElement(const Super& data)
    : Super(data)
{
}

Triangulation::SurfaceElement::Data Triangulation::SurfaceElement::data() const
{
    return {(*this)[0].point(), (*this)[1].point(), (*this)[2].point()};
}

Triangulation::SurfaceElement Triangulation::FiniteElement::face(Index i) const
{
    SurfaceElement face;
    for (Index j = 0; j < SN; ++j) {
        face[i] = (*this)[(i + j + 1) % N];
    }
    return face;
}

Triangulation::FiniteElement::FiniteElement(const Super& data)
    : Super(data)
{
}

Triangulation::FiniteElement::Data Triangulation::FiniteElement::data() const
{
    return {(*this)[0].point(), (*this)[1].point(), (*this)[2].point(),
            (*this)[3].point()};
}

bool Triangulation::on_third(const SurfaceElement& e) const
{
    const auto z = dim_[2];

    return isnear(e[0].point()[2], z) && isnear(e[1].point()[2], z)
        && isnear(e[2].point()[2], z);
}

Triangulation::OnFirst Triangulation::on_first(const NodePtr& n) const
{
    OnFirst result {false, false};

    auto& p = n.point();
    // Sigma_1
    if (isnear(p[0], 0) || isnear(p[0], dim_[0])) {
        result.on_sigma_1 = true;
    }
    // Sigma_2
    if (isnear(p[1], 0) || isnear(p[1], dim_[1])) {
        result.on_sigma_2 = true;
    }

    return result;
}

Triangulation::OnFirst Triangulation::coord_on_first(OnFirst node,
                                                     Index coord) const
{

    return {node.on_sigma_1 && (coord == Coord::X || coord == Coord::Z),
            node.on_sigma_2 && (coord == Coord::Y)};
}

bool Triangulation::is_boundary(const SurfaceElement& e) const
{
    auto d = e.data();

    for (Index i = 0; i < SN; ++i) {
        auto& v = dim_[i];
        bool one
            = isnear(d[0][i], 0) && isnear(d[1][i], 0) && isnear(d[2][i], 0);
        bool two
            = isnear(d[0][i], v) && isnear(d[1][i], v) && isnear(d[2][i], v);

        if (one || two) {
            return true;
        }
    }
    return false;
}

void Triangulation::extract_triangles()
{
    unordered_set<SurfaceElement> unique;

    for (auto& k : elems_) {
        for (size_t i = 0; i < N; ++i) {
            auto f = k.face(i);
            if (is_boundary(f)) {
                unique.insert(f);
            }
        }
    }

    triangles_.resize(unique.size());
    copy(begin(unique), end(unique), begin(triangles_));
}

Triangulation Triangulation::cuboid(array<double, DIM> dim, size_t scale)
{
    Triangulation result(dim);

    array<Index, DIM> m;
    array<double, DIM> s;
    auto mind = *min_element(begin(dim), end(dim));

    if (scale == 0) {
        scale = 1;
    }

    for (Index i = 0; i < DIM; ++i) {
        m[i] = scale * Index(dim[i] / mind);
        s[i] = dim[i] / m[i];
    }

    result.elems_.reserve(6 * m[0] * m[1] * m[2]);
    for (Index i = 0; i < m[0]; ++i) {
        auto x0 = i * s[0];
        auto x1 = (i + 1) * s[0];

        for (Index j = 0; j < m[1]; ++j) {
            auto y0 = j * s[1];
            auto y1 = (j + 1) * s[1];

            for (Index k = 0; k < m[2]; ++k) {
                auto z0 = k * s[2];
                auto z1 = (k + 1) * s[2];

                auto p = result.append_nodes(
                    combine({x0, x1}, {y0, y1}, {z0, z1}));
                result.append_elem(p[0], p[1], p[2], p[6]);
                result.append_elem(p[0], p[1], p[4], p[6]);
                result.append_elem(p[1], p[3], p[6], p[7]);
                result.append_elem(p[1], p[5], p[6], p[7]);
                result.append_elem(p[1], p[2], p[3], p[6]);
                result.append_elem(p[1], p[4], p[5], p[6]);
            }
        }
    }

    return result;
}
