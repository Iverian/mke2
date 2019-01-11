#include <debug.hpp>
#include <triangulation.hpp>
#include <util.hpp>

#include <algorithm>

using namespace std;

Triangulation::Triangulation(const array<double, 3>& dim)
    : dim_(dim)
    , nodes_()
    , elems_()
    , first_()
    , third_()
{
}

const Triangulation::NodeContainer& Triangulation::nodes() const
{
    return nodes_;
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
    auto [it, flag] = nodes_.insert({p, nodes_.size()});
    auto result = &(*it);
    if (flag) {
        if (isnear(p[0], 0) || isnear(p[0], dim_[0])) {
            first_[0].push_back(result);
        }
        if (isnear(p[1], 0) || isnear(p[1], dim_[1])) {
            first_[1].push_back(result);
        }
    }

    return &(*it);
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

void Triangulation::append_elem(const FiniteElement& e)
{
    elems_.push_back(e);
    for (Index i = 0; i < N; ++i) {
        auto f = face(e, i);
        if (on_third(f)) {
            third_.emplace_back(move(f));
        }
    }
}

Triangulation::FiniteElementData
Triangulation::data(const FiniteElement& e) const
{
    return {e[0]->first, e[1]->first, e[2]->first, e[3]->first};
}

Triangulation::SurfaceElementData
Triangulation::data(const SurfaceElement& e) const
{
    return {e[0]->first, e[1]->first, e[2]->first};
}

Triangulation::SurfaceElement Triangulation::face(const FiniteElement& e,
                                                  Index i) const
{
    return {e[(i + 1) % N], e[(i + 2) % N], e[(i + 3) % N]};
}

bool Triangulation::on_third(const SurfaceElement& e) const
{
    return isnear(e[0]->first[2], dim_[2]) && isnear(e[1]->first[2], dim_[2])
        && isnear(e[2]->first[2], dim_[2]);
}

bool Triangulation::on_first(const NodePtr& n) const
{
    auto& p = n->first;
    // Sigma_1
    if (isnear(p[0], 0) || isnear(p[0], dim_[0])) {
        return true;
    }
    // Sigma_2
    if (isnear(p[1], 0) || isnear(p[1], dim_[1])) {
        return true;
    }
    return false;
}

Triangulation Triangulation::cuboid(const array<double, 3>& dim, size_t scale)
{
    Triangulation result(dim);

    array<size_t, 3> m;
    array<double, 3> s;
    auto d = *min_element(begin(dim), end(dim));

    scale = max(1ull, scale);
    for (size_t i = 0; i < 3; ++i) {
        m[i] = scale * size_t(dim[i] / d);
        s[i] = dim[i] / m[i];
    }

    result.elems_.reserve(6 * m[0] * m[1] * m[2]);
    for (size_t i = 0; i < m[0]; ++i) {
        auto x0 = i * s[0];
        auto x1 = (i + 1) * s[0];

        for (size_t j = 0; j < m[1]; ++j) {
            auto y0 = j * s[1];
            auto y1 = (j + 1) * s[1];

            for (size_t k = 0; k < m[2]; ++k) {
                auto z0 = k * s[2];
                auto z1 = (k + 1) * s[2];

                auto p = result.append_nodes(
                    combine({x0, x1}, {y0, y1}, {z0, z1}));
                result.append_elem({p[0], p[1], p[2], p[6]});
                result.append_elem({p[0], p[1], p[4], p[6]});
                result.append_elem({p[1], p[3], p[6], p[7]});
                result.append_elem({p[1], p[5], p[6], p[7]});
                result.append_elem({p[1], p[2], p[3], p[6]});
                result.append_elem({p[1], p[4], p[5], p[6]});
            }
        }
    }
    return result;
}
