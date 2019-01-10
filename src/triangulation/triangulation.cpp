#include <debug.hpp>
#include <point3d.hpp>
#include <triangulation.hpp>
#include <util.hpp>

#include <algorithm>
#include <cmath>

using namespace std;

void Triang3d::append(const array<Point3d, N> vertices)
{
    nodes_.push_back(make_shared<Node>(*this));
    auto result = nodes_.back();
    auto new_index = nodes_.size() - 1;

    for (Index i = 0; i < N; ++i) {
        auto j = vertices_.insert({vertices[i], {}}).first;
        result->vindices_[i] = &(*j);
        j->second.emplace_back(new_index);
    }
}

#define _adj(i, k) (cur.vindices_[((i) + (k)) % N]->second)

void Triang3d::set_adjacent()
{
    auto s = nodes_.size();
    for (Index i = 0; i < s; ++i) {
        auto& cur = *nodes_[i];
        for (Index j = 0; j < N; ++j) {
            cur.adjacent_[j]
                = find_adjacent(_adj(j, 1), _adj(j, 2), _adj(j, 3), i);
        }
    }
}

#undef _adj

Triang3d::Index Triang3d::find_adjacent(const IndexSet& a, const IndexSet& b,
                                        const IndexSet& c, Index cur)
{
    IndexSet tmp;
    IndexSet result;

    set_intersection(begin(a), end(a), begin(b), end(b), back_inserter(tmp));
    set_intersection(begin(c), end(c), begin(tmp), end(tmp),
                     back_inserter(result));

    auto s = result.size();
    check_if(s == 1 || s == 2, "Некорректное число соседей");
    if (s == 1) {
        return npos;
    } else {
        return (result.front() == cur) ? result.back() : result.front();
    }
}

void Triang3d::reserve(Index node_count)
{
    nodes_.reserve(node_count);
}

void Triang3d::shrink_to_fit()
{
    nodes_.shrink_to_fit();
}

const Triang3d::NodeContainer& Triang3d::nodes() const
{
    return nodes_;
}

Triang3d::Triang3d::Node::Node(const Triang3d& parent)
    : parent_(&parent)
    , vindices_()
    , adjacent_{npos, npos, npos, npos}
{
}

const Point3d& Triang3d::Node::v(Index i) const
{
    return vindices_[i]->first;
}

Triang3d::NodePtr Triang3d::Node::adj(Index i) const
{
    auto index = iadj(i);
    return index != npos ? parent_->nodes_[index] : nullptr;
}

Triang3d::Index Triang3d::Node::iadj(Index i) const
{
    return adjacent_[i];
}

vector<Triang3d::Index> Triang3d::Node::boundary() const
{
    vector<Index> result;
    result.reserve(N);
    for (Index i = 0; i < N; ++i) {
        if (adjacent_[i] == npos) {
            result.push_back(i);
        }
    }
    result.shrink_to_fit();
    return result;
}

array<Point3d, Triang3d::N - 1> Triang3d::Node::face(Index i) const
{
    return {v((i + 1) % N), v((i + 2) % N), v((i + 3) % N)};
}

Point3d Triang3d::Node::normal(Index i) const
{
    auto p = v((i + 1) % N);
    auto m = cross(v((i + 2) % N) - p, v((i + 3) % N) - p);
    if (dot(m, v(i) - p) > 0) {
        m = -m;
    }
    return m;
}

DenseMatrix Triang3d::Node::vertex_view() const
{
    DenseMatrix result({N, 3});
    for (size_t i = 0; i < N; ++i) {
        auto& cur = v(i);
        for (size_t j = 0; j < 3; ++j) {
            result(i, j) = cur[j];
        }
    }
    return result;
}

Triang3d triangulate_cuboid(array<double, 3> dim, size_t scale)
{
    Triang3d result;

    array<size_t, 3> m;
    array<double, 3> s;
    auto d = *min_element(begin(dim), end(dim));

    scale = max(1ull, scale);
    for (size_t i = 0; i < 3; ++i) {
        m[i] = scale * size_t(dim[i] / d);
        s[i] = dim[i] / m[i];
    }

    result.reserve(6 * m[0] * m[1] * m[2]);
    for (size_t i = 0; i < m[0]; ++i) {
        auto x0 = i * s[0];
        auto x1 = (i + 1) * s[0];
        for (size_t j = 0; j < m[1]; ++j) {
            auto y0 = j * s[1];
            auto y1 = (j + 1) * s[1];
            for (size_t k = 0; k < m[2]; ++k) {
                auto z0 = k * s[2];
                auto z1 = (k + 1) * s[2];

                auto p = combine({x0, x1}, {y0, y1}, {z0, z1});
                result.append({p[0], p[1], p[2], p[6]});
                result.append({p[0], p[1], p[4], p[6]});
                result.append({p[1], p[3], p[6], p[7]});
                result.append({p[1], p[5], p[6], p[7]});
                result.append({p[1], p[2], p[3], p[6]});
                result.append({p[1], p[4], p[5], p[6]});
            }
        }
    }
    result.shrink_to_fit();
    result.set_adjacent();

    return result;
}

ostream& operator<<(ostream& os, const Triang3d& obj)
{
    os << "[";
    Triang3d::Index i, j, s = obj.nodes_.size();
    for (i = 0; i < s; ++i) {
        auto& node = *obj.nodes_[i];
        os << "{ \"index\": " << i << ", \"vertices\": [";
        for (j = 0; j < Triang3d::N; ++j) {
            os << node.v(j) << (j + 1 != Triang3d::N ? ", " : "], ");
        }
        os << "\"adjacent\": [";
        for (j = 0; j < Triang3d::N; ++j) {
            auto k = node.iadj(j);
            if (k == Triang3d::npos) {
                os << "-1";
            } else {
                os << k;
            }
            os << (j + 1 != Triang3d::N ? ", " : "]");
        }
        os << "}" << (i + 1 != s ? ", " : "]");
    }
    return os;
}