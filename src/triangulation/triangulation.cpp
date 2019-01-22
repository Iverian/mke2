#include <debug.hpp>
#include <global_indices.hpp>
#include <triangulation.hpp>
#include <util.hpp>

#include <algorithm>
#include <fstream>
#include <iterator>
#include <string>
#include <unordered_set>

using namespace std;

Triangulation::Triangulation()
    : dim_ {0, 0, 0}
    , size_(0)
    , nodes_()
    , elems_()
    , first_()
    , third_()
{
}

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

const vector<pair<Triangulation::Index, Triangulation::SurfaceElement>>&
Triangulation::triangles() const
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
        for (Index i = 0; i < 2; ++i) {
            if (first_flag[i]) {
                first_[i].push_back(result);
            }
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
    return SurfaceElement(
        {(*this)[(i + 1) % N], (*this)[(i + 2) % N], (*this)[(i + 3) % N]});
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
    OnFirst result {{false, false}};

    auto& p = n.point();

    for (Index i = 0; i < 2; ++i) {
        if (iszero(p[i]) || isnear(p[i], dim_[i])) {
            result[i] = true;
        }
    }

    return result;
}

Triangulation::OnFirst Triangulation::on_first(OnFirst cond, Index coord) const
{
    return {{cond[0] && (coord == Index(Coord::X) || coord == Index(Coord::Z)),
             cond[1] && (coord == Index(Coord::Y))}};
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
    auto p = elems_.size();
    for (Index i = 0; i < p; ++i) {
        auto& k = elems_[i];
        for (Index j = 0; j < N; ++j) {
            auto f = k.face(j);
            if (is_boundary(f)) {
                triangles_.emplace_back(i, move(f));
            }
        }
    }
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
                /*
                 *  i xyz P
                 *  0 000 A
                 *  1 001 A1
                 *  2 010 D
                 *  3 011 D1
                 *  4 100 B
                 *  5 101 B1
                 *  6 110 C
                 *  7 111 C1
                 */
                auto p = result.append_nodes(
                    combine({x0, x1}, {y0, y1}, {z0, z1}));

                result.append_elem(p[0], p[1], p[4], p[6]);
                result.append_elem(p[0], p[1], p[2], p[6]);
                result.append_elem(p[1], p[4], p[5], p[6]);
                result.append_elem(p[1], p[5], p[6], p[7]);
                result.append_elem(p[1], p[3], p[6], p[7]);
                result.append_elem(p[1], p[2], p[3], p[6]);
            }
        }
    }

    return result;
}

ostream& operator<<(ostream& os, const Triangulation& obj)
{
    size_t i, j;
    const auto m = obj.nodes_.size();
    const auto n = obj.elems_.size();
    const auto q = obj.third_.size();

    vector<Point3d> nodes(m);
    for (auto& n : obj.nodes_) {
        nodes[n.second] = n.first;
    }
    os << "{";
    os << " \"dim\": [ " << obj.dim_[0] << ", " << obj.dim_[1] << ", "
       << obj.dim_[2] << "],";

    os << " \"nodes\": [ ";
    for (i = 0; i < m; ++i) {
        os << nodes[i] << (i + 1 != m ? ", " : "],");
    }

    os << " \"elems\":[ ";
    for (i = 0; i < n; ++i) {
        os << obj.elems_[i] << (i + 1 != n ? ", " : "],");
    }

    os << " \"first\":[ ";
    for (i = 0; i < 2; ++i) {
        const auto p = obj.first_[i].size();
        os << "[";
        for (j = 0; j < p; ++j) {
            os << obj.first_[i][j].index() << (j + 1 != p ? ", " : "]");
        }
        os << (i + 1 != 2 ? ", " : "],");
    }

    os << " \"third\":[ ";
    for (i = 0; i < q; ++i) {
        os << obj.third_[i] << (i + 1 != q ? ", " : "]");
    }
    os << "}";

    return os;
}

Triangulation Triangulation::from_msh(const char* filename,
                                      array<double, DIM> dim)
{
    ifstream is(filename);
    Triangulation result(dim);

    string buf;
    // skip useless strings
    for (int i = 0; i < 4; ++i) {
        getline(is, buf);
    }
    check_ifd(buf == "$Nodes", "Parsing error");

    Index node_count = 0;
    is >> node_count;

    vector<NodePtr> nodes(node_count);

    for (Index i = 0; i < node_count; ++i) {
        Index ind = 0;
        double x = 0., y = 0., z = 0.;

        is >> ind >> x >> y >> z;
        auto p = result.append_node({x, y, z});
        nodes[i] = p;
    }

    // skip useless strings
    for (int i = 0; i < 3; ++i) {
        getline(is, buf);
    }
    check_ifd(buf == "$Elements", "Parsing error");

    Index elem_count = 0;
    is >> elem_count;

    for (Index i = 0; i < elem_count; ++i) {
        Index ind = 0, ele_type = 0;
        is >> ind >> ele_type;
        if (ele_type == 4) {
            Index ntags, a, b, c, d = 0;
            is >> ntags;
            for (Index j = 0; j < ntags; ++j) {
                Index tag = 0;
                is >> tag;
            }

            is >> a >> b >> c >> d;
            result.append_elem(nodes[a - 1], nodes[b - 1], nodes[c - 1],
                               nodes[d - 1]);
        } else {
            getline(is, buf);
        }
    }

    return result;
}

ostream& operator<<(ostream& os, const Triangulation::SurfaceElement& obj)
{
    return os << "[" << obj[0].index() << ", " << obj[1].index() << ", "
              << obj[2].index() << "]";
}

ostream& operator<<(ostream& os, const Triangulation::FiniteElement& obj)
{
    return os << "[" << obj[0].index() << ", " << obj[1].index() << ", "
              << obj[2].index() << ", " << obj[3].index() << "]";
}