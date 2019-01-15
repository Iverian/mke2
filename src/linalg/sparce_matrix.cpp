#include <cmath>
#include <debug.hpp>
#include <dense_matrix.hpp>
#include <sparce_matrix.hpp>
#include <util.hpp>
#include <vec.hpp>

#include <omp.h>

#include <algorithm>

using namespace std;

SparceMatrix::SparceMatrix()
    : data_()
    , diag_()
    , indptr_(1)
    , indices_()
    , m_(0)
{
}

#define _(i, j) ((i) * (m_) + (j))

SparceMatrix::SparceMatrix(Index side, const DataContainer& vals)
    : data_()
    , diag_(side)
    , indptr_(side + 1)
    , indices_()
    , nnz_(0)
    , m_(side)
{
    if (!vals.empty()) {
        for (Index i = 0; i < m_; ++i) {
            for (Index j = 0; j < m_; ++j) {
                auto& v = vals[_(i, j)];
                if (!isnear(v, 0)) {
                    if (i == j) {
                        diag_[i] = v;
                    } else {
                        data_.emplace_back(v);
                        indices_.emplace_back(j);
                        ++indptr_[i + 1];
                    }
                    ++nnz_;
                }
            }
        }

        data_.shrink_to_fit();
        indices_.shrink_to_fit();
        if (!data_.empty()) {
            for (Index i = 0; i < m_; ++i) {
                indptr_[i + 1] += indptr_[i];
            }
        }
    }
}

#undef _

SparceMatrix::Index SparceMatrix::size() const
{
    return m_ * m_;
}

SparceMatrix::Shape SparceMatrix::shape() const
{
    return {m_, m_};
}

SparceMatrix::Index SparceMatrix::non_zero() const
{
    return nnz_;
}

const SparceMatrix::DataContainer& SparceMatrix::data() const noexcept
{
    return data_;
}

const SparceMatrix::DataContainer& SparceMatrix::diag() const noexcept
{
    return diag_;
}

const SparceMatrix::IndexContainer& SparceMatrix::indptr() const noexcept
{
    return indptr_;
}

const SparceMatrix::IndexContainer& SparceMatrix::indices() const noexcept
{
    return indices_;
}

SparceMatrix::Value SparceMatrix::operator()(Index i, Index j) const
{
    check_if(i < m_ && j < m_, "Index out of range");

    Value result = 0;
    if (i == j) {
        result = diag_[i];
    } else {
        auto b = begin(indices_);
        auto ifirst = b + indptr_[i];
        auto ilast = b + indptr_[i + 1];
        if (ifirst < ilast) {
            auto p = lower_bound(ifirst, ilast, j);
            if (p != ilast && *p == j) {
                result = data_[Index(p - b)];
            }
        }
    }

    return result;
}

SparceMatrix::Value SparceMatrix::add(Index i, Index j, Value val)
{
    return fetch_modify(i, j, [&val](auto x) { return x + val; });
}

SparceMatrix::Value SparceMatrix::sub(Index i, Index j, Value val)
{
    return fetch_modify(i, j, [&val](auto x) { return x - val; });
}

SparceMatrix::Value SparceMatrix::set(Index i, Index j, Value val)
{
    return fetch_modify(i, j, [&val](auto) { return val; });
}

void SparceMatrix::remove_zeroes()
{
    DataContainer new_data;
    IndexContainer new_indptr(m_ + 1, 0);
    IndexContainer new_indices;

    for (Index i = 0; i < m_; ++i) {
        auto k = indptr_[i];
        auto last = indptr_[i + 1];
        for (; k < last; ++k) {
            auto v = data_[k];
            if (!isnear(v, 0)) {
                new_data.emplace_back(v);
                new_indices.emplace_back(indices_[k]);
                ++new_indptr[i + 1];
            } else {
                --nnz_;
            }
        }
    }

    new_data.shrink_to_fit();
    new_indices.shrink_to_fit();
    for (Index i = 0; i < m_; ++i) {
        new_indptr[i + 1] += new_indptr[i];
    }

    data_ = move(new_data);
    indptr_ = move(new_indptr);
    indices_ = move(new_indices);
}

Vec operator*(const SparceMatrix& lhs, const Vec& rhs)
{
    Vec result(lhs.m_, 0.);
    dot(result, lhs, rhs);
    return result;
}

void cax_y(Vec& res, double a, const Vec& x, const Vec& y);
void cax_by_z(Vec& res, double a, const Vec& x, double b, const Vec& y,
              const Vec& z);
double cdist(const Vec& lhs, const Vec& rhs);

Vec solve(const SparceMatrix& lhs, const Vec& rhs, Vec x0)
{
    static constexpr size_t max_iter = 100000;

    auto [m, n] = lhs.shape();
    check_if(rhs.size() == n, "Incompatible shapes");

    size_t step = 0;

    if (x0.empty()) {
        x0.resize(m, 0);
    }
    auto rsqr = sqr(rhs);
    auto r = rhs - lhs * x0;
    auto r0 = r;
    auto x = move(x0);

    // Vec p(m, 0.);
    Vec v(m, 0.);
    Vec s(m, 0.);
    Vec t(m, 0.);
    Vec y(m, 0.);
    Vec z(m, 0.);
    Vec q(m, 0.);
    Vec p(m, 0.);

    double res = 0, u = 1, a = 1, w = 1, b = 0;
    do {
        auto up = u;
        u = dot(r0, r);

        b = (u / up) * (a / w);
        // p = r + b * (p - w * v)
        cax_by_z(p, b, p, -w * b, v, r);
        // v = lhs * p
        dot(v, lhs, p);

        a = u / dot(r0, v);

        // s = r - a * v
        cax_y(s, -a, v, r);
        // t = lhs * s
        dot(t, lhs, s);
        w = dot(t, s) / sqr(t);

        // x = a * p + w * s + x
        cax_by_z(x, a, p, w, s, x);
        // r = s - w * t
        cax_y(r, -w, t, s);

        res = sqrt(sqr(r) / rsqr);
        if (isnear(res, 0, Tolerance::TRIPLE)) {
            break;
        }
    } while (++step < max_iter);

    auto diff = sqrt(sqr(lhs * x - rhs));
    if (step != max_iter) {
        cout << "Iteration converged";
    } else {
        cout << "Iteration did not converge";
    }
    cout << ": res=" << res << ", |Ax - b|=" << diff << ", step=" << step
         << endl;

    return x;
}

Vec solve_p(const SparceMatrix& lhs, const Vec& rhs, Vec x0)
{
    static constexpr size_t max_iter = 100000;

    auto [m, n] = lhs.shape();
    check_if(rhs.size() == n, "Incompatible shapes");
    check_if(!any_of(begin(lhs.diag_), end(lhs.diag_),
                     [](auto& x) { return isnear(x, 0); }),
             "Zero on diag");

    size_t step = 0;

    if (x0.empty()) {
        x0.resize(m, 0);
    }
    auto rsqr = sqr(rhs);
    auto r = rhs - lhs * x0;
    auto r0 = r;
    auto x = move(x0);

    // Vec p(m, 0.);
    Vec v(m, 0.);
    Vec s(m, 0.);
    Vec t(m, 0.);
    Vec y(m, 0.);
    Vec z(m, 0.);
    Vec q(m, 0.);

    mdot_diag(z, lhs, r);
    auto p = z;

    double res = 0, u = 1, a = 1, w = 1, b = 0;
    do {
        auto up = u;
        u = dot(r0, r);

        auto b = (u / up) * (a / w);
        // p = r + b * (p - w * v)
        cax_by_z(p, b, p, -w * b, v, r);
        // y = diag(lhs)^{-1} * p
        mdot_diag(y, lhs, p);
        // v = lhs * y
        dot(v, lhs, y);

        a = u / dot(r0, v);

        // s = r - a * v
        cax_y(s, -a, v, r);
        // z = diag(lhs)^{-1} * s
        mdot_diag(z, lhs, s);
        // t = lhs * s
        dot(t, lhs, z);
        // q = diag(lhs)^{-1} * t
        mdot_diag(q, lhs, t);
        w = dot(z, q) / sqr(q);

        // x = a * p + w * s + x
        cax_by_z(x, a, y, w, z, x);
        // r = s - w * t
        cax_y(r, -w, t, s);

        res = sqrt(sqr(r) / rsqr);
        if (isnear(res, 0, Tolerance::DOUBLE)) {
            break;
        }
    } while (++step < max_iter);

    res = sqrt(sqr(lhs * x - rhs));
    if (isnear(res, 0, Tolerance::SINGLE)) {
        cout << "Iteration converged: res = " << res << ", step = " << step
             << endl;
    } else {
        cout << "Iteration did not converge: res = " << res
             << ", step = " << step << endl;
    }

    return x;
}

void mdot_diag(Vec& result, const SparceMatrix& lhs, const Vec& rhs)
{
    for (auto i = 0; i < lhs.m_; ++i) {
        result[i] = rhs[i] / lhs.diag_[i];
    }
}

void dot(Vec& result, const SparceMatrix& lhs, const Vec& rhs)
{
    check_if(lhs.m_ == rhs.size(), "Incompatible shapes");

    SparceMatrix::Index i, m = lhs.m_;
    size_t ifirst, ilast, k;
    double u;

    for (i = 0; i < m; ++i) {
        u = lhs.diag_[i] * rhs[i];
        ifirst = lhs.indptr_[i];
        ilast = lhs.indptr_[i + 1];
        if (ifirst < ilast) {
            for (k = ifirst; k < ilast; ++k) {
                u += lhs.data_[k] * rhs[lhs.indices_[k]];
            }
        }
        result[i] = u;
    }
}

void cax_y(Vec& res, double a, const Vec& x, const Vec& y)
{
    SparceMatrix::Index i, n = res.size();

    for (i = 0; i < n; ++i) {
        res[i] = a * x[i] + y[i];
    }
}

void cax_by_z(Vec& res, double a, const Vec& x, double b, const Vec& y,
              const Vec& z)
{
    SparceMatrix::Index i, n = res.size();

    for (i = 0; i < n; ++i) {
        res[i] = a * x[i] + b * y[i] + z[i];
    }
}

ostream& operator<<(ostream& os, const SparceMatrix& obj)
{
    auto m = obj.m_;
    SparceMatrix::Index i, j;

    os << "{\"shape\": " << obj.shape() << ", \"data\": [";
    for (i = 0; i < m; ++i) {
        auto pos = obj.indptr_[i];
        auto end = obj.indptr_[i + 1];

        os << "[";
        for (j = 0; j < m; ++j) {
            if (i == j) {
                os << obj.diag_[i];
            }
            if (pos < end && obj.indices_[pos] == j) {
                os << obj.data_[pos++];
            } else {
                os << 0;
            }
            os << (j + 1 != m ? ", " : "]");
        }
        os << (i + 1 != m ? ", " : "]");
    }
    os << "}";

    return os;
}
