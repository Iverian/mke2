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

SparceMatrix::Value SparceMatrix::operator()(Index i, Index j) const
{
    auto ptr = find(i, j);
    return (ptr != nullptr) ? *ptr : 0;
}

SparceMatrix::Value SparceMatrix::fetch_add(Index i, Index j, Value val)
{
    return fetch_modify(i, j, [&val](auto x) { return x + val; });
}

SparceMatrix::Value SparceMatrix::fetch_sub(Index i, Index j, Value val)
{
    return fetch_modify(i, j, [&val](auto x) { return x - val; });
}

SparceMatrix::Value SparceMatrix::fetch_set(Index i, Index j, Value val)
{
    return fetch_modify(i, j, [&val](auto) { return val; });
}

void SparceMatrix::clean_up()
{
    DataContainer new_data;
    IndexContainer new_indptr(m_ + 1, 0);
    IndexContainer new_indices;

    Index pos = 0;
    for (Index i = 0; i < m_; ++i) {
        for (; pos < indptr_[i + 1]; ++pos) {
            auto& v = data_[pos];
            if (!isnear(v, 0)) {
                new_data.emplace_back(v);
                new_indices.emplace_back(indices_[pos]);
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

SparceMatrix::Value* SparceMatrix::find(Index i, Index j)
{
    auto result = const_cast<const SparceMatrix*>(this)->find(i, j);
    return const_cast<Value*>(result);
}

const SparceMatrix::Value* SparceMatrix::find(Index i, Index j) const
{
    check_if(i < m_ && j < m_, "Index out of range");

    const Value* result = nullptr;
    if (i == j) {
        result = &diag_[i];
    } else {
        auto b = begin(indices_);
        auto ifirst = b + indptr_[i];
        auto ilast = b + indptr_[i + 1];
        if (ifirst < ilast) {
            auto p = lower_bound(ifirst, ilast, j);
            if (p != end(indices_) && *p == j) {
                result = &data_[Index(p - b)];
            }
        }
    }

    return result;
}

Vec operator*(const SparceMatrix& lhs, const Vec& rhs)
{
    Vec result(lhs.m_, 0.);
    dot(result.data(), lhs, rhs);
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

    double dst = 0, res = 0, u = 1, a = 1, w = 1, b = 0;
    do {
        dot(v.data(), lhs, p);
        w = dot(r, z);
        a = w / dot(p, v);
        cax_y(x, a, p, x);
        copy(begin(r), end(r), begin(q));
        cax_y(r, -a, v, r);
        res = sqrt(sqr(r) / rsqr);
        if (isnear(res, 0, Tolerance::DOUBLE)) {
            break;
        }
        mdot_diag(z, lhs, r);
        b = dot(r, z) / w;
        cax_y(p, b, p, z);
        // auto up = u;
        // u = dot(r0, r);

        // auto b = (u / up) * (a / w);
        // // p = r + b * (p - w * v)
        // cax_by_z(p, b, p, -w * b, v, r);
        // mdot_diag(y, lhs, p);
        // // v = lhs * p
        // dot(v.data(), lhs, y);

        // a = u / dot(r0, v);

        // // s = r - a * v
        // cax_y(s, -a, v, r);
        // mdot_diag(z, lhs, s);
        // // t = lhs * s
        // dot(t.data(), lhs, z);

        // w = dot(t, s) / sqr(t);
        // // x = a * p + w * s + x
        // cax_by_z(x, a, p, w, s, x);
        // // r = s - w * t
        // cax_y(r, -w, t, s);
        // res = sqrt(sqr(r) / rsqr);
        // if (isnear(res, 0, Tolerance::DOUBLE)) {
        //     break;
        // }
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

void dot(double* const result, const SparceMatrix& lhs, const Vec& rhs)
{
    check_if(lhs.m_ == rhs.size(), "Incompatible shapes");

    ptrdiff_t i;
    const ptrdiff_t m = ptrdiff_t(lhs.m_);
    size_t ifirst, ilast, k;
    double u;

#pragma omp parallel shared(result, lhs, rhs, m)
    {
#pragma omp for private(i, k, u, ifirst, ilast)
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
}

void cax_y(Vec& res, double a, const Vec& x, const Vec& y)
{
    int i;
    const int n = int(res.size());

#pragma omp parallel shared(res, a, x, y, n)
    {
#pragma omp for private(i)
        for (i = 0; i < n; ++i) {
            res[i] = a * x[i] + y[i];
        }
    }
}

void cax_by_z(Vec& res, double a, const Vec& x, double b, const Vec& y,
              const Vec& z)
{
    int i;
    const int n = int(res.size());

#pragma omp parallel shared(res, a, x, b, y, z, n)
    {
#pragma omp for private(i)
        for (i = 0; i < n; ++i) {
            res[i] = a * x[i] + b * y[i] + z[i];
        }
    }
}

// ostream& operator<<(ostream& os, const SparceMatrix& obj)
// {
//     auto p = obj.indptr_.size() - 1;
//     auto m = obj.shape_.m;
//     SparceMatrix::Index i, j;

//     os << "{\"shape\": " << obj.shape_ << ", \"data\": [";
//     for (i = 0; i < p; ++i) {
//         auto pos = obj.indptr_[i];
//         auto end = obj.indptr_[i + 1];

//         os << "[";
//         for (j = 0; j < m; ++j) {
//             if (pos < end && obj.indices_[pos] == j) {
//                 os << obj.data_[pos++];
//             } else {
//                 os << 0;
//             }
//             os << (j + 1 != m ? ", " : "]");
//         }
//         os << (i + 1 != p ? ", " : "]");
//     }
//     os << "}";

//     return os;
// }
