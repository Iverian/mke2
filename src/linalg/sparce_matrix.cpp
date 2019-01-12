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
    , indptr_(1)
    , indices_()
    , shape_{0, 0}
{
}

SparceMatrix::SparceMatrix(const Shape& shape)
    : data_()
    , indptr_()
    , indices_()
    , shape_(shape)
{
    indptr_.resize(shape.m + 1);
}

SparceMatrix::SparceMatrix(const DenseMatrix& mat)
    : SparceMatrix(mat.shape())
{
    for (Index i = 0; i < shape_.m; ++i) {
        for (Index j = 0; j < shape_.n; ++j) {
            auto v = mat(i, j);
            if (!isnear(v, 0)) {
                data_.emplace_back(v);
                indices_.emplace_back(j);
                ++indptr_[i + 1];
            }
        }
    }

    data_.shrink_to_fit();
    indices_.shrink_to_fit();
    if (!data_.empty()) {
        for (auto it = next(begin(indptr_)); it != end(indptr_); ++it) {
            *it += *prev(it);
        }
    }
}

SparceMatrix::Index SparceMatrix::size() const
{
    return shape_.m * shape_.n;
}

SparceMatrix::Shape SparceMatrix::shape() const
{
    return shape_;
}

SparceMatrix::Index SparceMatrix::non_zero() const
{
    return data_.size();
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
    IndexContainer new_indptr(shape_.m + 1, 0);
    IndexContainer new_indices;

    Index pos = 0;
    for (Index i = 0; i < shape_.m; ++i) {
        for (; pos < indptr_[i + 1]; ++pos) {
            auto& v = data_[pos];
            if (!isnear(v, 0)) {
                new_data.emplace_back(v);
                new_indices.emplace_back(indices_[pos]);
                ++new_indptr[i + 1];
            }
        }
    }

    new_data.shrink_to_fit();
    new_indices.shrink_to_fit();
    for (Index i = 0; i < shape_.m; ++i) {
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
    check_if(i < shape_.m && j < shape_.n, "Index out of range");

    const Value* result = nullptr;

    auto b = begin(indices_);
    auto ifirst = b + indptr_[i];
    auto ilast = b + indptr_[i + 1];
    if (ifirst < ilast) {
        auto p = lower_bound(ifirst, ilast, j);
        if (p != end(indices_) && *p == j) {
            result = &data_[Index(p - b)];
        }
    }

    return result;
}

void dot(double* const result, const SparceMatrix& lhs, const Vec& rhs)
{
    check_if(lhs.shape_.n == rhs.size(), "Incompatible shapes");

    ptrdiff_t i;
    const ptrdiff_t m = ptrdiff_t(lhs.shape_.m);
    size_t ifirst, ilast, k;
    double u;

#pragma omp parallel shared(result, lhs, rhs)
    {
#pragma omp for private(i, k, u, ifirst, ilast)
        for (i = 0; i < m; ++i) {
            u = 0;
            ifirst = lhs.indptr_[i];
            ilast = lhs.indptr_[i + 1];
            if (ifirst < ilast) {
                for (k = ifirst; k < ilast; ++k) {
                    u += lhs.data_[k] * rhs[lhs.indices_[k]];
                }
            }
#pragma omp flush(result)
            result[i] = u;
        }
    }
}

Vec operator*(const SparceMatrix& lhs, const Vec& rhs)
{
    Vec result(lhs.shape_.m, 0.);
    dot(result.data(), lhs, rhs);
    return result;
}

ostream& operator<<(ostream& os, const SparceMatrix& obj)
{
    auto p = obj.indptr_.size() - 1;
    auto m = obj.shape_.m;
    SparceMatrix::Index i, j;

    os << "{\"shape\": " << obj.shape_ << ", \"data\": [";
    for (i = 0; i < p; ++i) {
        auto pos = obj.indptr_[i];
        auto end = obj.indptr_[i + 1];

        os << "[";
        for (j = 0; j < m; ++j) {
            if (pos < end && obj.indices_[pos] == j) {
                os << obj.data_[pos++];
            } else {
                os << 0;
            }
            os << (j + 1 != m ? ", " : "]");
        }
        os << (i + 1 != p ? ", " : "]");
    }
    os << "}";

    return os;
}

Vec& admul_0(Vec& result, const Vec& x, const Vec& y, double c)
{
    ptrdiff_t i;
    const ptrdiff_t n = ptrdiff_t(result.size());

#pragma omp parallel shared(result, x, y, c)
    {
#pragma omp for private(i)
        for (i = 0; i < n; ++i) {
            result[i] = x[i] + c * y[i];
        }
    }
    return result;
}

Vec& admul_1(Vec& result, const Vec& x, const Vec& y, const Vec& z, double c,
             double d)
{
    ptrdiff_t i;
    const ptrdiff_t n = ptrdiff_t(result.size());

#pragma omp parallel shared(result, x, y, z, c, d)
    {
#pragma omp for private(i)
        for (i = 0; i < n; ++i) {
            result[i] = x[i] + c * (y[i] + d * z[i]);
        }
    }

    return result;
}

double dist(const Vec& lhs, const Vec& rhs)
{
    double result = 0;
    ptrdiff_t i;
    const ptrdiff_t n = ptrdiff_t(lhs.size());

#pragma omp parallel shared(lhs, rhs, result)
    {
        double u = 0;
#pragma omp for private(i)
        for (i = 0; i < n; ++i) {
            u += sqr(lhs[i] - rhs[i]);
        }
#pragma omp flush(result)
#pragma omp critical
        result += u;
    }

    return sqrt(result);
}

Vec solve_cg(const SparceMatrix& lhs, const Vec& rhs, Vec x0, const double tol)
{
    static constexpr size_t max_step = 10000;

    check_if(rhs.size() == lhs.shape_.n, "Incompatible shapes");

    size_t step = 0;
    size_t m = lhs.shape_.m;
    double eps = 0;

    if (x0.empty()) {
        x0.resize(m, 0);
    }
    auto rsqr = sqr(rhs);

    auto r = rhs - lhs * x0;
    auto z = r;
    auto x = move(x0);

    Vec p(m, 0.);
    Vec v(m, 0.);
    Vec s(m, 0.);
    Vec t(m, 0.);
    Vec h(m, 0.);

    double up, b, uc = 1, a = 1, w = 1;
    do {
        up = uc;
        uc = dot(z, r);

        b = (uc / up) * (a / w);
        // p = r + b * (p - w * v)
        admul_1(p, r, p, v, b, -w);

        // v = lhs * p
        dot(v.data(), lhs, p);
        a = uc / dot(z, v);

        // h = x + a * p
        admul_0(h, x, p, a);
        // if (dist(x, h) < tol) {
        //     cdbg << "Exited at ??" << endl;
        //     return h;
        // }

        // s = r - a * v
        admul_0(s, r, v, -a);

        // t = lhs * s
        dot(t.data(), lhs, s);
        w = dot(t, s) / sqr(t);
        admul_0(x, h, s, w);

        // r = s - w * t
        admul_0(r, s, t, -w);

        if (eps = sqrt(sqr(r) / rsqr); eps < tol) {
            break;
        }
    } while (++step, step < max_step);

    check_if(eps < tol,
             "Iteration didn't converge after %llu steps with eps = %.5g",
             step, eps);

    coutd << "Iteration converged with step = " << step << " and eps = " << eps
          << endl;

    return x;
}
