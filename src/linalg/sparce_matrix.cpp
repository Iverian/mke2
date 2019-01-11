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

SparceMatrix::SparceMatrix(const DataContainer& data,
                           const IndexContainer& indptr,
                           const IndexContainer& indices, Shape shape)
    : data_(data)
    , indptr_(indptr)
    , indices_(indices)
    , shape_(shape)
{
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

SparceMatrix::Value SparceMatrix::fetch_set(Index i, Index j, Value val)
{
    Value result = 0;
    if (!isnear(val, 0)) {
        result = fetch_modify(i, j, [&val](auto) { return val; });
    }
    return result;
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
    const Value* result = nullptr;

    if (index_in_range(i, j)) {
        auto b = begin(indices_);
        auto ifirst = b + indptr_[i];
        auto ilast = b + indptr_[i + 1];
        if (ifirst < ilast) {
            auto it = lower_bound(ifirst, ilast, j);
            if (it != ilast) {
                result = &data_[it - b];
            }
        }
    }

    return result;
}

bool SparceMatrix::index_in_range(Index i, Index j) const
{
    return i < shape_.m && j < shape_.n;
}

Vec operator*(const SparceMatrix& lhs, const Vec& rhs)
{
    check_if(lhs.shape_.n == rhs.size(), "Incompatible shapes");

    ptrdiff_t i, m = ptrdiff_t(lhs.shape_.m);
    Vec result(m, 0.);

#pragma omp parallel for
    for (i = 0; i < m; ++i) {
        double u = 0;
        auto ifirst = lhs.indptr_[i];
        auto ilast = lhs.indptr_[i + 1];
        if (ifirst < ilast) {
            for (size_t k = ifirst; k < ilast; ++k) {
                u += lhs.data_[k] * rhs[lhs.indices_[k]];
            }
        }
#pragma omp flush(result)
        result[i] = u;
    }

    return result;
}

ostream& operator<<(ostream& os, const SparceMatrix& obj)
{
    auto p = obj.indptr_.size() - 1;
    auto m = obj.shape_.m;
    SparceMatrix::Index i, j;

    os << "{ \"shape\": " << obj.shape_ << ", \"data\": [";
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

Vec& admul(Vec& result, const Vec& x, const Vec& y, double c)
{
    ptrdiff_t i, n = ptrdiff_t(result.size());

#pragma omp parallel for
    for (i = 0; i < n; ++i) {
        result[i] = x[i] + c * y[i];
    }

    return result;
}

Vec solve_cg(const SparceMatrix& lhs, const Vec& rhs, Vec x0)
{
    static constexpr size_t max_step = 10000;
    static constexpr double accuracy = 1e-5;

    size_t step = 0;
    size_t m = lhs.shape_.m;
    double eps = 0;

    check_if(rhs.size() == lhs.shape_.n, "Incompatible shapes");

    if (x0.empty()) {
        x0.resize(m, 0);
    }
    auto rsqr = sqr(rhs);

    auto r = rhs - lhs * x0;
    auto x = x0;
    auto z = r;

    double a, c, d;
    do {
        auto t = lhs * z;
        c = sqr(r);
        a = c / dot(t, z);

        admul(x, x, z, a);
        admul(r, r, t, -a);
        d = sqr(r);
        if (eps = sqrt(d / rsqr); eps < accuracy) {
            break;
        }

        admul(z, r, z, d / c);
    } while (++step, step < max_step);

    check_if(eps < accuracy,
             "Iteration didn't converge after %llu steps with eps=%.5g", step,
             eps);

    return x;
}
