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

pair<SparceMatrix::Value*, bool> SparceMatrix::insert(Index i, Index j,
                                                      Value val)
{
    Value* result = nullptr;
    bool inserted = false;

    if (index_in_range(i, j)) {
        auto b = begin(indices_);
        auto ifirst = b + indptr_[i];
        auto ilast = b + indptr_[i + 1];

        auto it = lower_bound(ifirst, ilast, j);
        if (*it == j) {
            result = &data_[it - b];
            *result = val;
        } else {
            auto new_index = indices_.emplace(it, j);
            result = &(*data_.emplace(begin(data_) + (it - b), val));
            inserted = true;

            for (auto j = i + 1; j < shape_.m + 1; ++j) {
                ++indptr_[j];
            }
        }
    }

    return make_pair(result, inserted);
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
    auto& m = lhs.shape_.m;

    SparceMatrix::Index i, k;
    Vec result(m, 0.);

#pragma omp parallel for
    for (i = 0; i < m; ++i) {
        double u = 0;
        auto ifirst = lhs.indptr_[i];
        auto ilast = lhs.indptr_[i + 1];
        if (ifirst < ilast) {
            for (k = ifirst; k < ilast; ++k) {
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
