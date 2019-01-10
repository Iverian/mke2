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
    Value result = 0;

    if (index_in_range(i, j)) {
        auto b = begin(indices_);
        auto ifirst = b + indptr_[i];
        auto ilast = b + indptr_[i + 1];
        auto pos = lower_bound(ifirst, ilast, j);

        if (pos != end(indices_) && *pos == j) {
            result = (data_[pos - b] += val);
        } else if (!isnear(val, 0)) {
            indices_.emplace(pos, j);
            data_.emplace(begin(data_) + (pos - b), val);
            result = val;

            for (auto k = i + 1; k < shape_.m; ++k) {
                ++indptr_[k];
            }
        }
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

Vec CG(const SparseMatrix& mat, const Vec& right_side)
{
    int max_step = 10000;
    int step = 0;
    double accuracy = 1e-5;

    Vec x0(right_side.size(), 0), 
        r0(right_side), 
        z0(right_side), 
        xK(right_side.size()), 
        rK(right_side.size()), 
        zk(right_side.size()),
        temp(right_side.size());
    double alpha, beta;

    do {

        temp = mat * z0;

        alpha = (r0.dot(r0)) / (z0.dot(temp));

        xK = x0 + (z0 * alpha);

        rK = r0 - (temp * alpha);

        if (((rK.dot(rK)) / (right_side.dot(right_side))) < accuracy) 
            break;

        beta = (rK.dot(rK)) / (r0.dot(r0));

        zK = rK + z0 * beta;

        step++;

    } while(step < max_step);

    return xK;
}