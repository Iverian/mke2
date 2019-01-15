#include <debug.hpp>
#include <dense_matrix.hpp>
#include <sparce_matrix.hpp>
#include <util.hpp>
#include <vec.hpp>

#include <omp.h>

#include <algorithm>
#include <cmath>

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

Vec operator*(const Vec& lhs, const SparceMatrix& rhs)
{
    Vec result(rhs.m_, 0.);
    dot(result, lhs, rhs);
    return result;
}

void dot(Vec& result, const SparceMatrix& lhs, const Vec& rhs)
{
    check_if(lhs.m_ == rhs.size(), "Incompatible shapes");

    for (SparceMatrix::Index i = 0; i < lhs.m_; ++i) {
        auto u = lhs.diag_[i] * rhs[i];
        auto k = lhs.indptr_[i];
        auto last = lhs.indptr_[i + 1];
        for (; k < last; ++k) {
            u += lhs.data_[k] * rhs[lhs.indices_[k]];
        }
        result[i] = u;
    }
}

void dot(Vec& result, const Vec& lhs, const SparceMatrix& rhs)
{
    check_if(lhs.size() == rhs.m_, "Incompatible shapes");

    for (SparceMatrix::Index i = 0; i < rhs.m_; ++i) {
        auto k = rhs.indptr_[i];
        auto last = rhs.indptr_[i + 1];
        for (; k < last; ++k) {
            auto j = rhs.indices_[k];
            result[j] += lhs[j] * rhs.data_[k];
        }
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
