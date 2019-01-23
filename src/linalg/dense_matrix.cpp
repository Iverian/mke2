#include <debug.hpp>
#include <dense_matrix.hpp>
#include <util.hpp>

#include <utility>

using namespace std;

DenseMatrix::DenseMatrix()
    : Super()
    , shape_{0, 0}
{
}

DenseMatrix::DenseMatrix(Index2d shape, const Super& data)
    : Super(data)
    , shape_(shape)
{
    if (data.empty()) {
        Super::resize(size());
    }
}

DenseMatrix::DenseMatrix(Index2d shape, Value val)
    : Super()
    , shape_(shape)
{
    resize(size(), val);
}

Index DenseMatrix::size() const
{
    return shape_.first * shape_.second;
}

Index2d DenseMatrix::shape() const
{
    return shape_;
}

DenseMatrix::Super& DenseMatrix::data()
{
    Super& result = *this;
    return result;
}

const DenseMatrix::Super& DenseMatrix::data() const
{
    const Super& result = *this;
    return result;
}

#define _(i, j) ((*this)[(shape_.second) * (i) + (j)])

Value& DenseMatrix::operator()(Index i, Index j)
{
    check_if(i < shape_.first && j < shape_.second, "Index out of range");

    return _(i, j);
}

const Value& DenseMatrix::operator()(Index i, Index j) const

{
    check_if(i < shape_.first && j < shape_.second, "Index out of range");

    return _(i, j);
}

DenseMatrix& DenseMatrix::swap_rows(Index i, Index j)
{
    check_if(i < shape_.first && j < shape_.first, "Index out of range");

    for (Index k = 0; k < shape_.second; ++k) {
        ::swap(_(i, k), _(j, k));
    }
    return *this;
}

DenseMatrix& DenseMatrix::swap_columns(Index i, Index j)
{
    check_if(i < shape_.second && j < shape_.second, "Index out of range");

    for (Index k = 0; k < shape_.first; ++k) {
        ::swap(_(k, i), _(k, j));
    }
    return *this;
}

DenseMatrix& DenseMatrix::append_row(const Vec& v)
{
    check_if(v.size() == shape_.second, "Incompatible shapes");

    resize(size() + shape_.second);
    for (Index j = 0; j < shape_.second; ++j) {
        (*this)(shape_.first, j) = v[j];
    }
    ++shape_.first;

    return *this;
}

DenseMatrix& DenseMatrix::append_column(const Vec& v)
{
    check_if(v.size() == shape_.first, "Incompatible shapes");

    reserve(size() + shape_.first);
    auto pos = begin() + shape_.second;
    for (Index i = 0; i < shape_.first; ++i) {
        insert(pos, v[i]);
        pos += shape_.second;
    }
    return *this;
}

DenseMatrix DenseMatrix::transpose() const
{
    DenseMatrix result({shape_.second, shape_.first});
    for (Index i = 0; i < shape_.first; ++i) {
        for (Index j = 0; j < shape_.second; ++j) {
            result(j, i) = (*this)(i, j);
        }
    }
    return result;
}

DenseMatrix& DenseMatrix::operator+=(const DenseMatrix& rhs)
{
    check_if(shape_ == rhs.shape_, "Incompatible shapes");

    for (Index i = 0; i < size(); ++i) {
        (*this)[i] += rhs[i];
    }
    return *this;
}

DenseMatrix& DenseMatrix::operator-=(const DenseMatrix& rhs)
{
    check_if(shape_ == rhs.shape_, "Incompatible shapes");

    for (Index i = 0; i < size(); ++i) {
        (*this)[i] -= rhs[i];
    }
    return *this;
}

DenseMatrix& DenseMatrix::operator*=(Value rhs)
{
    for (Index i = 0; i < size(); ++i) {
        (*this)[i] *= rhs;
    }
    return *this;
}

DenseMatrix& DenseMatrix::operator/=(Value rhs)
{
    for (Index i = 0; i < size(); ++i) {
        (*this)[i] /= rhs;
    }
    return *this;
}

#undef _

bool operator==(const DenseMatrix& lhs, const DenseMatrix& rhs)
{
    auto result = lhs.shape_ == rhs.shape_;
    if (result) {
        Value max = 0;
        auto s = lhs.size();
        for (Index i = 0; i < s; ++i) {
            if (auto cur = fabs(rhs[i] - lhs[i]); cur > max) {
                max = cur;
            }
        }
        result = iszero(max, Tolerance::DOUBLE);
    }
    return result;
}

bool operator!=(const DenseMatrix& lhs, const DenseMatrix& rhs)
{
    return !(lhs == rhs);
}

ostream& operator<<(ostream& os, const DenseMatrix& obj)
{
    os << "{\"shape\": " << obj.shape_ << ",\"data\": [";
    for (Index i = 0; i < obj.shape_.first; ++i) {
        os << "[";
        for (Index j = 0; j < obj.shape_.second; ++j) {
            os << obj(i, j) << (j + 1 != obj.shape_.second ? ", " : "]");
        }
        os << (i + 1 != obj.shape_.first ? ", " : "]");
    }
    return os << "}";
}

DenseMatrix operator+(const DenseMatrix& lhs, const DenseMatrix& rhs)
{
    auto result = lhs;
    return (result += rhs);
}

DenseMatrix operator-(const DenseMatrix& lhs, const DenseMatrix& rhs)
{
    auto result = lhs;
    return (result -= rhs);
}

DenseMatrix operator*(const DenseMatrix& lhs, const DenseMatrix& rhs)
{
    check_if(lhs.shape_.second == rhs.shape_.first, "Incompatible shapes");

    DenseMatrix result({lhs.shape_.first, rhs.shape_.second});
    Index i, j, k;
    for (i = 0; i < lhs.shape_.first; ++i) {
        for (k = 0; k < lhs.shape_.second; ++k) {
            for (j = 0; j < rhs.shape_.second; ++j) {
                result(i, j) += lhs(i, k) * rhs(k, j);
            }
        }
    }
    return result;
}

DenseMatrix operator*(const DenseMatrix& lhs, Value rhs)
{
    auto result = lhs;
    return (result *= rhs);
}

DenseMatrix operator*(Value lhs, const DenseMatrix& rhs)
{
    auto result = rhs;
    return (result *= lhs);
}

DenseMatrix operator/(const DenseMatrix& lhs, Value rhs)
{
    auto result = lhs;
    return (result /= rhs);
}

Vec operator*(const DenseMatrix& lhs, const Vec& rhs)
{
    check_if(lhs.shape_.second == rhs.size(), "Incompatible shapes");

    Vec result(lhs.shape_.first, 0.);
    for (Index i = 0; i < lhs.shape_.first; ++i) {
        for (Index j = 0; j < lhs.shape_.second; ++j) {
            result[i] += lhs(i, j) * rhs[j];
        }
    }
    return result;
}

Vec operator*(const Vec& lhs, const DenseMatrix& rhs)
{
    check_if(lhs.size() == rhs.shape_.first, "Incompatible shapes");

    Vec result(rhs.shape_.second, 0.);
    for (Index i = 0; i < rhs.shape_.first; ++i) {
        auto v = lhs[i];
        for (Index j = 0; j < rhs.shape_.second; ++j) {
            result[j] += v * rhs(i, j);
        }
    }
    return result;
}

DenseMatrix kroneker_product(const DenseMatrix& lhs, const DenseMatrix& rhs)
{
    auto [m, n] = lhs.shape_;
    auto [p, q] = rhs.shape_;
    DenseMatrix result({m * p, n * q});

    for (Index i = 0; i < m; ++i) {
        for (Index j = 0; j < n; ++j) {
            for (Index r = 0; r < p; ++r) {
                for (Index s = 0; s < q; ++s) {
                    result(i * p + r, j * q + s) = lhs(i, j) * rhs(r, s);
                }
            }
        }
    }

    return result;
}

DenseMatrix DenseMatrix::eye(Index dim)
{
    DenseMatrix result({dim, dim});
    for (Index i = 0; i < dim; ++i) {
        result(i, i) = 1.;
    }
    return result;
}