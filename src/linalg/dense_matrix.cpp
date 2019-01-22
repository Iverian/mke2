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

DenseMatrix::DenseMatrix(Shape shape, const Super& data)
    : Super(data)
    , shape_(shape)
{
    if (data.empty()) {
        Super::resize(size());
    }
}

DenseMatrix::DenseMatrix(Shape shape, Value val)
    : Super()
    , shape_(shape)
{
    resize(size(), val);
}

DenseMatrix::Index DenseMatrix::size() const
{
    return shape_.m * shape_.n;
}

DenseMatrix::Shape DenseMatrix::shape() const
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

#define _(i, j) ((*this)[(shape_.n) * (i) + (j)])

DenseMatrix::Value& DenseMatrix::operator()(Index i, Index j)
{
    check_if(i < shape_.m && j < shape_.n, "Index out of range");

    return _(i, j);
}

const DenseMatrix::Value& DenseMatrix::operator()(Index i, Index j) const

{
    check_if(i < shape_.m && j < shape_.n, "Index out of range");

    return _(i, j);
}

DenseMatrix& DenseMatrix::swap_rows(Index i, Index j)
{
    check_if(i < shape_.m && j < shape_.m, "Index out of range");

    for (Index k = 0; k < shape_.n; ++k) {
        ::swap(_(i, k), _(j, k));
    }
    return *this;
}

DenseMatrix& DenseMatrix::swap_columns(Index i, Index j)
{
    check_if(i < shape_.n && j < shape_.n, "Index out of range");

    for (Index k = 0; k < shape_.m; ++k) {
        ::swap(_(k, i), _(k, j));
    }
    return *this;
}

DenseMatrix& DenseMatrix::append_row(const Vec& v)
{
    check_if(v.size() == shape_.n, "Incompatible shapes");

    resize(size() + shape_.n);
    for (Index j = 0; j < shape_.n; ++j) {
        (*this)(shape_.m, j) = v[j];
    }
    ++shape_.m;

    return *this;
}

DenseMatrix& DenseMatrix::append_column(const Vec& v)
{
    check_if(v.size() == shape_.m, "Incompatible shapes");

    reserve(size() + shape_.m);
    auto pos = begin() + shape_.n;
    for (Index i = 0; i < shape_.m; ++i) {
        insert(pos, v[i]);
        pos += shape_.n;
    }
    return *this;
}

DenseMatrix DenseMatrix::transpose() const
{
    DenseMatrix result({shape_.n, shape_.m});
    for (Index i = 0; i < shape_.m; ++i) {
        for (Index j = 0; j < shape_.n; ++j) {
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
        double max = 0;
        auto s = lhs.size();
        for (DenseMatrix::Index i = 0; i < s; ++i) {
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
    for (DenseMatrix::Index i = 0; i < obj.shape_.m; ++i) {
        os << "[";
        for (DenseMatrix::Index j = 0; j < obj.shape_.n; ++j) {
            os << obj(i, j) << (j + 1 != obj.shape_.n ? ", " : "]");
        }
        os << (i + 1 != obj.shape_.m ? ", " : "]");
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
    check_if(lhs.shape_.n == rhs.shape_.m, "Incompatible shapes");

    DenseMatrix result({lhs.shape_.m, rhs.shape_.n});
    DenseMatrix::Index i, j, k;
    for (i = 0; i < lhs.shape_.m; ++i) {
        for (k = 0; k < lhs.shape_.n; ++k) {
            for (j = 0; j < rhs.shape_.n; ++j) {
                result(i, j) += lhs(i, k) * rhs(k, j);
            }
        }
    }
    return result;
}

DenseMatrix operator*(const DenseMatrix& lhs, DenseMatrix::Value rhs)
{
    auto result = lhs;
    return (result *= rhs);
}

DenseMatrix operator*(DenseMatrix::Value lhs, const DenseMatrix& rhs)
{
    auto result = rhs;
    return (result *= lhs);
}

DenseMatrix operator/(const DenseMatrix& lhs, DenseMatrix::Value rhs)
{
    auto result = lhs;
    return (result /= rhs);
}

Vec operator*(const DenseMatrix& lhs, const Vec& rhs)
{
    check_if(lhs.shape_.n == rhs.size(), "Incompatible shapes");

    Vec result(lhs.shape_.m, 0.);
    for (DenseMatrix::Index i = 0; i < lhs.shape_.m; ++i) {
        for (DenseMatrix::Index j = 0; j < lhs.shape_.n; ++j) {
            result[i] += lhs(i, j) * rhs[j];
        }
    }
    return result;
}

Vec operator*(const Vec& lhs, const DenseMatrix& rhs)
{
    check_if(lhs.size() == rhs.shape_.m, "Incompatible shapes");

    Vec result(rhs.shape_.n, 0.);
    for (DenseMatrix::Index i = 0; i < rhs.shape_.m; ++i) {
        auto v = lhs[i];
        for (DenseMatrix::Index j = 0; j < rhs.shape_.n; ++j) {
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

    for (DenseMatrix::Index i = 0; i < m; ++i) {
        for (DenseMatrix::Index j = 0; j < n; ++j) {
            for (DenseMatrix::Index r = 0; r < p; ++r) {
                for (DenseMatrix::Index s = 0; s < q; ++s) {
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