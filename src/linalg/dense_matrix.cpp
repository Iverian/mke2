#include <dense_matrix.h>
#include <util.h>

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
}

DenseMatrix::Index DenseMatrix::size() const
{
    return shape_.m * shape_.n;
}

DenseMatrix::Shape DenseMatrix::shape() const
{
    return shape_;
}

#define _(i, j) ((*this)[(shape_.m) * (i) + (j)])

DenseMatrix::Value& DenseMatrix::operator()(Index i, Index j) noexcept
{
    return _(i, j);
}

const DenseMatrix::Value& DenseMatrix::operator()(Index i, Index j) const
    noexcept
{
    return _(i, j);
}

void DenseMatrix::swap_rows(Index i, Index j)
{
    for (Index k = 0; k < shape_.n; ++k) {
        ::swap(_(i, k), _(j, k));
    }
}

void DenseMatrix::swap_columns(Index i, Index j)
{
    for (Index k = 0; k < shape_.m; ++k) {
        ::swap(_(k, i), _(k, j));
    }
}

#undef _

bool operator==(const DenseMatrix& lhs, const DenseMatrix& rhs)
{
    auto result = lhs.shape_ == rhs.shape_;
    if (result) {
        auto s = lhs.size();
        for (DenseMatrix::Index i = 0; i < s; ++i) {
            if (!isnear(lhs[i], rhs[i], Tolerance::DOUBLE)) {
                result = false;
                break;
            }
        }
    }
    return result;
}

bool operator!=(const DenseMatrix& lhs, const DenseMatrix& rhs)
{
    return !(lhs == rhs);
}

ostream& operator<<(ostream& os, const DenseMatrix& obj)
{
    os << "{ \"shape\": " << obj.shape_ << ", \"data\": [";
    for (DenseMatrix::Index i = 0; i < obj.shape_.m; ++i) {
        for (DenseMatrix::Index j = 0; j < obj.shape_.n; ++j) {
            os << obj(i, j) << (j + 1 != obj.shape_.n ? ", " : "]");
        }
        os << (i + 1 != obj.shape_.m ? ", " : "]");
    }
    return os << "}";
}

Vec operator*(const DenseMatrix& lhs, const Vec& rhs)
{
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
    Vec result(rhs.shape_.n, 0.);
    for (DenseMatrix::Index i = 0; i < rhs.shape_.m; ++i) {
        auto v = lhs[i];
        for (DenseMatrix::Index j = 0; j < rhs.shape_.n; ++j) {
            result[j] += v * rhs(i, j);
        }
    }
    return result;
}