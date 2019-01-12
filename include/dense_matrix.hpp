#ifndef MKE2_INCLUDE_DENSE_MATRIX_HPP_
#define MKE2_INCLUDE_DENSE_MATRIX_HPP_

#include "abstract_matrix.hpp"
#include "vec.hpp"

#include <ostream>
#include <vector>

class DenseMatrix : public AbstractMatrix, std::vector<AbstractMatrix::Value> {
    using Super = std::vector<AbstractMatrix::Value>;

public:
    using Super::begin;
    using Super::cbegin;
    using Super::cend;
    using Super::const_iterator;
    using Super::end;
    using Super::iterator;
    using Super::operator[];

    static DenseMatrix eye(Index dim);

    DenseMatrix();
    DenseMatrix(Shape shape, const Super& data = Super());
    DenseMatrix(Shape shape, Value val);
    template <class InputIt>
    DenseMatrix(Shape shape, InputIt first, InputIt last);

    DenseMatrix(const DenseMatrix&) = default;
    DenseMatrix(DenseMatrix&&) noexcept = default;
    DenseMatrix& operator=(const DenseMatrix&) = default;
    DenseMatrix& operator=(DenseMatrix&&) noexcept = default;

    Index size() const override;
    Shape shape() const override;
    Value& operator()(Index i, Index j) noexcept;
    const Value& operator()(Index i, Index j) const noexcept;

    Super& data();
    const Super& data() const;

    DenseMatrix& swap_rows(Index i, Index j);
    DenseMatrix& swap_columns(Index i, Index j);
    DenseMatrix& append_row(const Vec& v);
    DenseMatrix& append_column(const Vec& v);

    DenseMatrix transpose() const;

    DenseMatrix& operator+=(const DenseMatrix& rhs);
    DenseMatrix& operator-=(const DenseMatrix& rhs);
    DenseMatrix& operator*=(Value rhs);
    DenseMatrix& operator/=(Value rhs);

    friend std::ostream& operator<<(std::ostream& os, const DenseMatrix& obj);

    friend bool operator==(const DenseMatrix& lhs, const DenseMatrix& rhs);
    friend bool operator!=(const DenseMatrix& lhs, const DenseMatrix& rhs);

    friend DenseMatrix operator+(const DenseMatrix& lhs,
                                 const DenseMatrix& rhs);
    friend DenseMatrix operator-(const DenseMatrix& lhs,
                                 const DenseMatrix& rhs);
    friend DenseMatrix operator*(const DenseMatrix& lhs,
                                 const DenseMatrix& rhs);
    friend DenseMatrix operator*(const DenseMatrix& lhs, Value rhs);
    friend DenseMatrix operator*(Value lhs, const DenseMatrix& rhs);
    friend DenseMatrix operator/(const DenseMatrix& lhs, Value rhs);
    friend Vec operator*(const DenseMatrix& lhs, const Vec& rhs);
    friend Vec operator*(const Vec& lhs, const DenseMatrix& rhs);

    friend DenseMatrix kroneker_product(const DenseMatrix& lhs,
                                        const DenseMatrix& rhs);

private:
    Shape shape_;
};

template <class InputIt>
DenseMatrix::DenseMatrix(Shape shape, InputIt first, InputIt last)
    : Super(first, last)
    , shape_(shape)
{
}

#endif // MKE2_INCLUDE_DENSE_MATRIX_HPP_