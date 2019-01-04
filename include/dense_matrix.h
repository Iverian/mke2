#ifndef MKE2_INCLUDE_DENSE_MATRIX_H_
#define MKE2_INCLUDE_DENSE_MATRIX_H_

#include "abstract_matrix.h"
#include "vec.h"

#include <vector>

class DenseMatrix : public AbstractMatrix, std::vector<AbstractMatrix::Value> {
    using Super = std::vector<double>;

public:
    using Super::begin;
    using Super::cbegin;
    using Super::cend;
    using Super::const_iterator;
    using Super::end;
    using Super::iterator;
    using Super::operator[];

    DenseMatrix();
    DenseMatrix(Shape shape, const Super& data = Super());

    Index size() const override;
    Shape shape() const override;
    Value& operator()(Index i, Index j) noexcept;
    const Value& operator()(Index i, Index j) const noexcept;

    void swap_rows(Index i, Index j);
    void swap_columns(Index i, Index j);

    friend std::ostream& operator<<(std::ostream& os, const DenseMatrix& obj);
    friend bool operator==(const DenseMatrix& lhs, const DenseMatrix& rhs);
    friend bool operator!=(const DenseMatrix& lhs, const DenseMatrix& rhs);

    // DenseMatrix& operator+=(const DenseMatrix& rhs);
    // DenseMatrix& operator-=(const DenseMatrix& rhs);
    // DenseMatrix& operator*=(const DenseMatrix& rhs);
    // DenseMatrix& operator/=(const DenseMatrix& rhs);
    // DenseMatrix& operator*=(Value rhs);
    // DenseMatrix& operator/=(Value rhs);

    // friend DenseMatrix operator+(const DenseMatrix& lhs,
    //                              const DenseMatrix& rhs);
    // friend DenseMatrix operator-(const DenseMatrix& lhs,
    //                              const DenseMatrix& rhs);
    // friend DenseMatrix operator*(const DenseMatrix& lhs,
    //                              const DenseMatrix& rhs);
    // friend DenseMatrix operator/(const DenseMatrix& lhs,
    //                              const DenseMatrix& rhs);
    // friend DenseMatrix operator*(const DenseMatrix& lhs, Value rhs);
    // friend DenseMatrix operator*(Value lhs, const DenseMatrix& rhs);
    // friend DenseMatrix operator/(const DenseMatrix& lhs, Value rhs);
    friend Vec operator*(const DenseMatrix& lhs, const Vec& rhs);
    friend Vec operator*(const Vec& lhs, const DenseMatrix& rhs);

private:
    Shape shape_;
};

#endif // MKE2_INCLUDE_DENSE_MATRIX_H_