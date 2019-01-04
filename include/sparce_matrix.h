#ifndef MKE2_INCLUDE_SPARCE_MATRIX_H_
#define MKE2_INCLUDE_SPARCE_MATRIX_H_

#include <array>
#include <vector>

#include "abstract_matrix.h"
#include "vec.h"

class DenseMatrix;

class SparceMatrix : public AbstractMatrix {
public:
    using DataContainer = std::vector<Value>;
    using IndexContainer = std::vector<Index>;

    SparceMatrix();
    SparceMatrix(const Shape& shape);
    SparceMatrix(const DenseMatrix& mat);
    SparceMatrix(const DataContainer& data, const IndexContainer& indptr,
                 const IndexContainer& indices, Shape shape);

    Index size() const override;
    Shape shape() const override;
    Index non_zero() const;

    std::pair<Value*, bool> insert(Index i, Index j, Value val = 0);
    Value* find(Index i, Index j);
    const Value* find(Index i, Index j) const;

    friend std::ostream& operator<<(std::ostream& os, const SparceMatrix& obj);

    friend Vec operator*(const SparceMatrix& lhs, const Vec& rhs);
    // friend Vec operator*(const Vec& lhs, const SparceMatrix& rhs);
protected:
    bool index_in_range(Index i, Index j) const;

private:
    DataContainer data_;
    IndexContainer indptr_;
    IndexContainer indices_;
    Shape shape_;
};

#endif // MKE2_INCLUDE_SPARCE_MATRIX_H_