#ifndef MKE2_INCLUDE_SPARCE_MATRIX_HPP_
#define MKE2_INCLUDE_SPARCE_MATRIX_HPP_

#include <array>
#include <vector>

#include "abstract_matrix.hpp"
#include "vec.hpp"

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

    Value operator()(Index i, Index j) const;
    Value fetch_add(Index i, Index j, Value val);
    void clean_up();

    Value* find(Index i, Index j);
    const Value* find(Index i, Index j) const;

    friend std::ostream& operator<<(std::ostream& os, const SparceMatrix& obj);
    friend Vec operator*(const SparceMatrix& lhs, const Vec& rhs);

protected:
    bool index_in_range(Index i, Index j) const;

private:
    DataContainer data_;
    IndexContainer indptr_;
    IndexContainer indices_;
    Shape shape_;
};

#endif // MKE2_INCLUDE_SPARCE_MATRIX_HPP_