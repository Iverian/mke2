#ifndef MKE2_INCLUDE_SPARCE_MATRIX_HPP_
#define MKE2_INCLUDE_SPARCE_MATRIX_HPP_

#include <array>
#include <vector>

#include "abstract_matrix.hpp"
#include "debug.hpp"
#include "util.hpp"
#include "vec.hpp"

class DenseMatrix;

class SparceMatrix : public AbstractMatrix {
public:
    using DataContainer = std::vector<Value>;
    using IndexContainer = std::vector<Index>;

    SparceMatrix();
    SparceMatrix(const Shape& shape);
    SparceMatrix(const DenseMatrix& mat);

    Index size() const override;
    Shape shape() const override;
    Index non_zero() const;

    Value operator()(Index i, Index j) const;
    Value fetch_add(Index i, Index j, Value val);
    Value fetch_set(Index i, Index j, Value val);
    void clean_up();

    Value* find(Index i, Index j);
    const Value* find(Index i, Index j) const;

    friend std::ostream& operator<<(std::ostream& os, const SparceMatrix& obj);
    friend Vec operator*(const SparceMatrix& lhs, const Vec& rhs);
    friend void dot(double* const result, const SparceMatrix& lhs,
                    const Vec& rhs);

    friend Vec solve_cg(const SparceMatrix& lhs, const Vec& rhs, Vec x0,
                        const double tol);

protected:
    template <class Callable>
    Value fetch_modify(Index i, Index j, Callable&& f);

private:
    DataContainer data_;
    IndexContainer indptr_;
    IndexContainer indices_;
    Shape shape_;
};

template <class Callable>
SparceMatrix::Value SparceMatrix::fetch_modify(Index i, Index j, Callable&& f)
{
    check_if(i < shape_.m && j < shape_.n, "Index out of range");

    Value result = 0;

    auto b = begin(indices_);
    auto p = lower_bound(b + indptr_[i], b + indptr_[i + 1], j);
    auto q = begin(data_) + Index(p - b);

    if (p != end(indices_) && *p == j) {
        result = *q;
        *q = f(*q);
    } else if (auto new_val = f(0); !isnear(new_val, 0)) {
        indices_.emplace(p, j);
        data_.emplace(q, new_val);

        for (auto r = i; r < shape_.m; ++r) {
            ++indptr_[r + 1];
        }
    }

    return result;
}

#endif // MKE2_INCLUDE_SPARCE_MATRIX_HPP_