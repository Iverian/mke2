#ifndef MKE2_INCLUDE_SPARCE_MATRIX_HPP_
#define MKE2_INCLUDE_SPARCE_MATRIX_HPP_

#include <array>
#include <vector>

#include "abstract_matrix.hpp"
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
    SparceMatrix(const DataContainer& data, const IndexContainer& indptr,
                 const IndexContainer& indices, Shape shape);

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
    bool index_in_range(Index i, Index j) const;

private:
    DataContainer data_;
    IndexContainer indptr_;
    IndexContainer indices_;
    Shape shape_;
};

template <class Callable>
SparceMatrix::Value SparceMatrix::fetch_modify(Index i, Index j, Callable&& f)
{
    Value result = 0;

    if (index_in_range(i, j)) {
        auto b = begin(indices_);
        auto ifirst = b + indptr_[i];
        auto ilast = b + indptr_[i + 1];
        auto pos = lower_bound(ifirst, ilast, j);
        auto k = Index(pos - b);

        if (pos != end(indices_) && *pos == j) {
            result = data_[k];
            data_[k] = f(data_[k]);
        } else {
            auto new_val = f(0);
            if (!isnear(new_val, 0)) {
                auto dpos = begin(data_) + k;

                indices_.emplace(pos, j);
                data_.emplace(dpos, new_val);

                for (auto r = i; r < shape_.m; ++r) {
                    ++indptr_[r + 1];
                }
            }
        }
    }

    return result;
}

#endif // MKE2_INCLUDE_SPARCE_MATRIX_HPP_