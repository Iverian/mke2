#ifndef MKE2_INCLUDE_SPARCE_MATRIX_HPP_
#define MKE2_INCLUDE_SPARCE_MATRIX_HPP_

#include <array>
#include <iostream>
#include <vector>

#include "abstract_matrix.hpp"
#include "debug.hpp"
#include "util.hpp"
#include "vec.hpp"

class DenseMatrix;
class SparceMatrix;

class SparceMatrix : public AbstractMatrix {
public:
    using DataContainer = std::vector<Value>;
    using IndexContainer = std::vector<Index>;

    SparceMatrix();
    SparceMatrix(Index side, const DataContainer& vals = DataContainer());

    Index size() const override;
    Shape shape() const override;
    Index non_zero() const;

    const DataContainer& data() const noexcept;
    const DataContainer& diag() const noexcept;
    const IndexContainer& indptr() const noexcept;
    const IndexContainer& indices() const noexcept;

    Value operator()(Index i, Index j) const;
    Value add(Index i, Index j, Value val);
    Value sub(Index i, Index j, Value val);
    Value set(Index i, Index j, Value val);
    void remove_zeroes();

    friend Vec operator*(const SparceMatrix& lhs, const Vec& rhs);
    friend Vec operator*(const Vec& lhs, const SparceMatrix& rhs);
    friend void dot(Vec& result, const SparceMatrix& lhs, const Vec& rhs);
    friend void dot(Vec& result, const Vec& lhs, const SparceMatrix& rhs);

    friend Vec solve(const SparceMatrix& lhs, const Vec& rhs, Vec x0);
    friend Vec solve_p(const SparceMatrix& lhs, const Vec& rhs, Vec x0);

    friend std::ostream& operator<<(std::ostream& os, const SparceMatrix& obj);

protected:
    template <class Callable>
    Value fetch_modify(Index i, Index j, Callable&& f);

private:
    DataContainer data_;
    DataContainer diag_;
    IndexContainer indptr_;
    IndexContainer indices_;
    Index nnz_;
    Index m_;
};

template <class Callable>
SparceMatrix::Value SparceMatrix::fetch_modify(Index i, Index j, Callable&& f)
{
    check_if(i < m_ && j < m_, "Index out of range");

    Value result = 0;

    if (i == j) {
        auto v = f(diag_[i]);
        auto a = isnear(diag_[i], 0);
        auto b = isnear(v, 0);
        result = diag_[i];
        diag_[i] = v;

        if (a && !b) {
            ++nnz_;
        } else if (!a && b) {
            --nnz_;
        }
    } else {
        auto b = std::begin(indices_);
        auto ifirst = b + indptr_[i];
        auto ilast = b + indptr_[i + 1];
        auto p = std::lower_bound(ifirst, ilast, j);
        auto q = std::begin(data_) + Index(p - b);

        if (p != ilast && *p == j) {
            result = *q;
            *q = f(*q);
        } else if (auto val = f(0); !isnear(val, 0)) {
            indices_.emplace(p, j);
            data_.emplace(q, val);
            ++nnz_;

            for (auto r = i; r < m_; ++r) {
                ++indptr_[r + 1];
            }
        }
    }

    return result;
}

#endif // MKE2_INCLUDE_SPARCE_MATRIX_HPP_