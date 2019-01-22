#ifndef MKE2_INCLUDE_SPARCE_MATRIX_HPP_
#define MKE2_INCLUDE_SPARCE_MATRIX_HPP_

#include <array>
#include <iostream>
#include <iterator>
#include <set>
#include <type_traits>
#include <vector>

#include "abstract_matrix.hpp"
#include "debug.hpp"
#include "dok_matrix.hpp"
#include "util.hpp"
#include "vec.hpp"

class CsrMatrix : public AbstractMatrix {
public:
    using DataContainer = std::vector<Value>;
    using IndexContainer = std::vector<Index>;

    CsrMatrix();
    CsrMatrix(Index side, const DataContainer& vals = DataContainer());
    explicit CsrMatrix(const DokMatrix& dok);

    Index size() const override;
    Shape shape() const override;
    Index non_zero() const;

    const DataContainer& data() const noexcept;
    const IndexContainer& indptr() const noexcept;
    const IndexContainer& indices() const noexcept;

    Value operator()(Index i, Index j) const;
    friend Vec operator*(const CsrMatrix& lhs, const Vec& rhs);
    friend Vec operator*(const Vec& lhs, const CsrMatrix& rhs);
    friend void dot(Vec& result, const CsrMatrix& lhs, const Vec& rhs);
    friend void dot(Vec& result, const Vec& lhs, const CsrMatrix& rhs);

    friend Vec solve(const CsrMatrix& lhs, const Vec& rhs, Vec x0);
    friend std::ostream& operator<<(std::ostream& os, const CsrMatrix& obj);

private:
    DataContainer data_;
    IndexContainer indptr_;
    IndexContainer indices_;
    Index m_;
};

#endif // MKE2_INCLUDE_SPARCE_MATRIX_HPP_