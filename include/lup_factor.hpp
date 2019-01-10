#ifndef MKE2_INCLUDE_LUP_DECOMPOSITION_HPP_
#define MKE2_INCLUDE_LUP_DECOMPOSITION_HPP_

#include "dense_matrix.hpp"
#include "vec.hpp"

class LupFactor {
public:
    using Index = AbstractMatrix::Index;
    using PivotType = std::vector<Index>;

    explicit LupFactor(const DenseMatrix& in);
    explicit LupFactor(DenseMatrix&& in);

    const DenseMatrix& mat() const;
    const PivotType& pivot() const;

    LupFactor& factor();
    Vec solve(const Vec& v) const;
    DenseMatrix inverse() const;
    double det() const;

private:
    Index m_;
    DenseMatrix mat_;
    PivotType pivot_;
};

#endif // MKE2_INCLUDE_LUP_DECOMPOSITION_HPP_