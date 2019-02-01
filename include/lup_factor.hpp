#ifndef MKE2_INCLUDE_LUP_FACTOR_HPP_
#define MKE2_INCLUDE_LUP_FACTOR_HPP_

#include "dense_matrix.hpp"
#include "vec.hpp"
#include <cmath>

class LupFactor {
public:
    using PivotType = std::vector<Index>;

    explicit LupFactor(const DenseMatrix& in);
    explicit LupFactor(DenseMatrix&& in);

    const DenseMatrix& mat() const;
    const PivotType& pivot() const;

    LupFactor& factor();
    Vec solve(const Vec& v) const;
    DenseMatrix inverse() const;
    Value det() const;

private:
    Index m_;
    DenseMatrix mat_;
    PivotType pivot_;
};

#endif // MKE2_INCLUDE_LUP_FACTOR_HPP_