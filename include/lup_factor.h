#ifndef MKE2_INCLUDE_LUP_DECOMPOSITION_H_
#define MKE2_INCLUDE_LUP_DECOMPOSITION_H_

#include "dense_matrix.h"
#include "vec.h"

class LupFactor {
public:
    using Index = AbstractMatrix::Index;
    using PivotType = std::vector<Index>;

    explicit LupFactor(const DenseMatrix& in);

    const DenseMatrix& mat() const;
    const PivotType& pivot() const;

    Vec solve(const Vec& v) const;
    DenseMatrix invert() const;
    double det() const;

private:
    Index m_;
    DenseMatrix mat_;
    PivotType pivot_;
};

#endif // MKE2_INCLUDE_LUP_DECOMPOSITION_H_