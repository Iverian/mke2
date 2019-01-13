#ifndef MKE2_INCLUDE_GLOBAL_EQ_BUILDER_HPP_
#define MKE2_INCLUDE_GLOBAL_EQ_BUILDER_HPP_

#include "abstract_local_eq.hpp"
#include "sparce_matrix.hpp"
#include "triangulation.hpp"

#include <tuple>

class GlobalEqBuilder {
public:
    using Index = AbstractMatrix::Index;

    explicit GlobalEqBuilder(const Triangulation& triang,
                             std::shared_ptr<AbstractLocalEq> gen);

    GlobalEqBuilder& get();

    SparceMatrix& mat();
    Vec& vec();

private:
    const Triangulation& t_;
    std::shared_ptr<AbstractLocalEq> g_;
    Index m_;
    SparceMatrix lhs_;
    Vec rhs_;
};

#endif // MKE2_INCLUDE_GLOBAL_EQ_BUILDER_HPP_