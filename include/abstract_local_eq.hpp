#ifndef MKE2_INCLUDE_ABSTRACT_LOCAL_EQ_HPP_
#define MKE2_INCLUDE_ABSTRACT_LOCAL_EQ_HPP_

#include "dense_matrix.hpp"
#include "triangulation.hpp"
#include "vec.hpp"

#include <functional>
#include <memory>

class AbstractLocalEq {
public:
    using Result = std::pair<DenseMatrix, Vec>;

    virtual ~AbstractLocalEq() = default;
    virtual Result
    get_internal(Triangulation::FiniteElementData elem) const = 0;
    virtual Result
    get_boundary(Triangulation::SurfaceElementData elem) const = 0;
};

#endif // MKE2_INCLUDE_ABSTRACT_LOCAL_EQ_HPP_