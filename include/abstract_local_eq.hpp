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
    get_internal(Triangulation::FiniteElement::Data elem) const = 0;
    virtual Result
    get_boundary(Triangulation::SurfaceElement::Data elem) const = 0;
};

using LocalEqGen = std::function<std::pair<DenseMatrix, Vec>(
    const Triangulation&, const Triangulation::FiniteElement&)>;

#endif // MKE2_INCLUDE_ABSTRACT_LOCAL_EQ_HPP_