#ifndef MKE2_INCLUDE_ABSTRACT_LOCAL_EQ_H_
#define MKE2_INCLUDE_ABSTRACT_LOCAL_EQ_H_

#include "dense_matrix.hpp"
#include "triangulation.hpp"
#include "vec.hpp"

#include <functional>
#include <memory>

class AbstractLocalEq;

using LocalEqGenerator
    = std::function<std::shared_ptr<AbstractLocalEq>(const Triang3d::Node&)>;

class AbstractLocalEq {
public:
    virtual ~AbstractLocalEq() = default;
    virtual DenseMatrix get_internal_mat() const = 0;
    virtual DenseMatrix get_boundary_mat() const = 0;
    virtual Vec get_internal_vec() const = 0;
    virtual Vec get_boundary_vec() const = 0;
};

#endif // MKE2_INCLUDE_ABSTRACT_LOCAL_EQ_H_