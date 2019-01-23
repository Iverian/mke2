#ifndef MKE2_INCLUDE_ABSTRACT_MATRIX_HPP_
#define MKE2_INCLUDE_ABSTRACT_MATRIX_HPP_

#include "common_types.hpp"

#include <ostream>

class AbstractMatrix {
public:
    virtual ~AbstractMatrix();

    virtual Index size() const;
    virtual Index2d shape() const = 0;
};

#endif // MKE2_INCLUDE_ABSTRACT_MATRIX_HPP_