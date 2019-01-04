#ifndef MKE2_INCLUDE_ABSTRACT_MATRIX_H_
#define MKE2_INCLUDE_ABSTRACT_MATRIX_H_

#include <cstdlib>
#include <ostream>

class AbstractMatrix {
public:
    using Value = double;
    using Index = size_t;

    struct Shape {
        Index m;
        Index n;
    };

    virtual ~AbstractMatrix();
    virtual Index size() const;
    virtual Shape shape() const = 0;
};

std::ostream& operator<<(std::ostream& os, const AbstractMatrix::Shape& obj);
bool operator==(const AbstractMatrix::Shape& lhs,
                const AbstractMatrix::Shape& rhs);
bool operator!=(const AbstractMatrix::Shape& lhs,
                const AbstractMatrix::Shape& rhs);

#endif // MKE2_INCLUDE_ABSTRACT_MATRIX_H_