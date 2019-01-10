#ifndef MKE2_INCLUDE_ABSTRACT_MATRIX_H_
#define MKE2_INCLUDE_ABSTRACT_MATRIX_H_

#include <ostream>

class AbstractMatrix {
public:
    using Value = double;
    using Index = size_t;

    struct Shape {
        Index m;
        Index n;

        friend std::ostream& operator<<(std::ostream& os, const Shape& obj);
        friend bool operator==(const Shape& lhs, const Shape& rhs);
        friend bool operator!=(const Shape& lhs, const Shape& rhs);
    };

    virtual ~AbstractMatrix();

    virtual Index size() const;
    virtual Shape shape() const = 0;
};

#endif // MKE2_INCLUDE_ABSTRACT_MATRIX_H_