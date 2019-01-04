#ifndef MKE2_INCLUDE_VEC_H_
#define MKE2_INCLUDE_VEC_H_

#include <vector>

// using Vec = std::vector<double>;

class Vec : public std::vector<double> {
    using Super = std::vector<double>;

public:
    using Super::Super;
    friend bool operator==(const Vec& lhs, const Vec& rhs);
    friend bool operator!=(const Vec& lhs, const Vec& rhs);
};

#endif // MKE2_INCLUDE_VEC_H_