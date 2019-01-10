#ifndef MKE2_INCLUDE_VEC_HPP_
#define MKE2_INCLUDE_VEC_HPP_

#include <vector>

class Vec : public std::vector<double> {
    using Super = std::vector<double>;

public:
    using Super::Super;
    friend bool operator==(const Vec& lhs, const Vec& rhs);
    friend bool operator!=(const Vec& lhs, const Vec& rhs);
    double dot(const Vec& rhs);
    Vec operator*(const double& rhs);
    Vec operator+(const Vec& rhs);
    Vec operator-(const Vec& rhs);
};

#endif // MKE2_INCLUDE_VEC_HPP_