#ifndef MKE2_INCLUDE_VEC_HPP_
#define MKE2_INCLUDE_VEC_HPP_

#include <vector>

class Vec : public std::vector<double> {
    using Super = std::vector<double>;

public:
    using Super::Super;

    Vec& operator+=(const Vec& rhs);
    Vec& operator-=(const Vec& rhs);
    Vec& operator*=(double rhs);
    Vec& operator/=(double rhs);

    friend bool operator==(const Vec& lhs, const Vec& rhs);
    friend bool operator!=(const Vec& lhs, const Vec& rhs);

    friend Vec operator+(const Vec& lhs, const Vec& rhs);
    friend Vec operator-(const Vec& lhs, const Vec& rhs);
    friend Vec operator*(const Vec& lhs, double rhs);
    friend Vec operator*(double lhs, const Vec& rhs);
    friend Vec operator/(const Vec& lhs, double rhs);

    friend double dot(const Vec& lhs, const Vec& rhs);
    friend double sqr(const Vec& obj);
};

#endif // MKE2_INCLUDE_VEC_HPP_