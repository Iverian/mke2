#ifndef MKE2_INCLUDE_POINT3D_HPP_
#define MKE2_INCLUDE_POINT3D_HPP_

#include <array>
#include <functional>
#include <vector>

#include "vec.hpp"

class Point3d;

std::vector<Point3d> combine(const std::initializer_list<double>& x,
                             const std::initializer_list<double>& y,
                             const std::initializer_list<double>& z);

class Point3d : public std::array<double, 3> {
    using Super = std::array<double, 3>;

public:
    using Super::Super;

    Point3d(double x, double y, double z);

    Point3d& operator+=(const Point3d& rhs);
    Point3d& operator-=(const Point3d& rhs);
    Point3d& operator*=(double rhs);
    Point3d& operator/=(double rhs);

    explicit operator Vec() const;

    friend std::ostream& operator<<(std::ostream& os, const Point3d& obj);

    friend bool operator==(const Point3d& lhs, const Point3d& rhs);
    friend bool operator!=(const Point3d& lhs, const Point3d& rhs);

    friend Point3d operator-(const Point3d& obj);
    friend Point3d operator+(const Point3d& lhs, const Point3d& rhs);
    friend Point3d operator-(const Point3d& lhs, const Point3d& rhs);
    friend Point3d operator*(const Point3d& lhs, double rhs);
    friend Point3d operator*(double lhs, const Point3d& rhs);
    friend Point3d operator/(const Point3d& lhs, double rhs);

    friend double dot(const Point3d& lhs, const Point3d& rhs);
    friend Point3d cross(const Point3d& lhs, const Point3d& rhs);
    friend double triple(const Point3d& a, const Point3d& b, const Point3d& c);
    friend double sqr(const Point3d& obj);
    friend double norm(const Point3d& obj);
    friend double dist(const Point3d& lhs, const Point3d& rhs);
    friend bool isnear(const Point3d& lhs, const Point3d& rhs);
};

namespace std {
template <>
struct hash<Point3d> {
    using argument_type = Point3d;

    size_t operator()(const argument_type& key) const
    {
        return (hasher_(key[0]) << 2) ^ (hasher_(key[1]) << 1)
            ^ hasher_(key[0]);
    }

private:
    hash<double> hasher_;
};
};

#endif // MKE2_INCLUDE_POINT3D_HPP_