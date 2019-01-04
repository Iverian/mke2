#ifndef MKE2_INCLUDE_POINT3D_H_
#define MKE2_INCLUDE_POINT3D_H_

#include <array>
#include <vector>

#include "vec.h"

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

    friend bool operator==(const Point3d& lhs, const Point3d& rhs);
    friend bool operator!=(const Point3d& lhs, const Point3d& rhs);

    friend Point3d operator-(const Point3d& obj);
    friend Point3d operator+(const Point3d& lhs, const Point3d& rhs);
    friend Point3d operator-(const Point3d& lhs, const Point3d& rhs);
    friend Point3d operator*(const Point3d& lhs, double rhs);
    friend Point3d operator*(double lhs, const Point3d& rhs);
    friend Point3d operator/(const Point3d& lhs, double rhs);

    friend double dot(const Point3d& lhs, const Point3d& rhs);
    friend double sqr(const Point3d& obj);
    friend double norm(const Point3d& obj);
    friend double dist(const Point3d& lhs, const Point3d& rhs);
    friend bool isnear(const Point3d& lhs, const Point3d& rhs);
};

#endif // MKE2_INCLUDE_POINT3D_H_