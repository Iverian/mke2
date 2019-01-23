#ifndef MKE2_INCLUDE_POINT3D_HPP_
#define MKE2_INCLUDE_POINT3D_HPP_

#include <array>
#include <functional>
#include <ostream>
#include <vector>

#include "vec.hpp"
#include "common_types.hpp"

class Point3d;

std::vector<Point3d> combine(const std::initializer_list<Value>& x,
                             const std::initializer_list<Value>& y,
                             const std::initializer_list<Value>& z);

class Point3d : public std::array<Value, 3> {
    static constexpr Index N = 3;
    using Super = std::array<Value, N>;

public:
    using Super::Super;

    Point3d(Value x, Value y, Value z);

    Point3d& operator+=(const Point3d& rhs);
    Point3d& operator-=(const Point3d& rhs);
    Point3d& operator*=(Value rhs);
    Point3d& operator/=(Value rhs);

    explicit operator Vec() const;

    friend std::ostream& operator<<(std::ostream& os, const Point3d& obj);

    friend bool operator==(const Point3d& lhs, const Point3d& rhs);
    friend bool operator!=(const Point3d& lhs, const Point3d& rhs);

    friend Point3d operator-(const Point3d& obj);
    friend Point3d operator+(const Point3d& lhs, const Point3d& rhs);
    friend Point3d operator-(const Point3d& lhs, const Point3d& rhs);
    friend Point3d operator*(const Point3d& lhs, Value rhs);
    friend Point3d operator*(Value lhs, const Point3d& rhs);
    friend Point3d operator/(const Point3d& lhs, Value rhs);

    friend Value dot(const Point3d& lhs, const Point3d& rhs);
    friend Point3d cross(const Point3d& lhs, const Point3d& rhs);
    friend Value triple(const Point3d& a, const Point3d& b, const Point3d& c);
    friend Value sqr(const Point3d& obj);
    friend Value norm(const Point3d& obj);
    friend Point3d unit(const Point3d& obj);
    friend Value dist(const Point3d& lhs, const Point3d& rhs);
    friend bool isnear(const Point3d& lhs, const Point3d& rhs);
};

namespace std {
template <>
struct hash<Point3d> {
    using argument_type = Point3d;

    size_t operator()(const argument_type& key) const
    {
        size_t result = 16651;
        for (auto& i : key) {
            result = (result << 1) ^ hasher_(i);
        }
        return result;
    }

private:
    hash<Value> hasher_;
};
};

#endif // MKE2_INCLUDE_POINT3D_HPP_