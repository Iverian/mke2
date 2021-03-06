#include <point3d.hpp>
#include <util.hpp>

#include <cmath>

using namespace std;

Point3d::Point3d(Value x, Value y, Value z)
    : Super {x, y, z}
{
}

Point3d& Point3d::operator+=(const Point3d& rhs)
{
    for (size_t i = 0; i < size(); ++i) {
        (*this)[i] += rhs[i];
    }
    return *this;
}

Point3d& Point3d::operator-=(const Point3d& rhs)
{
    for (size_t i = 0; i < size(); ++i) {
        (*this)[i] -= rhs[i];
    }
    return *this;
}
Point3d& Point3d::operator*=(Value rhs)
{
    for (size_t i = 0; i < size(); ++i) {
        (*this)[i] += rhs;
    }
    return *this;
}
Point3d& Point3d::operator/=(Value rhs)
{
    for (size_t i = 0; i < size(); ++i) {
        (*this)[i] /= rhs;
    }
    return *this;
}

Point3d::operator Vec() const
{
    return Vec(::begin(*this), ::end(*this));
}

ostream& operator<<(ostream& os, const Point3d& obj)
{
    return os << "[" << obj[0] << ", " << obj[1] << ", " << obj[2] << "]";
}

bool operator==(const Point3d& lhs, const Point3d& rhs)
{
    Value max = 0;
    for (size_t i = 0; i < lhs.size(); ++i) {
        if (auto cur = fabs(rhs[i] - lhs[i]); cur > max) {
            max = cur;
        }
    }
    return iszero(max, Tolerance::DOUBLE);
}

bool operator!=(const Point3d& lhs, const Point3d& rhs)
{
    return !(lhs == rhs);
}

Point3d operator-(const Point3d& obj)
{
    auto result = obj;
    return (result *= -1.);
}

Point3d operator+(const Point3d& lhs, const Point3d& rhs)
{
    auto result = lhs;
    return (result += rhs);
}

Point3d operator-(const Point3d& lhs, const Point3d& rhs)
{
    auto result = lhs;
    return (result -= rhs);
}

Point3d operator*(const Point3d& lhs, Value rhs)
{
    auto result = lhs;
    return (result *= rhs);
}

Point3d operator*(Value lhs, const Point3d& rhs)
{
    auto result = rhs;
    return (result *= lhs);
}

Point3d operator/(const Point3d& lhs, Value rhs)
{
    auto result = lhs;
    return (result /= rhs);
}

Point3d cross(const Point3d& lhs, const Point3d& rhs)
{
    return Point3d(lhs[1] * rhs[2] - lhs[2] * rhs[1],
                   lhs[2] * rhs[0] - lhs[0] * rhs[2],
                   lhs[0] * rhs[1] - lhs[1] * rhs[0]);
}

Value dot(const Point3d& lhs, const Point3d& rhs)
{
    Value result = 0;
    for (size_t i = 0; i < lhs.size(); ++i) {
        result += lhs[i] * rhs[i];
    }
    return result;
}

Value triple(const Point3d& a, const Point3d& b, const Point3d& c)
{
    return dot(a, cross(b, c));
}

Value sqr(const Point3d& obj)
{
    return dot(obj, obj);
}

Value norm(const Point3d& obj)
{
    return sqrt(sqr(obj));
}

Point3d unit(const Point3d& obj)
{
    return obj / norm(obj);
}

Value dist(const Point3d& lhs, const Point3d& rhs)
{
    return norm(rhs - lhs);
}

bool isnear(const Point3d& lhs, const Point3d& rhs)
{
    return isnear(dist(lhs, rhs), 0);
}

vector<Point3d> combine(const initializer_list<Value>& x,
                        const initializer_list<Value>& y,
                        const initializer_list<Value>& z)
{
    vector<Point3d> result;

    for (auto& i : x) {
        for (auto& j : y) {
            for (auto& k : z) {
                result.emplace_back(i, j, k);
            }
        }
    }

    result.shrink_to_fit();
    return result;
}