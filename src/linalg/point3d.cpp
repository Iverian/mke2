#include <point3d.h>
#include <util.h>

#include <cmath>

using namespace std;

Point3d::Point3d(double x, double y, double z)
    : Super{x, y, z}
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
Point3d& Point3d::operator*=(double rhs)
{
    for (size_t i = 0; i < size(); ++i) {
        (*this)[i] += rhs;
    }
    return *this;
}
Point3d& Point3d::operator/=(double rhs)
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

bool operator==(const Point3d& lhs, const Point3d& rhs)
{
    auto result = true;
    for (size_t i = 0; i < lhs.size(); ++i) {
        if (!isnear(lhs[i], rhs[i], Tolerance::DOUBLE)) {
            result = false;
            break;
        }
    }
    return result;
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

Point3d operator*(const Point3d& lhs, double rhs)
{
    auto result = lhs;
    return (result *= rhs);
}

Point3d operator*(double lhs, const Point3d& rhs)
{
    auto result = rhs;
    return (result *= lhs);
}

Point3d operator/(const Point3d& lhs, double rhs)
{
    auto result = lhs;
    return (result /= rhs);
}

double dot(const Point3d& lhs, const Point3d& rhs)
{
    double result = 0;
    for (size_t i = 0; i < lhs.size(); ++i) {
        result += lhs[i] * rhs[i];
    }
    return result;
}

double sqr(const Point3d& obj)
{
    return dot(obj, obj);
}

double norm(const Point3d& obj)
{
    return sqrt(sqr(obj));
}

double dist(const Point3d& lhs, const Point3d& rhs)
{
    return norm(rhs - lhs);
}

bool isnear(const Point3d& lhs, const Point3d& rhs)
{
    return isnear(dist(lhs, rhs), 0);
}

vector<Point3d> combine(const initializer_list<double>& x,
                        const initializer_list<double>& y,
                        const initializer_list<double>& z)
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