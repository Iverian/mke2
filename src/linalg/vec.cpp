#include <debug.hpp>
#include <util.hpp>
#include <vec.hpp>

#include <iterator>

using namespace std;

Vec& Vec::operator+=(const Vec& rhs)
{
    check_if(size() == rhs.size(), "Dimensions are not equal");

    for (size_t i = 0; i < size(); ++i) {
        (*this)[i] += rhs[i];
    }
    return *this;
}

Vec& Vec::operator-=(const Vec& rhs)
{
    check_if(size() == rhs.size(), "Dimensions are not equal");

    for (size_t i = 0; i < size(); ++i) {
        (*this)[i] -= rhs[i];
    }
    return *this;
}
Vec& Vec::operator*=(double rhs)
{
    for (size_t i = 0; i < size(); ++i) {
        (*this)[i] *= rhs;
    }
    return *this;
}

Vec& Vec::operator/=(double rhs)
{
    for (size_t i = 0; i < size(); ++i) {
        (*this)[i] /= rhs;
    }
    return *this;
}

bool operator==(const Vec& lhs, const Vec& rhs)
{
    return lhs.size() == rhs.size()
        && iszero(cdist(lhs, rhs), Tolerance::DOUBLE);
}

double cdist(const Vec& lhs, const Vec& rhs)
{
    check_if(lhs.size() == rhs.size(), "Dimensions are not equal");

    double cur, max = 0;
    for (size_t i = 0; i < lhs.size(); ++i) {
        cur = fabs(rhs[i] - lhs[i]);
        if (cur > max) {
            max = cur;
        }
    }
    return max;
}

bool operator!=(const Vec& lhs, const Vec& rhs)
{
    return !(lhs == rhs);
}

Vec operator+(const Vec& lhs, const Vec& rhs)
{
    auto result = lhs;
    return (result += rhs);
}

Vec operator-(const Vec& lhs, const Vec& rhs)
{
    auto result = lhs;
    return (result -= rhs);
}

Vec operator*(const Vec& lhs, double rhs)
{
    auto result = lhs;
    return (result *= rhs);
}

Vec operator*(double lhs, const Vec& rhs)
{
    auto result = rhs;
    return (result *= lhs);
}

Vec operator/(const Vec& lhs, double rhs)
{
    auto result = lhs;
    return (result /= rhs);
}

double dot(const Vec& lhs, const Vec& rhs)
{
    check_if(lhs.size() == rhs.size(), "Dimensions are not equal");

    double result = 0;
    for (size_t i = 0; i < lhs.size(); ++i) {
        result += lhs[i] * rhs[i];
    }
    return result;
}

double sqr(const Vec& obj)
{
    return dot(obj, obj);
}

ostream& operator<<(ostream& os, const Vec& obj)
{
    os << "[";
    for (auto i = begin(obj); i != end(obj); ++i) {
        os << *i << (next(i) != end(obj) ? ", " : "]");
    }
    return os;
}