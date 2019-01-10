#include <util.hpp>
#include <vec.hpp>

bool operator==(const Vec& lhs, const Vec& rhs)
{
    auto result = lhs.size() == rhs.size();
    if (result) {
        for (size_t i = 0; i < lhs.size(); ++i) {
            if (!isnear(lhs[i], rhs[i], Tolerance::DOUBLE)) {
                result = false;
                break;
            }
        }
    }
    return result;
}

bool operator!=(const Vec& lhs, const Vec& rhs)
{
    return !(lhs == rhs);
}

double Vec::dot(const Vec& rhs)
{
    check_if((*this).size() == rhs.size(),
        "Dimensions are not equal");

    double res = 0;

    for (size_t i = 0; i < rhs.size(); i++) {
        res += (*this)[i] * rhs[i];
    }

    return res;
}

Vec Vec::operator*(const double& rhs)
{
    Vec res((*this).size());

    for (size_t i = 0; i < res.size(); i++) {
        res[i] = rhs * (*this)[i];
    }

    return res;
}

Vec Vec::operator+(const Vec& rhs)
{
    check_if((*this).size() == rhs.size(),
        "Dimensions are not equal");
    
    Vec res(rhs.size());

    for (size_t i = 0; i < rhs.size(); i++) {
        res[i] = (*this)[i] + rhs[i];
    }

    return res;
}

Vec Vec::operator-(const Vec& rhs)
{
    check_if((*this).size() == rhs.size(),
        "Dimensions are not equal");
    
    Vec res(rhs.size());

    for (size_t i = 0; i < rhs.size(); i++) {
        res[i] = (*this)[i] - rhs[i];
    }

    return res;
}