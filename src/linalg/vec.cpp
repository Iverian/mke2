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