#include <util.hpp>

#include <cmath>

using namespace std;

bool isnear(Value lhs, Value rhs, Tolerance t)
{
    static constexpr Value tol[] = {1e-2, 1e-5, 1e-10, 1e-20};

    return fabs(lhs - rhs) < tol[size_t(t)];
}

bool iszero(Value x, Tolerance t)
{
    return isnear(x, 0., t);
}

Value sqr(Value x)
{
    return x * x;
}
