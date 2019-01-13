#include <util.hpp>

#include <cmath>

using namespace std;

bool isnear(double lhs, double rhs, Tolerance t)
{
    static constexpr double tol[] = {1e-2, 1e-5, 1e-8, 1e-13};

    return fabs(lhs - rhs) < tol[size_t(t)];
}

double sqr(double x)
{
    return x * x;
}
