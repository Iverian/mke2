#include <util.hpp>

#include <cmath>

using namespace std;

bool isnear(double lhs, double rhs, Tolerance t)
{
    return fabs(lhs - rhs) < tol[size_t(t)];
}

double sqr(double x)
{
    return x * x;
}
