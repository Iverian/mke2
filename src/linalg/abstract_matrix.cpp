#include <abstract_matrix.h>

using namespace std;

AbstractMatrix::~AbstractMatrix() = default;

AbstractMatrix::Index AbstractMatrix::size() const
{
    auto s = shape();
    return s.m * s.n;
}

ostream& operator<<(ostream& os, const AbstractMatrix::Shape& obj)
{
    return os << "[" << obj.m << ", " << obj.n << "]";
}

bool operator==(const AbstractMatrix::Shape& lhs,
                const AbstractMatrix::Shape& rhs)
{
    return lhs.m == rhs.m && lhs.n == rhs.n;
}

bool operator!=(const AbstractMatrix::Shape& lhs,
                const AbstractMatrix::Shape& rhs)
{
    return !(lhs == rhs);
}