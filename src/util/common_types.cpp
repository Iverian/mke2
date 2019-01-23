#include <common_types.hpp>

using namespace std;

ostream& operator<<(ostream& os, const Index2d& obj)
{
    return os << "[" << obj.first << ", " << obj.second << "]";
}
