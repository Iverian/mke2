#include <mv2_export.hpp>
#include <algorithm>

using namespace std;

void mv2_export(ostream& os, const Triangulation& t, const Vec& values)
{
    auto m = t.nodes().size();
    os << m << " " << 3 << " " << 3 << " Ux Uy Uz" << endl;

    for (size_t j = 0; j < t.nodes().size(); j++)
    {
        auto n = find_if(t.nodes().begin(), t.nodes().end(),[j](auto& item){
            auto [a,b] = item;
            return (b == j);
        });
        auto [p, i] = *n;
        auto& x = values[i];
        auto& y = values[m + i];
        auto& z = values[2 * m + i];
        
        //temp = to_string(i) + " " + to_string(p[0]) + " " + to_string(p[1]) + " " + to_string(p[2]) 
        //    + " " + to_string(x) + " " + to_string(y) + " " + to_string(z) + "\n";
        os << i << " " << p[0] << " " << p[1] << " " << p[2] << " " << x << " "
           << y << " " << z << endl;
    }
    //for (auto& n : t.nodes()) {
    //    auto [p, i] = n;
    //    auto& x = values[i];
    //    auto& y = values[m + i];
    //    auto& z = values[2 * m + i];
    //    
        //temp = to_string(i) + " " + to_string(p[0]) + " " + to_string(p[1]) + " " + to_string(p[2]) 
        //    + " " + to_string(x) + " " + to_string(y) + " " + to_string(z) + "\n";
    //    os << i << " " << p[0] << " " << p[1] << " " << p[2] << " " << x << " "
    //       << y << " " << z << endl;
        //s.push_back(temp);
    //}
}