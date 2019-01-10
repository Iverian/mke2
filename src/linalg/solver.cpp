#include <sparse_matrix.hpp>
#include <vec.hpp>

using namespace std;

Vec CG(const SparseMatrix& mat, const Vec& right_side)
{
    int max_step = 10000;
    int step = 0;
    double accuracy = 1e-5;

    Vec x0(right_side.size(), 0), 
        r0(right_side), 
        z0(right_side), 
        xK(right_side.size()), 
        rK(right_side.size()), 
        zk(right_side.size()),
        temp(right_side.size());
    double alpha, beta;

    do {

        temp = mat * z0;

        alpha = (r0.dot(r0)) / (z0.dot(temp));

        xK = x0 + (z0 * alpha);

        rK = r0 - (temp * alpha);

        if (((rK.dot(rK)) / (right_side.dot(right_side))) < accuracy) 
            break;

        beta = (rK.dot(rK)) / (r0.dot(r0));

        zK = rK + z0 * beta;

        step++;

    } while(step < max_step);
}

