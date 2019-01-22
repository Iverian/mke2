#include <csr_matrix.hpp>
#include <debug.hpp>

#include <cmath>

static constexpr size_t max_iter = 10000;

using namespace std;

void cax_y(Vec& res, double a, const Vec& x, const Vec& y);
void cax_by_z(Vec& res, double a, const Vec& x, double b, const Vec& y,
              const Vec& z);
double cdist(const Vec& lhs, const Vec& rhs);

Vec solve(const CsrMatrix& lhs, const Vec& rhs, Vec x0)
{
    auto [m, n] = lhs.shape();
    check_if(rhs.size() == n, "Incompatible shapes");

    size_t step = 0;

    if (x0.empty()) {
        x0.resize(m, 0);
    }
    auto rsqr = sqr(rhs);
    auto r = rhs - lhs * x0;
    auto r0 = r;
    auto x = move(x0);

    // Vec p(m, 0.);
    Vec v(m, 0.);
    Vec s(m, 0.);
    Vec t(m, 0.);
    Vec y(m, 0.);
    Vec z(m, 0.);
    Vec q(m, 0.);
    Vec p(m, 0.);

    double mdiff = numeric_limits<double>::max();
    size_t mind = 0;
    bool success = false;

    double pdiff = 0, diff = 0, res = 0, u = 1, a = 1, w = 1, b = 0;
    do {
        auto up = u;
        u = dot(r0, r);

        b = (u / up) * (a / w);
        // p = r + b * (p - w * v)
        cax_by_z(p, b, p, -w * b, v, r);
        // v = lhs * p
        dot(v, lhs, p);

        a = u / dot(r0, v);

        // s = r - a * v
        cax_y(s, -a, v, r);
        // t = lhs * s
        dot(t, lhs, s);
        w = dot(t, s) / sqr(t);

        cax_y(r, -w, t, s);
        res = sqrt(sqr(r) / rsqr);
        if (!isfinite(res)) {
            break;
        }
        // if (iszero(res, Tolerance::TRIPLE)) {
        //     success = true;
        //     break;
        // }

        // x = a * p + w * s + x
        cax_by_z(x, a, p, w, s, x);
        // r = s - w * t

        dot(y, lhs, x);
        pdiff = diff;
        diff = cdist(y, rhs);
        if (diff < mdiff) {
            mdiff = diff;
            mind = step;
        }
        if (iszero(diff)) {
            success = true;
            break;
        }
    } while (++step < max_iter);

    // diff = sqrt(sqr(lhs * x - rhs));
    if (success) {
        cout << "Iteration converged";
    } else {
        cout << "Iteration did not converge";
    }
    cout << ": res=" << res << ", |Ax - b|=" << diff << ", step=" << step
         << ", mind=" << mind << endl;

    return x;
}

void cax_y(Vec& res, double a, const Vec& x, const Vec& y)
{
    CsrMatrix::Index i, n = res.size();

    for (i = 0; i < n; ++i) {
        res[i] = a * x[i] + y[i];
    }
}

void cax_by_z(Vec& res, double a, const Vec& x, double b, const Vec& y,
              const Vec& z)
{
    CsrMatrix::Index i, n = res.size();

    for (i = 0; i < n; ++i) {
        res[i] = a * x[i] + b * y[i] + z[i];
    }
}
