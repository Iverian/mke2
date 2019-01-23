#include <constants.hpp>
#include <csr_matrix.hpp>
#include <debug.hpp>


#include <cmath>

using namespace std;

void cax_y(Vec& res, double a, const Vec& x, const Vec& y);
void cax_by_z(Vec& res, double a, const Vec& x, double b, const Vec& y,
              const Vec& z);
void dot_diag(Vec& res, const CsrMatrix& lhs, const Vec& rhs);

Vec solve(const CsrMatrix& lhs, const Vec& rhs, Vec x0)
{
    auto [m, n] = lhs.shape();
    check_if(rhs.size() == n, "Incompatible shapes");

    size_t step = 0;

    if (x0.empty()) {
        x0.resize(m, 0);
    }
    auto rnorm = cnorm(rhs);
    auto r = rhs - lhs * x0;
    auto r0 = r;
    auto x = move(x0);

    Vec v(m, 0.);
    Vec s(m, 0.);
    Vec t(m, 0.);
    Vec p(m, 0.);
    Vec ax(m, 0.);

    bool success = false;

    double diff = 0, res = 0, u = 1, a = 1, w = 1, b = 0;
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

        // r = s - w * t
        cax_y(r, -w, t, s);
        res = cnorm(r) / rnorm;
        if (!isnormal(res)) {
            break;
        }

        // x = a * p + w * s + x
        cax_by_z(x, a, p, w, s, x);

        dot(ax, lhs, x);
        diff = cdist(ax, rhs);
        if (iszero(diff)) {
            success = true;
            break;
        }
    } while (++step < cnst::max_iter);

    if (success) {
        cout << "Iteration converged";
    } else {
        cout << "Iteration did not converge";
    }
    cout << ": res = " << res << ", |Ax - b| = " << diff << ", step = " << step
         << endl;

    return x;
}

Vec psolve(const CsrMatrix& lhs, const Vec& rhs, Vec x0)
{
    auto [m, n] = lhs.shape();

    check_if(!any_of(begin(lhs.diag()), end(lhs.diag()),
                     [](auto& x) { return iszero(x); }),
             "Zero on diag");
    check_if(rhs.size() == n, "Incompatible shapes");

    size_t step = 0;

    if (x0.empty()) {
        x0.resize(m, 0);
    }
    auto rnorm = cnorm(rhs);
    auto r = rhs - lhs * x0;
    auto r0 = r;
    auto x = move(x0);

    Vec v(m, 0.);
    Vec s(m, 0.);
    Vec t(m, 0.);
    Vec y(m, 0.);
    Vec z(m, 0.);
    Vec q(m, 0.);
    Vec p(m, 0.);

    Vec ax(m, 0.);

    bool success = false;

    double diff = 0, res = 0, u = 1, a = 1, w = 1, b = 0;
    do {
        auto up = u;
        u = dot(r0, r);

        b = (u / up) * (a / w);
        // p = r + b * (p - w * v)
        cax_by_z(p, b, p, -w * b, v, r);
        // y = K^{-1} p
        dot_diag(y, lhs, p);
        // v = lhs * y
        dot(v, lhs, y);

        a = u / dot(r0, v);

        // s = r - a * v
        cax_y(s, -a, v, r);
        // z = K^{-1} s
        dot_diag(z, lhs, s);
        // t = lhs * z
        dot(t, lhs, z);
        // q = K^{-1} * t
        dot_diag(q, lhs, t);
        w = dot(q, z) / sqr(q);

        // r = s - w * t
        cax_y(r, -w, t, s);
        res = cnorm(r) / rnorm;
        if (!isnormal(res)) {
            break;
        }

        // x = a * y + w * z + x
        cax_by_z(x, a, y, w, z, x);

        dot(ax, lhs, x);
        diff = cdist(ax, rhs);
        if (iszero(diff)) {
            success = true;
            break;
        }
    } while (++step < cnst::max_iter);

    if (success) {
        cout << "Iteration converged";
    } else {
        cout << "Iteration did not converge";
    }
    cout << ": res = " << res << ", |Ax - b| = " << diff << ", step = " << step
         << endl;

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

void dot_diag(Vec& res, const CsrMatrix& lhs, const Vec& rhs)
{
    CsrMatrix::Index i, n = res.size();

    for (i = 0; i < n; ++i) {
        res[i] = rhs[i] / lhs.diag()[i];
    }
}