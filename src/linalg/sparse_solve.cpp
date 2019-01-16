#include <debug.hpp>
#include <sparse_matrix.hpp>

#include <cmath>

static constexpr size_t max_iter = 100000;

using namespace std;

void cax_y(Vec& res, double a, const Vec& x, const Vec& y);
void cax_by_z(Vec& res, double a, const Vec& x, double b, const Vec& y,
              const Vec& z);
double cdist(const Vec& lhs, const Vec& rhs);
void mdot_diag(Vec& result, const SparseMatrix& lhs, const Vec& rhs);

Vec solve(const SparseMatrix& lhs, const Vec& rhs, Vec x0)
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

    double res = 0, u = 1, a = 1, w = 1, b = 0;
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

        // x = a * p + w * s + x
        cax_by_z(x, a, p, w, s, x);
        // r = s - w * t
        cax_y(r, -w, t, s);

        res = sqrt(sqr(r) / rsqr);
        if (isnear(res, 0, Tolerance::TRIPLE)) {
            break;
        }
    } while (++step < max_iter);

    auto diff = sqrt(sqr(lhs * x - rhs));
    if (step != max_iter) {
        cout << "Iteration converged";
    } else {
        cout << "Iteration did not converge";
    }
    cout << ": res=" << res << ", |Ax - b|=" << diff << ", step=" << step
         << endl;

    return x;
}

Vec solve_p(const SparseMatrix& lhs, const Vec& rhs, Vec x0)
{
    auto [m, n] = lhs.shape();
    check_if(rhs.size() == n, "Incompatible shapes");
    check_if(!any_of(begin(lhs.diag_), end(lhs.diag_),
                     [](auto& x) { return isnear(x, 0); }),
             "Zero on diag");

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

    mdot_diag(z, lhs, r);
    auto p = z;

    double res = 0, u = 1, a = 1, w = 1, b = 0;
    do {
        auto up = u;
        u = dot(r0, r);

        auto b = (u / up) * (a / w);
        // p = r + b * (p - w * v)
        cax_by_z(p, b, p, -w * b, v, r);
        // y = diag(lhs)^{-1} * p
        mdot_diag(y, lhs, p);
        // v = lhs * y
        dot(v, lhs, y);

        a = u / dot(r0, v);

        // s = r - a * v
        cax_y(s, -a, v, r);
        // z = diag(lhs)^{-1} * s
        mdot_diag(z, lhs, s);
        // t = lhs * s
        dot(t, lhs, z);
        // q = diag(lhs)^{-1} * t
        mdot_diag(q, lhs, t);
        w = dot(z, q) / sqr(q);

        // x = a * p + w * s + x
        cax_by_z(x, a, y, w, z, x);
        // r = s - w * t
        cax_y(r, -w, t, s);

        res = sqrt(sqr(r) / rsqr);
        if (isnear(res, 0, Tolerance::DOUBLE)) {
            break;
        }
    } while (++step < max_iter);

    res = sqrt(sqr(lhs * x - rhs));
    if (isnear(res, 0, Tolerance::SINGLE)) {
        cout << "Iteration converged: res = " << res << ", step = " << step
             << endl;
    } else {
        cout << "Iteration did not converge: res = " << res
             << ", step = " << step << endl;
    }

    return x;
}

void mdot_diag(Vec& result, const SparseMatrix& lhs, const Vec& rhs)
{
    for (SparseMatrix::Index i = 0; i < lhs.shape().m; ++i) {
        result[i] = rhs[i] / lhs.diag()[i];
    }
}

void cax_y(Vec& res, double a, const Vec& x, const Vec& y)
{
    SparseMatrix::Index i, n = res.size();

    for (i = 0; i < n; ++i) {
        res[i] = a * x[i] + y[i];
    }
}

void cax_by_z(Vec& res, double a, const Vec& x, double b, const Vec& y,
              const Vec& z)
{
    SparseMatrix::Index i, n = res.size();

    for (i = 0; i < n; ++i) {
        res[i] = a * x[i] + b * y[i] + z[i];
    }
}
