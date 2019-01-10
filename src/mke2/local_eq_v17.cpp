#include "local_eq_v17.hpp"
#include "constants.hpp"

#include <lup_factor.hpp>
#include <util.hpp>

#include <iterator>

using namespace std;

DenseMatrix LocalEqV17::cm({6, 6}, begin(cnst::C), end(cnst::C));

LocalEqV17::LocalEqV17(const Triang3d::Node& node)
    : p_(node)
    , v_()
    , s0_()
{
    v_ = triple(p_.v(0) - p_.v(3), p_.v(1) - p_.v(3), p_.v(2) - p_.v(3));
    s0_ = get_s0_mat();
}

DenseMatrix LocalEqV17::get_internal_mat() const
{
    return v_ * (get_gk_mat() / 6 - sqr(cnst::omega) * s0_);
}

DenseMatrix LocalEqV17::get_boundary_mat() const
{
    return kroneker_product((5 * cnst::mu / 120) * get_bskeleton(),
                            DenseMatrix::eye(3));
}

Vec LocalEqV17::get_internal_vec() const
{
    Vec g(12, 0.);
    for (size_t i = 0; i < 4; ++i) {
        g[3 * i + 2] = -cnst::p;
    }
    return v_ * s0_ * g;
}

Vec LocalEqV17::get_boundary_vec() const
{
    return Vec(12, 0.);
}

DenseMatrix LocalEqV17::get_dq_mat() const
{
    auto r = LupFactor(p_.vertex_view().append_row(Vec(4, 1.))).factor();
    array<Vec, 3> f
        = {r.solve({1, 0, 0}), r.solve({0, 1, 0}), r.solve({0, 0, 1})};

    DenseMatrix result({6, 12});
    for (size_t i = 0; i < 4; ++i) {
        auto j = 3 * i;
        result(0, j + 0) = f[0][i];
        result(1, j + 1) = f[1][i];
        result(2, j + 2) = f[2][i];
        result(3, j + 0) = f[1][i];
        result(3, j + 1) = f[0][i];
        result(4, j + 0) = f[2][i];
        result(4, j + 2) = f[0][i];
        result(5, j + 1) = f[2][i];
        result(5, j + 2) = f[1][i];
    }

    return result;
}

DenseMatrix LocalEqV17::get_gk_mat() const
{
    auto dq = get_dq_mat();
    auto dqt = dq.transpose();

    return dq * cm * dqt;
}

DenseMatrix LocalEqV17::get_s0_mat() const
{
    DenseMatrix result({4, 4}, 1);

    for (DenseMatrix::Index i = 0; i < 4; ++i) {
        result(i, i) = 2;
    }
    result *= (cnst::rho / 120);
    result = kroneker_product(result, DenseMatrix::eye(3));

    return result;
}

DenseMatrix LocalEqV17::get_bskeleton() const
{
    DenseMatrix result({4, 4});
    for (auto& i : p_.boundary()) {
        if (p_.normal(i) == Point3d(0, 0, 1)) {
            result += get_face_bmat(i);
        }
    }
    return result;
}

DenseMatrix LocalEqV17::get_face_bmat(size_t index) const
{
    DenseMatrix result({4, 4});
    auto m = (index + 1) % 4;
    auto n = (index + 2) % 4;
    auto p = (index + 3) % 4;
    auto s = norm(cross(p_.v(m) - p_.v(p), p_.v(n) - p_.v(p)));

    result(m, m) = result(n, n) = result(p, p) = 2;
    result(m, n) = result(m, p) = result(n, p) = 1;
    result *= s;

    return result;
}