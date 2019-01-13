#include "local_eq_v17.hpp"
#include "constants.hpp"

#include <debug.hpp>
#include <lup_factor.hpp>
#include <util.hpp>

#include <iterator>

using namespace std;

const DenseMatrix LocalEqV17::cm({6, 6}, begin(cnst::C), end(cnst::C));

struct InternalProxy {
    explicit InternalProxy(Triangulation::FiniteElementData& pp)
        : p(pp)
        , s()
        , v()
    {
        v = triple(p[0] - p[3], p[1] - p[3], p[2] - p[3]);
        s = get_s_mat();
    }

    DenseMatrix get_dq_mat() const
    {
        DenseMatrix mat({4, 4}, 1);
        for (auto i = 0; i < 3; ++i) {
            for (auto j = 0; j < 4; ++j) {
                mat(i, j) = p[j][i];
            }
        }
        auto r = LupFactor(mat).factor();

        array<Vec, 3> f = {r.solve({1, 0, 0, 0}), r.solve({0, 1, 0, 0}),
                           r.solve({0, 0, 1, 0})};

        DenseMatrix result({6, 12});
        for (auto i = 0; i < 4; ++i) {
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

    DenseMatrix get_gk_mat() const
    {
        auto dq = get_dq_mat();
        auto dqt = dq.transpose();

        return dqt * LocalEqV17::cm * dq;
    }

    DenseMatrix get_s_mat() const
    {
        DenseMatrix result({4, 4}, 1);

        for (DenseMatrix::Index i = 0; i < 4; ++i) {
            result(i, i) = 2;
        }
        result *= (cnst::rho / 120);

        result = kroneker_product(result, DenseMatrix::eye(3));

        return result;
    }

    LocalEqV17::Result get()
    {
        Vec g(12, 0.);
        for (auto i = 0; i < 4; ++i) {
            g[3 * i + 2] = -cnst::p;
        }

        return {v * (get_gk_mat() / 6 - sqr(cnst::omega) * s), v * s * g};
    }

private:
    Triangulation::FiniteElementData& p;
    DenseMatrix s;
    double v;
};

struct SurfaceProxy {

    explicit SurfaceProxy(Triangulation::SurfaceElementData& pp)
        : p(pp)
        , s()
        , v()
    {
        v = norm(cross(p[0] - p[2], p[1] - p[2]));
        s = get_s_mat();
    }

    DenseMatrix get_s_mat()
    {
        DenseMatrix result({3, 3}, 1);
        for (auto i = 0; i < 3; ++i) {
            result(i, i) = 2;
        }
        result *= 5 * cnst::mu / 120;
        result = kroneker_product(result, DenseMatrix::eye(3));
        return result;
    }

    LocalEqV17::Result get()
    {
        return {v * s, Vec(9, 0.)};
    }

private:
    Triangulation::SurfaceElementData& p;
    DenseMatrix s;
    double v;
};

LocalEqV17::Result
LocalEqV17::get_internal(Triangulation::FiniteElementData elem) const
{
    return InternalProxy(elem).get();
}

LocalEqV17::Result
LocalEqV17::get_boundary(Triangulation::SurfaceElementData elem) const
{
    return SurfaceProxy(elem).get();
}
