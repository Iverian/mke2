#include "local_eq_v17.hpp"
#include "constants.hpp"

#include <debug.hpp>
#include <lup_factor.hpp>
#include <util.hpp>

#include <iterator>
#include <tuple>

#include "global_indices.hpp"

using namespace std;

using Index = Triangulation::Index;

class V17InternalGen {
    const Triangulation::FiniteElement::Data data;

public:
    explicit V17InternalGen(const Triangulation::FiniteElement::Data& data);

    double vk() const;
    DenseMatrix gk() const;
    DenseMatrix internal(double val) const;

private:
    DenseMatrix dq() const;
};

class V17SurfaceGen {
    const Triangulation::SurfaceElement::Data data;

public:
    explicit V17SurfaceGen(const Triangulation::SurfaceElement::Data& data);

    DenseMatrix boundary(double val) const;
};

class V17Gen {
    const Triangulation& t;
    const Triangulation::FiniteElement elem;
    const V17InternalGen igen;

public:
    V17Gen(const Triangulation& t, const Triangulation::FiniteElement& elem);

    double vk() const;
    DenseMatrix gk() const;
    DenseMatrix internal(double val) const;
    DenseMatrix boundary(double val) const;
    DenseMatrix boundary(const Triangulation::SurfaceElement& face, Index i,
                         double val) const;
};

LocalEqV17::Result
LocalEqV17::get_internal(Triangulation::FiniteElement::Data elem) const
{
    V17InternalGen g(elem);

    auto vk = g.vk();
    auto e3 = DenseMatrix::eye(Triangulation::DIM);
    auto s0 = kroneker_product(g.internal(vk * cnst::rho / 120), e3);

    auto mat = (vk / 6) * g.gk() - sqr(cnst::omega) * s0;
    auto vec = Vec(Triangulation::DIM * Triangulation::N, 0.);

    return {mat, vec};
}

LocalEqV17::Result
LocalEqV17::get_boundary(Triangulation::SurfaceElement::Data elem) const
{
    V17SurfaceGen g(elem);

    Vec tk(Triangulation::DIM * Triangulation::SN, 0.);
    for (Index i = 0; i < Triangulation::SN; ++i) {
        tk[_s(i, Coord::Z)] = -cnst::p;
    }

    auto e3 = DenseMatrix::eye(Triangulation::DIM);
    auto mat = kroneker_product(g.boundary(1. / 24), e3);
    auto vec = mat * tk;

    return {mat * cnst::mu, vec};
}

V17InternalGen::V17InternalGen(const Triangulation::FiniteElement::Data& data_)
    : data(data_)
{
}

double V17InternalGen::vk() const
{
    return triple(data[0] - data[3], data[1] - data[3], data[2] - data[3]);
}

DenseMatrix V17InternalGen::dq() const
{
    DenseMatrix mat({Triangulation::N, Triangulation::N}, 1);
    for (Index i = 0; i < Triangulation::DIM; ++i) {
        for (Index j = 0; j < Triangulation::N; ++j) {
            mat(i, j) = data[j][i];
        }
    }
    auto r = LupFactor(mat).factor();

    array<Vec, Triangulation::DIM> f = {
        r.solve({1, 0, 0, 0}), r.solve({0, 1, 0, 0}), r.solve({0, 0, 1, 0})};

    DenseMatrix result(
        {2 * Triangulation::DIM, Triangulation::N * Triangulation::DIM});
    for (Index i = 0; i < Triangulation::N; ++i) {
        auto j = Triangulation::DIM * i;
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

DenseMatrix V17InternalGen::gk() const
{
    static DenseMatrix cm({2 * Triangulation::DIM, 2 * Triangulation::DIM},
                          begin(cnst::C), end(cnst::C));

    auto m = dq();
    auto mt = m.transpose();
    return mt * cm * m;
}

DenseMatrix V17InternalGen::internal(double val) const
{
    DenseMatrix result({Triangulation::N, Triangulation::N}, val);
    for (Index i = 0; i < Triangulation::N; ++i) {
        result(i, i) *= 2;
    }
    return result;
}

V17SurfaceGen::V17SurfaceGen(const Triangulation::SurfaceElement::Data& data_)
    : data(data_)
{
}

DenseMatrix V17SurfaceGen::boundary(double val) const
{
    auto s = norm(cross(data[0] - data[2], data[1] - data[2]));

    DenseMatrix result({Triangulation::SN, Triangulation::SN}, val * s);
    for (Index i = 0; i < Triangulation::SN; ++i) {
        result(i, i) *= 2;
    }
    return result;
}

// Генератор, совмещающий матрицу по границе и по объему
LocalEqGen::result_type v17(const Triangulation& t,
                            const Triangulation::FiniteElement& elem)
{
    V17Gen gen(t, elem);

    auto vk = gen.vk();
    auto e3 = DenseMatrix::eye(Triangulation::DIM);
    auto s0 = kroneker_product(gen.internal(vk * cnst::rho / 120), e3);
    auto sb = kroneker_product(gen.boundary(1. / 24), e3);

    Vec tk(Triangulation::DIM * Triangulation::N, 0.);
    for (Index i = 0; i < Triangulation::N; ++i) {
        tk[_v(i, Coord::Z)] = -cnst::p;
    }

    auto mat = vk * gen.gk() / 6 - sqr(cnst::omega) * s0 + cnst::mu * sb;
    auto vec = sb * tk;

    return make_pair(mat, vec);
}

V17Gen::V17Gen(const Triangulation& t_,
               const Triangulation::FiniteElement& elem_)
    : t(t_)
    , elem(elem_)
    , igen(elem_.data())
{
}

double V17Gen::vk() const
{
    return igen.vk();
}

DenseMatrix V17Gen::internal(double val) const
{
    return igen.internal(val);
}

DenseMatrix V17Gen::boundary(double val) const
{
    DenseMatrix result({4, 4});
    for (Index i = 0; i < Triangulation::N; ++i) {
        auto f = elem.face(i);
        if (t.on_third(f)) {
            result += boundary(f, i, val);
        }
    }
    return result;
}

DenseMatrix V17Gen::boundary(const Triangulation::SurfaceElement& face,
                             Index i, double val) const
{
    auto data = face.data();

    DenseMatrix result({4, 4}, 0.);
    auto m = (i + 1) % Triangulation::N;
    auto n = (i + 2) % Triangulation::N;
    auto p = (i + 3) % Triangulation::N;
    auto c = norm(cross(data[0] - data[2], data[1] - data[2]));

    result(m, m) = result(n, n) = result(p, p) = 2 * c * val;
    result(m, n) = result(m, p) = result(n, p) = result(n, m) = result(p, m)
        = result(p, n) = 1 * c * val;

    return result;
}

DenseMatrix V17Gen::gk() const
{
    return igen.gk();
}

// const DenseMatrix LocalEqV17::cm({6, 6}, begin(cnst::C), end(cnst::C));

// struct InternalProxy {
//     explicit InternalProxy(Triangulation::FiniteElementData& pp)
//         : p(pp)
//         , s()
//         , v()
//     {
//         v = triple(p[0] - p[3], p[1] - p[3], p[2] - p[3]);
//         s = get_s_mat();
//     }

//     DenseMatrix get_dq_mat() const
//     {
// DenseMatrix mat({4, 4}, 1);
// for (auto i = 0; i < 3; ++i) {
//     for (auto j = 0; j < 4; ++j) {
//         mat(i, j) = p[j][i];
//     }
// }
// auto r = LupFactor(mat).factor();

// array<Vec, 3> f = {r.solve({1, 0, 0, 0}), r.solve({0, 1, 0, 0}),
//                    r.solve({0, 0, 1, 0})};

// DenseMatrix result({6, 12});
// for (auto i = 0; i < 4; ++i) {
//     auto j = 3 * i;
//     result(0, j + 0) = f[0][i];
//     result(1, j + 1) = f[1][i];
//     result(2, j + 2) = f[2][i];
//     result(3, j + 0) = f[1][i];
//     result(3, j + 1) = f[0][i];
//     result(4, j + 0) = f[2][i];
//     result(4, j + 2) = f[0][i];
//     result(5, j + 1) = f[2][i];
//     result(5, j + 2) = f[1][i];
// }

// return result;
//     }

//     DenseMatrix get_gk_mat() const
//     {
//         auto dq = get_dq_mat();
//         auto dqt = dq.transpose();

//         return dqt * LocalEqV17::cm * dq;
//     }

//     DenseMatrix get_s_mat() const
//     {
//         DenseMatrix result({4, 4}, 1);

//         for (DenseMatrix::Index i = 0; i < 4; ++i) {
//             result(i, i) = 2;
//         }
//         result *= (cnst::rho / 120);

//         result = kroneker_product(result, DenseMatrix::eye(3));

//         return result;
//     }

//     LocalEqV17::Result get()
//     {
//         return {v * (get_gk_mat() / 6 - sqr(cnst::omega) * s), Vec(12, 0.)};
//     }

// private:
//     Triangulation::FiniteElementData& p;
//     DenseMatrix s;
//     double v;
// };

// struct SurfaceProxy {

//     explicit SurfaceProxy(Triangulation::SurfaceElementData& pp)
//         : p(pp)
//         , s()
//         , v()
//     {
//         v = norm(cross(p[0] - p[2], p[1] - p[2]));
//         s = get_s_mat();
//     }

//     DenseMatrix get_s_mat()
//     {
//         DenseMatrix result({3, 3}, 1);
//         for (auto i = 0; i < 3; ++i) {
//             result(i, i) = 2;
//         }
//         // result *= 5 * cnst::mu / 120;
//         result = kroneker_product(result, DenseMatrix::eye(3));
//         return result;
//     }

//     LocalEqV17::Result get()
//     {
//         Vec t(9, 0.);
//         for (auto i = 0; i < 3; ++i) {
//             t[i * 3 + 2] = -cnst::p;
//         }

//         return {v * s * (5 * cnst::mu / 120), v * s  * t};
//     }

// private:
//     Triangulation::SurfaceElementData& p;
//     DenseMatrix s;
//     double v;
// };

// LocalEqV17::Result
// LocalEqV17::get_internal(Triangulation::FiniteElementData elem) const
// {
//     return InternalProxy(elem).get();
// }

// LocalEqV17::Result
// LocalEqV17::get_boundary(Triangulation::SurfaceElementData elem) const
// {
//     return SurfaceProxy(elem).get();
// }
