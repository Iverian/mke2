#include <dok_matrix.hpp>
#include <global_eq_builder.hpp>
#include <global_indices.hpp>

using namespace std;

using Index = Triangulation::Index;

pair<CsrMatrix, Vec> build_global_system(const Triangulation& t,
                                         LocalEqGen gen)
{
    auto m = Index(t.nodes().size());
    DokMatrix mat(Triangulation::DIM * m);
    Vec vec(Triangulation::DIM * m, 0);

    // Учет ГУ 1го рода
    for (auto& fc : t.first()) {
        for (auto& n : fc) {
            auto gn = n.index();
            auto cn = t.on_first(n);

            for (Index p = 0; p < Triangulation::DIM; ++p) {
                auto cnp = t.on_first(cn, p);

                if (cnp) {
                    auto gnp = _g(gn, p, m);

                    vec[gnp] = 0.;
                    mat.set(gnp, gnp, 1.);
                }
            }
        }
    }

    // Основная часть формирования СЛАУ
    for (auto& k : t.elems()) {
        const auto [kmat, kvec] = gen(t, k);

        for (Index i = 0; i < Triangulation::N; ++i) {
            const auto gi = k[i].index();
            const auto ci = t.on_first(k[i]);

            for (Index p = 0; p < Triangulation::DIM; ++p) {
                const auto cip = t.on_first(ci, p);

                if (!cip) {
                    vec[_g(gi, p, m)] += kvec[_v(i, p)];
                }
            }

            for (Index j = 0; j < Triangulation::N; ++j) {
                const auto gj = k[j].index();
                const auto cj = t.on_first(k[j]);

                for (Index p = 0; p < Triangulation::DIM; ++p) {
                    const auto cip = t.on_first(ci, p);

                    if (!cip) {
                        const auto gip = _g(gi, p, m);

                        for (Index q = 0; q < Triangulation::DIM; ++q) {
                            const auto gjq = _g(gj, q, m);
                            const auto val = kmat(_v(i, p), _v(j, q));

                            const auto cjq = t.on_first(cj, q);
                            if (!cjq) {
                                mat.add(gip, gjq, val);
                            } else {
                                // vec[gip] -= val * vec[gjq];
                            }
                        }
                    }
                }
            }
        }
    }

    return {CsrMatrix(mat), vec};
}

pair<DenseMatrix, Vec> build_global_system_dense(const Triangulation& t,
                                                 LocalEqGen gen)
{
    auto m = Index(t.nodes().size());
    DenseMatrix mat({Triangulation::DIM * m, Triangulation::DIM * m}, 0.);
    Vec vec(Triangulation::DIM * m, 0.);

    // Учет ГУ 1го рода
    for (auto& fc : t.first()) {
        for (auto& n : fc) {
            auto gn = n.index();
            auto cn = t.on_first(n);

            for (Index p = 0; p < Triangulation::DIM; ++p) {
                auto cnp = t.on_first(cn, p);

                if (cnp) {
                    auto gnp = _g(gn, p, m);

                    vec[gnp] = 0.;
                    mat(gnp, gnp) = 1.;
                }
            }
        }
    }

    // Основная часть формирования СЛАУ
    for (auto& k : t.elems()) {
        const auto [kmat, kvec] = gen(t, k);

        for (Index i = 0; i < Triangulation::N; ++i) {
            const auto gi = k[i].index();
            const auto ci = t.on_first(k[i]);

            for (Index p = 0; p < Triangulation::DIM; ++p) {
                const auto cip = t.on_first(ci, p);

                if (!cip) {
                    vec[_g(gi, p, m)] += kvec[_v(i, p)];
                }
            }

            for (Index j = 0; j < Triangulation::N; ++j) {
                const auto gj = k[j].index();
                const auto cj = t.on_first(k[j]);

                for (Index p = 0; p < Triangulation::DIM; ++p) {
                    const auto cip = t.on_first(ci, p);

                    if (!cip) {
                        const auto gip = _g(gi, p, m);

                        for (Index q = 0; q < Triangulation::DIM; ++q) {
                            const auto gjq = _g(gj, q, m);
                            const auto val = kmat(_v(i, p), _v(j, q));

                            const auto cjq = t.on_first(cj, q);
                            if (!cjq) {
                                mat(gip, gjq) += val;
                            } else {
                                // vec[gip] -= val * vec[gjq];
                            }
                        }
                    }
                }
            }
        }
    }

    return {mat, vec};
}
