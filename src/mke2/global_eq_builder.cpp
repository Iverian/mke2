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
                            }
                        }
                    }
                }
            }
        }
    }

    return {mat, vec};
}

// if (!t.on_first(k[i])) {
//     for (Index j = 0; j < Triangulation::N; ++j) {
//         auto gj = t.get_index(k[j]);

//         if (!t.on_first(k[j])) {
//             for (Index p = 0; p < Triangulation::DIM; ++p) {
//                 for (Index q = 0; q < Triangulation::DIM; ++q) {
//                     auto v = r.first(_v(i, p), _v(j, q));
//                     mat.add(_g(gi, p, m), _g(gj, q, m), v);
//                 }
//             }
//         } else {
//             for (Index p = 0; p < Triangulation::DIM; ++p) {
//                 for (Index q = 0; q < Triangulation::DIM; ++q) {
//                     vec[_g(gi, p, m)]
//                         -= r.first(_v(i, p), _v(j, q))
//                         * vec[_g(gj, q, m)];
//                 }
//             }
//         }
//     }
// }

// GlobalEqBuilder::GlobalEqBuilder(const Triangulation& triang,
//                                  shared_ptr<AbstractLocalEq> gen)
//     : t_(triang)
//     , g_(gen)
//     , m_(Index(t_.nodes().size()))
//     , lhs_(3 * m_)
//     , rhs_(3 * m_)
// {
// }

// GlobalEqBuilder& GlobalEqBuilder::get()
// {
//     // сборка правой части глобальной СЛАУ
//     for (auto& k : t_.elems()) {
//         auto v = g_->get_internal(t_.data(k)).second;
//         for (auto i = 0; i < 4; ++i) {
//             auto p = k[i]->second;
//             for (auto p = 0; p < 3; ++p) {
//                 rhs_[p * m_ + p] += v[p * Triangulation::N + i];
//             }
//         }
//     }

//     // ГУ третьего рода
//     for (auto& s : t_.third()) {
//         auto r = g_->get_boundary(t_.data(s));
//         for (auto i = 0; i < 3; ++i) {
//             auto p = s[i]->second;
//             for (auto p = 0; p < 3; ++p) {
//                 rhs_[p * m_ + p] += r.second[p * Triangulation::SN + i];
//             }

//             for (auto j = 0; j < 3; ++j) {
//                 auto q = s[j]->second;
//                 for (auto p = 0; p < 3; ++p) {
//                     lhs_.add(p * m_ + p, p * m_ + q,
//                                    r.first[p * Triangulation::SN + i,
//                                            p * Triangulation::SN + j]);
//                 }
//             }
//         }
//     }

//     // ГУ первого рода Sigma_1 (x, z)
//     for (auto& n : t_.first()[0]) {
//         for (auto i = 0; i < 3; ++i) {
//             auto p = n->second;
//             rhs_[p] = 0;
//             rhs_[2 * m_ + p] = 0;
//             for (auto p = 0; p < 3; ++p) {
//                 lhs_.set(p * m_ + p, p * m_ + p, 1);
//             }
//         }
//     }

//     // ГУ первого рода Sigma_2 (y)
//     for (auto& n : t_.first()[1]) {
//         for (auto i = 0; i < 3; ++i) {
//             auto p = n->second;
//             rhs_[m_ + p] = 0;
//             for (auto p = 0; p < 3; ++p) {
//                 lhs_.set(p * m_ + p, p * m_ + p, 1);
//             }
//         }
//     }

//     // Сборка матрицы глобальной СЛАУ
//     for (auto& k : t_.elems()) {
//         auto r = g_->get_internal(t_.data(k));
//         for (auto i = 0; i < 4; ++i) {
//             auto p = k[i]->second;

//             if (!t_.on_first(k[i])) {
//                 for (auto j = 0; j < 4; ++j) {
//                     auto q = k[j]->second;
//                     if (!t_.on_first(k[j])) {
//                         for (auto p = 0; p < 3; ++p) {
//                             lhs_.add(p * m_ + p, p * m_ + q,
//                                            r.first(p * Triangulation::N + i,
//                                                    p * Triangulation::N +
//                                                    j));
//                         }
//                     } else {
//                         for (auto p = 0; p < 3; ++p) {
//                             rhs_[p * m_ + p]
//                                 -= r.first(p * Triangulation::N + i,
//                                            p * Triangulation::N + j)
//                                 * rhs_[p * m_ + q];
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     lhs_.remove_zeroes();
//     return *this;
// }

// CsrMatrix& GlobalEqBuilder::mat()
// {
//     return lhs_;
// }

// Vec& GlobalEqBuilder::vec()
// {
//     return rhs_;
// }
