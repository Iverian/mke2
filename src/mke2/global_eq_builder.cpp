#include <global_eq_builder.hpp>

#include "global_indices.hpp"

using namespace std;

using Index = Triangulation::Index;

pair<SparceMatrix, Vec> build_global_system(const Triangulation& t,
                                            AbstractLocalEq& gen)
{
    auto m = Index(t.nodes().size());

    SparceMatrix mat(Triangulation::DIM * m);
    Vec vec(Triangulation::DIM * m, 0.);

    // Сборка правой части СЛАУ
    for (auto& k : t.elems()) {
        auto b = gen.get_internal(t.data(k)).second;

        for (Index i = 0; i < Triangulation::N; ++i) {
            auto gi = t.get_index(k[i]);

            for (Index p = 0; p < Triangulation::DIM; ++p) {
                vec[_g(gi, p, m)] += b[_v(i, p)];
            }
        }
    }

    // Учет ГУ 3го рода
    for (auto& s : t.third()) {
        auto r = gen.get_boundary(t.data(s));

        for (Index i = 0; i < Triangulation::SN; ++i) {
            auto gi = t.get_index(s[i]);

            for (Index j = 0; j < Triangulation::SN; ++j) {
                auto gj = t.get_index(s[j]);

                for (Index p = 0; p < Triangulation::DIM; ++p) {
                    for (Index q = 0; q < Triangulation::DIM; ++q) {
                        auto v = r.first(_s(i, p), _s(j, q));
                        mat.add(_g(gi, p, m), _g(gj, q, m), v);
                    }
                }
            }

            for (Index p = 0; p < Triangulation::DIM; ++p) {
                vec[_g(gi, p, m)] += r.second[_s(i, p)];
            }
        }
    }
    // Учет ГУ 1го рода
    // Sigma_1 u_x = u_z = 0
    for (auto& n : t.first()[0]) {
        auto gn = t.get_index(n);

        for (Index p : {Coord::X, Coord::Z}) {
            vec[_g(gn, p, m)] = 0;
            for (Index q : {Coord::X, Coord::Z}) {
                mat.set(_g(gn, p, m), _g(gn, q, m), 1);
            }
        }
    }
    // Sigma_2 u_y = 0
    for (auto& n : t.first()[1]) {
        auto gn = t.get_index(n);

        vec[_g(gn, Coord::Y, m)] = 0;
        mat.set(_g(gn, Coord::Y, m), _g(gn, Coord::Y, m), 1);
    }

    // Основной цикл сборки СЛАУ
    for (auto& k : t.elems()) {
        auto r = gen.get_internal(t.data(k));

        for (Index i = 0; i < Triangulation::N; ++i) {
            auto gi = t.get_index(k[i]);

            if (!t.on_first(k[i])) {
                for (Index j = 0; j < Triangulation::N; ++j) {
                    auto gj = t.get_index(k[j]);

                    if (!t.on_first(k[j])) {
                        for (Index p = 0; p < Triangulation::DIM; ++p) {
                            for (Index q = 0; q < Triangulation::DIM; ++q) {
                                auto v = r.first(_v(i, p), _v(j, q));
                                mat.add(_g(gi, p, m), _g(gj, q, m), v);
                            }
                        }
                    } else {
                        for (Index p = 0; p < Triangulation::DIM; ++p) {
                            for (Index q = 0; q < Triangulation::DIM; ++q) {
                                vec[_g(gi, p, m)]
                                    -= r.first(_v(i, p), _v(j, q))
                                    * vec[_g(gj, q, m)];
                            }
                        }
                    }
                }
            }
        }
    }

    mat.remove_zeroes();
    return {mat, vec};
}

pair<SparceMatrix, Vec> build_global_system(const Triangulation& t,
                                            LocalEqGen gen)
{
    auto m = Index(t.nodes().size());
    SparceMatrix mat(Triangulation::DIM * m);
    Vec vec(Triangulation::DIM * m, 0);

    for (auto& k : t.elems()) {
        auto b = gen(t, k).second;

        for (Index i = 0; i < Triangulation::N; ++i) {
            auto gi = t.get_index(k[i]);

            for (Index p = 0; p < Triangulation::DIM; ++p) {
                vec[_g(gi, p, m)] += b[_v(i, p)];
            }
        }
    }

    for (auto& n : t.first()[0]) {
        auto gn = t.get_index(n);

        for (Index p : {Coord::X, Coord::Z}) {
            vec[_g(gn, p, m)] = 0;
            for (Index q : {Coord::X, Coord::Z}) {
                mat.set(_g(gn, p, m), _g(gn, q, m), 1);
            }
        }
    }

    for (auto& n : t.first()[1]) {
        auto gn = t.get_index(n);

        vec[_g(gn, Coord::Y, m)] = 0;
        mat.set(_g(gn, Coord::Y, m), _g(gn, Coord::Y, m), 1);
    }

    for (auto& k : t.elems()) {
        auto r = gen(t, k);

        for (Index i = 0; i < Triangulation::N; ++i) {
            auto gi = t.get_index(k[i]);

            if (!t.on_first(k[i])) {
                for (Index j = 0; j < Triangulation::N; ++j) {
                    auto gj = t.get_index(k[j]);

                    if (!t.on_first(k[j])) {
                        for (Index p = 0; p < Triangulation::DIM; ++p) {
                            for (Index q = 0; q < Triangulation::DIM; ++q) {
                                mat.add(_g(gi, p, m), _g(gj, q, m),
                                        r.first(_v(i, p), _v(j, q)));
                            }
                        }
                    } else {
                        for (Index p = 0; p < Triangulation::DIM; ++p) {
                            for (Index q = 0; q < Triangulation::DIM; ++q) {
                                vec[_g(gi, p, m)]
                                    -= r.first(_v(i, p), _v(j, q))
                                    * vec[_g(gj, q, m)];
                            }
                        }
                    }
                }
            }
        }
    }

    mat.remove_zeroes();
    return {mat, vec};
}

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

// SparceMatrix& GlobalEqBuilder::mat()
// {
//     return lhs_;
// }

// Vec& GlobalEqBuilder::vec()
// {
//     return rhs_;
// }
