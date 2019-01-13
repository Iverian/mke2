#include <global_eq_builder.hpp>

using namespace std;

GlobalEqBuilder::GlobalEqBuilder(const Triangulation& triang,
                                 shared_ptr<AbstractLocalEq> gen)
    : t_(triang)
    , g_(gen)
    , m_(Index(t_.nodes().size()))
    , lhs_(3 * m_)
    , rhs_(3 * m_)
{
}

GlobalEqBuilder& GlobalEqBuilder::get()
{
    // сборка правой части глобальной СЛАУ
    for (auto& k : t_.elems()) {
        auto v = g_->get_internal(t_.data(k)).second;
        for (auto i = 0; i < 4; ++i) {
            auto p = k[i]->second;
            for (auto l = 0; l < 3; ++l) {
                rhs_[l * m_ + p] += v[l * Triangulation::N + i];
            }
        }
    }

    // ГУ третьего рода
    for (auto& s : t_.third()) {
        auto r = g_->get_boundary(t_.data(s));
        for (auto i = 0; i < 3; ++i) {
            auto p = s[i]->second;
            for (auto l = 0; l < 3; ++l) {
                rhs_[l * m_ + p] += r.second[l * Triangulation::SN + i];
            }

            for (auto j = 0; j < 3; ++j) {
                auto q = s[j]->second;
                for (auto l = 0; l < 3; ++l) {
                    lhs_.fetch_add(l * m_ + p, l * m_ + q,
                                   r.first[l * Triangulation::SN + i,
                                           l * Triangulation::SN + j]);
                }
            }
        }
    }

    // ГУ первого рода Sigma_1 (x, z)
    for (auto& n : t_.first()[0]) {
        for (auto i = 0; i < 3; ++i) {
            auto p = n->second;
            rhs_[p] = 0;
            rhs_[2 * m_ + p] = 0;
            for (auto l = 0; l < 3; ++l) {
                lhs_.fetch_set(l * m_ + p, l * m_ + p, 1);
            }
        }
    }

    // ГУ первого рода Sigma_2 (y)
    for (auto& n : t_.first()[1]) {
        for (auto i = 0; i < 3; ++i) {
            auto p = n->second;
            rhs_[m_ + p] = 0;
            for (auto l = 0; l < 3; ++l) {
                lhs_.fetch_set(l * m_ + p, l * m_ + p, 1);
            }
        }
    }

    // Сборка матрицы глобальной СЛАУ
    for (auto& k : t_.elems()) {
        auto r = g_->get_internal(t_.data(k));
        for (auto i = 0; i < 4; ++i) {
            auto p = k[i]->second;

            if (!t_.on_first(k[i])) {
                for (auto j = 0; j < 4; ++j) {
                    auto q = k[j]->second;
                    if (!t_.on_first(k[j])) {
                        for (auto l = 0; l < 3; ++l) {
                            lhs_.fetch_add(l * m_ + p, l * m_ + q,
                                           r.first(l * Triangulation::N + i,
                                                   l * Triangulation::N + j));
                        }
                    } else {
                        for (auto l = 0; l < 3; ++l) {
                            rhs_[l * m_ + p]
                                -= r.first(l * Triangulation::N + i,
                                           l * Triangulation::N + j)
                                * rhs_[l * m_ + q];
                        }
                    }
                }
            }
        }
    }

    lhs_.clean_up();
    return *this;
}

SparceMatrix& GlobalEqBuilder::mat()
{
    return lhs_;
}

Vec& GlobalEqBuilder::vec()
{
    return rhs_;
}
