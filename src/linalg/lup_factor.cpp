#include <debug.hpp>
#include <lup_factor.hpp>
#include <util.hpp>

#include <cmath>

using namespace std;

LupFactor::LupFactor(const DenseMatrix& in)
    : m_()
    , mat_(in)
    , pivot_()
{
}

LupFactor::LupFactor(DenseMatrix&& in)
    : m_()
    , mat_(in)
    , pivot_()
{
}

LupFactor& LupFactor::factor()
{
    auto s = mat_.shape();
    check_if(s.m == s.n,
             "Не могу посчитать LUP разложение неквадратной матрицы");

    m_ = s.m;
    pivot_ = PivotType(m_ + 1);
    for (Index i = 0; i < m_ + 1; ++i) {
        pivot_[i] = i;
    }

    for (Index i = 0; i < m_; ++i) {
        double mmax = 0.0;
        auto imax = i;

        for (Index k = i; k < m_; ++k)
            if (auto mabs = fabs(mat_(k, i)); mabs > mmax) {
                mmax = mabs;
                imax = k;
            }

        check_if(!isnear(mmax, 0), "Не могу разложить вырожденную матрицу");

        if (imax != i) {
            swap(pivot_[i], pivot_[imax]);
            mat_.swap_rows(i, imax);
            ++pivot_[m_];
        }

        for (Index j = i + 1; j < m_; ++j) {
            mat_(j, i) /= mat_(i, i);
            for (Index k = i + 1; k < m_; ++k) {
                mat_(j, k) -= mat_(j, i) * mat_(i, k);
            }
        }
    }

    return *this;
}

const DenseMatrix& LupFactor::mat() const
{
    return mat_;
}

const LupFactor::PivotType& LupFactor::pivot() const
{
    return pivot_;
}

Vec LupFactor::solve(const Vec& v) const
{
    Vec result(v.size());

    for (Index i = 0; i < m_; ++i) {
        result[i] = v[pivot_[i]];
        for (Index k = 0; k < i; ++k) {
            result[i] -= mat_(i, k) * result[k];
        }
    }

    for (Index i = m_ - 1; i != Index(-1); --i) {
        for (Index k = i + 1; k < m_; ++k) {
            result[i] -= mat_(i, k) * result[k];
        }
        result[i] /= mat_(i, i);
    }

    return result;
}

DenseMatrix LupFactor::inverse() const
{
    DenseMatrix result(mat_.shape());

    for (Index j = 0; j < m_; ++j) {
        for (Index i = 0; i < m_; ++i) {
            if (pivot_[i] == j) {
                result(i, j) = 1;
            } else {
                result(i, j) = 0;
            }

            for (Index k = 0; k < i; ++k) {
                result(i, j) -= mat_(i, k) * result(k, j);
            }
        }

        for (Index i = m_ - 1; i != Index(-1); --i) {
            for (Index k = i + 1; k < m_; ++k) {
                result(i, j) -= mat_(i, k) * result(k, j);
            }
            result(i, j) /= mat_(i, i);
        }
    }

    return result;
}

double LupFactor::det() const
{
    auto result = mat_(0, 0);

    for (Index i = 1; i < m_; ++i) {
        result *= mat_(i, i);
    }

    if ((pivot_[m_] - m_) % 2 == 1) {
        result *= -1;
    }

    return result;
}