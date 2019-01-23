#ifndef MKE2_INCLUDE_DOK_MATRIX_HPP_
#define MKE2_INCLUDE_DOK_MATRIX_HPP_

#include "abstract_matrix.hpp"
#include "debug.hpp"
#include "util.hpp"

#include <map>

class DokMatrix : public AbstractMatrix {
public:
    using DataContainer = std::map<Index2d, Value>;

    explicit DokMatrix(Index side = 0);
    Index size() const override;
    Index2d shape() const override;
    Index non_zero() const noexcept;

    Value operator()(Index i, Index j) const;
    Value add(Index i, Index j, Value val);
    Value sub(Index i, Index j, Value val);
    Value set(Index i, Index j, Value val);
    const DataContainer& data() const noexcept;

private:
    template <class Callable>
    Value modify(Index2d pos, Callable&& f);

    DataContainer data_;
    Index m_;
};

template <class Callable>
Value DokMatrix::modify(Index2d ind, Callable&& f)
{
    check_if(ind.first < m_ && ind.second < m_, "Index out of range");

    Value result = 0;
    if (auto new_val = f(0); !iszero(new_val)) {
        auto [pos, inserted] = data_.insert({ind, new_val});
        if (!inserted) {
            result = pos->second;
            pos->second = f(pos->second);
        }
    }
    return result;
}

#endif // MKE2_INCLUDE_DOK_MATRIX_HPP_