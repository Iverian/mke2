#include <dok_matrix.hpp>

using namespace std;

DokMatrix::DokMatrix(Index side)
    : data_()
    , m_(side)
{
}

DokMatrix::Index DokMatrix::size() const
{
    return m_ * m_;
}

DokMatrix::Shape DokMatrix::shape() const
{
    return {m_, m_};
}

DokMatrix::Index DokMatrix::non_zero() const noexcept
{
    return data_.size();
}

DokMatrix::Value DokMatrix::operator()(Index i, Index j) const
{
    check_if(i < m_ && j < m_, "Index out of range");

    auto pos = data_.find({i, j});

    return pos != end(data_) ? pos->second : Value(0);
}

DokMatrix::Value DokMatrix::add(Index i, Index j, Value val)
{
    return modify({i, j}, [&val](auto x) { return x + val; });
}

DokMatrix::Value DokMatrix::sub(Index i, Index j, Value val)
{
    return modify({i, j}, [&val](auto x) { return x - val; });
}

DokMatrix::Value DokMatrix::set(Index i, Index j, Value val)
{
    return modify({i, j}, [&val](auto) { return val; });
}

const DokMatrix::DataContainer& DokMatrix::data() const noexcept
{
    return data_;
}
