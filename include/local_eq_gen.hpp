#ifndef MKE2_INCLUDE_LOCAL_EQ_GEN_HPP_
#define MKE2_INCLUDE_LOCAL_EQ_GEN_HPP_

#include <dense_matrix.hpp>
#include <triangulation.hpp>
#include <vec.hpp>

#include <functional>

using LocalEqGen = std::function<std::pair<DenseMatrix, Vec>(
    const Triangulation&, const Triangulation::FiniteElement&)>;

LocalEqGen::result_type gen_local(const Triangulation& t,
                                  const Triangulation::FiniteElement& elem);

#endif // MKE2_INCLUDE_LOCAL_EQ_GEN_HPP_