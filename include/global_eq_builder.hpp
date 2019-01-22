#ifndef MKE2_INCLUDE_GLOBAL_EQ_BUILDER_HPP_
#define MKE2_INCLUDE_GLOBAL_EQ_BUILDER_HPP_

#include "local_eq_gen.hpp"
#include "csr_matrix.hpp"
#include "triangulation.hpp"

std::pair<DenseMatrix, Vec> build_global_system_dense(const Triangulation& t,
                                                      LocalEqGen gen);

std::pair<CsrMatrix, Vec> build_global_system(const Triangulation& t,
                                              LocalEqGen gen);

#endif // MKE2_INCLUDE_GLOBAL_EQ_BUILDER_HPP_