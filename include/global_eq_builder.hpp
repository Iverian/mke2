#ifndef MKE2_INCLUDE_GLOBAL_EQ_BUILDER_HPP_
#define MKE2_INCLUDE_GLOBAL_EQ_BUILDER_HPP_

#include "abstract_local_eq.hpp"
#include "csr_matrix.hpp"
#include "triangulation.hpp"

std::pair<CsrMatrix, Vec> build_global_system(const Triangulation& t,
                                              LocalEqGen gen);

std::pair<CsrMatrix, Vec> build_global_system(const Triangulation& t,
                                              AbstractLocalEq& gen);

#endif // MKE2_INCLUDE_GLOBAL_EQ_BUILDER_HPP_