#ifndef MKE2_INCLUDE_LOCAL_EQ_GEN_HPP
#define MKE2_INCLUDE_LOCAL_EQ_GEN_HPP

#include "abstract_local_eq.hpp"

class LocalEqV17 : public AbstractLocalEq {
public:
    Result
    get_internal(Triangulation::FiniteElement::Data elem) const override;
    Result
    get_boundary(Triangulation::SurfaceElement::Data elem) const override;
};

LocalEqGen::result_type v17(const Triangulation& t,
                            const Triangulation::FiniteElement& elem);

#endif // MKE2_INCLUDE_LOCAL_EQ_GEN_HPP