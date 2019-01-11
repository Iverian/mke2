#ifndef MKE2_INCLUDE_LOCAL_EQ_V17_HPP
#define MKE2_INCLUDE_LOCAL_EQ_V17_HPP

#include "abstract_local_eq.hpp"

class LocalEqV17 : public AbstractLocalEq {
public:
    static const DenseMatrix cm;

    Result get_internal(Triangulation::FiniteElementData elem) const override;
    Result get_boundary(Triangulation::SurfaceElementData elem) const override;

protected:
    DenseMatrix get_dq_mat() const;
    DenseMatrix get_gk_mat() const;
    DenseMatrix get_s0_mat() const;
    DenseMatrix get_bskeleton() const;
    DenseMatrix get_face_bmat(size_t index) const;
};

#endif // MKE2_INCLUDE_LOCAL_EQ_V17_HPP