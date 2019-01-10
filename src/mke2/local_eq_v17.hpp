#ifndef MKE2_SRC_MKE2_VAR17_H_
#define MKE2_SRC_MKE2_VAR17_H_

#include <abstract_local_eq.hpp>

class LocalEqV17 : public AbstractLocalEq {
public:
    LocalEqV17(const Triang3d::Node& node);

    DenseMatrix get_internal_mat() const override;
    DenseMatrix get_boundary_mat() const override;
    Vec get_internal_vec() const override;
    Vec get_boundary_vec() const override;

private:
    DenseMatrix get_dq_mat() const;
    DenseMatrix get_gk_mat() const;
    DenseMatrix get_s0_mat() const;
    DenseMatrix get_bskeleton() const;
    DenseMatrix get_face_bmat(size_t index) const;

    static DenseMatrix cm;

    Triang3d::Node p_;

    double v_;
    DenseMatrix s0_;
};

#endif // MKE2_SRC_MKE2_VAR17_H_