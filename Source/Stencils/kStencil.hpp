#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowFieldKE.hpp"

namespace Stencils {
  class kStencil: public FieldStencil<TurbulentFlowFieldKE> /*change the flowfield here???*/ {
  private:
    RealType localVelocity_[27 * 3];
    RealType localViscosity_[27 * 3];
    RealType localk_[27 * 3];
    RealType localMeshsize_[27 * 3];

  public:
    kStencil(const Parameters& parameters);
    ~kStencil() override = default;

    void apply(TurbulentFlowFieldKE& flowField, int i, int j) override;
    void apply(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
  };
} // namespace Stencils
