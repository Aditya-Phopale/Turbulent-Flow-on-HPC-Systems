#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"
#include "TurbulentFlowFieldKE.hpp"

namespace Stencils {
  class epsilonStencil: public FieldStencil<TurbulentFlowFieldKE> /*change the flowfield here???*/ {
  private:
    RealType localVelocity_[27 * 3];
    RealType localViscosity_[27 * 3];
    RealType localEpsilon_[27 * 3];
    RealType localMeshsize_[27 * 3];
    RealType localk_[27 * 3];

  public:
    epsilonStencil(const Parameters& parameters);
    ~epsilonStencil() override = default;

    void apply(TurbulentFlowFieldKE& flowField, int i, int j) override;
    void apply(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
  };
} // namespace Stencils
