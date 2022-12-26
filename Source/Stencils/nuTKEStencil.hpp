#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowFieldKE.hpp"

namespace Stencils {
  class nuTStencil: public FieldStencil<TurbulentFlowFieldKE> {
  private:
    RealType localVelocity_[27 * 3];
    RealType localMeshsize_[27 * 3];

  public:
    nuTStencil(const Parameters& parameters);
    ~nuTStencil() override = default;

    void apply(TurbulentFlowFieldKE& flowField, int i, int j) override;
    void apply(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
  };
} // namespace Stencils
