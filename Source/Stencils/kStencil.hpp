#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {
  class kStencil: public FieldStencil<TurbulentFlowField> /*change the flowfield here???*/ {
  private:
    RealType localVelocity_[27 * 3];
    RealType localViscosity_[27 * 3];
    RealType localk_[27 * 3];
    RealType localMeshsize_[27 * 3];

  public:
    kStencil(const Parameters& parameters);
    ~kStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };
} // namespace Stencils
