#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {
  class epsilonStencil: public FieldStencil<TurbulentFlowField> /*change the flowfield here???*/ {
  private:
    RealType localVelocity_[27 * 3];
    RealType localViscosity_[27 * 3];
    RealType localEpsilon_[27 * 3];
    RealType localMeshsize_[27 * 3];

  public:
    epsilonStencil(const Parameters& parameters);
    ~epsilonStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };
} // namespace Stencils
