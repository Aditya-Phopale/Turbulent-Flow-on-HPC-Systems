#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowFieldKE.hpp"

#include "Stencils/nuTStencilKE.hpp"

namespace Stencils {
  class timeStepStencilKE: public FieldStencil<TurbulentFlowFieldKE> {
  private:
    RealType Mindt;

  public:
    timeStepStencilKE(const Parameters& parameters);
    ~timeStepStencilKE() override = default;

    void apply(TurbulentFlowFieldKE& flowField, int i, int j) override;
    void apply(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;

    RealType getDt();

    void reset();
  };
} // namespace Stencils
