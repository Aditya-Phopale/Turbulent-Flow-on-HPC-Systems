#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

#include "Stencils/nuTStencil.hpp"

namespace Stencils {
  class timeStepStencil: public FieldStencil<TurbulentFlowField> {
  private:
    RealType Mindt;

  public:
    timeStepStencil(const Parameters& parameters);
    ~timeStepStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;

    RealType getDt();

    void reset();
  };
} // namespace Stencils
