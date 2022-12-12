#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

#include "Stencils/nuTStencil.hpp"

namespace Stencils {
  class timeStepStencil: public FieldStencil<FlowField> {
  private:
    RealType Mindt;

  public:
    timeStepStencil(const Parameters& parameters);
    ~timeStepStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;

    RealType getDt();

    void reset();
  };
} // namespace Stencils
