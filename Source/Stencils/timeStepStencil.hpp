#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"
#include "TurbulentFlowFieldKE.hpp"

#include "Stencils/nuTStencil.hpp"
#include "Stencils/nuTStencilKE.hpp"

namespace Stencils {
  template <class FlowFieldType>
  class timeStepStencil: public FieldStencil<FlowFieldType> {
  private:
    RealType Mindt;

  public:
    timeStepStencil(const Parameters& parameters);
    ~timeStepStencil() override = default;

    void apply(FlowFieldType& flowField, int i, int j) override;
    void apply(FlowFieldType& flowField, int i, int j, int k) override;

    RealType getDt();

    void reset();
  };
} // namespace Stencils
