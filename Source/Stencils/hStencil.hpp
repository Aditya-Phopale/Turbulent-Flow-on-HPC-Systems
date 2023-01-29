#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {

  template <class FlowFieldType>
  class hStencil: public FieldStencil<FlowFieldType> {
  public:
    hStencil(const Parameters& parameters);
    ~hStencil() override = default;

    void apply(FlowFieldType& flowField, int i, int j) override;
    void apply(FlowFieldType& flowField, int i, int j, int k) override;
  };
} // namespace Stencils
