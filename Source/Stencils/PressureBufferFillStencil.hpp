#pragma once

#include "BoundaryStencil.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "Iterators.hpp"

namespace Stencils {

  class PressureBufferFillStencil: public FieldStencil<FlowField> {
  private:
  //
  public:
    PressureBufferFillStencil(const Parameters& parameters);
    ~PressureBufferFillStencil() override = default;

    // void apply(FlowField& flowField, int i, int j) override;
    // void apply(FlowField& flowField, int i, int j, int k) override;
    void fill();
  };
} // namespace Stencils