#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {
  class nuTStencil: public FieldStencil<FlowField> {
  private:
    RealType localVelocity_[27 * 3];
    RealType localMeshsize_[27 * 3];

  public:
    nuTStencil(const Parameters& parameters);
    ~nuTStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;
  };
} // namespace Stencils
