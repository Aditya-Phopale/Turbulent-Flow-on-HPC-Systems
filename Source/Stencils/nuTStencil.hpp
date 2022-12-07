#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {
  class nuTStencil: public FieldStencil<TurbulentFlowField> {
  public:
    nuTStencil(const Parameters& parameters);
    ~nuTStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };
} // namespace Stencils
