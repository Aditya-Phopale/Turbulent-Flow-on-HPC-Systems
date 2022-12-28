#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowFieldKE.hpp"

namespace Stencils {
  class nuTKEStencil: public FieldStencil<TurbulentFlowFieldKE> {

  public:
    nuTKEStencil(const Parameters& parameters);
    ~nuTKEStencil() override = default;

    void apply(TurbulentFlowFieldKE& flowField, int i, int j) override;
    void apply(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
  };
} // namespace Stencils
