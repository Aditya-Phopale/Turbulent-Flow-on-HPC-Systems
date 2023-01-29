#pragma once

#include "FieldStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowFieldKE.hpp"

namespace Stencils {
  class nuTStencilKE: public FieldStencil<TurbulentFlowFieldKE> {

  protected:
    const Parameters& parameters_;

  public:
    nuTStencilKE(const Parameters& parameters);
    ~nuTStencilKE() override = default;

    void apply(TurbulentFlowFieldKE& flowField, int i, int j) override;
    void apply(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
  };
} // namespace Stencils
