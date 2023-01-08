#pragma once

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "TurbulentFlowFieldKE.hpp"

namespace Stencils {

  /** Initialises the backward facing step scenario, i.e. sets the flag field.
   */
  class InitkFlowFieldStencil: public FieldStencil<TurbulentFlowFieldKE> {

  public:
    InitkFlowFieldStencil(const Parameters& parameters):
      FieldStencil<TurbulentFlowFieldKE>(parameters) {}
    ~InitkFlowFieldStencil() override = default;

    void apply(TurbulentFlowFieldKE& flowField, int i, int j) override {
      flowField.getk().getScalar(i, j) = 0.003;
      // flowField.getk().getScalar(i, j) = 1e-5;
      //  1.5 * (parameters_.turbulent.I * parameters_.walls.vectorLeft[0]);
    }
    void apply(TurbulentFlowFieldKE& flowField, int i, int j, int k) override {
      flowField.getk().getScalar(i, j, k) = 1.5 * (parameters_.turbulent.I * parameters_.walls.vectorLeft[0]);
    }
  };

} // namespace Stencils