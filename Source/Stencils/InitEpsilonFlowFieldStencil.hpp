#pragma once

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "TurbulentFlowFieldKE.hpp"

namespace Stencils {

  /** Initialises the backward facing step scenario, i.e. sets the flag field.
   */
  class InitEpsilonFlowFieldStencil: public FieldStencil<TurbulentFlowFieldKE> {

  public:
    InitEpsilonFlowFieldStencil(const Parameters& parameters):
      FieldStencil<TurbulentFlowFieldKE>(parameters) {}
    ~InitEpsilonFlowFieldStencil() override = default;

    void apply(TurbulentFlowFieldKE& flowField, int i, int j) override {
      // flowField.geteps().getScalar(i, j) = 0.09 * 0.00015;
      flowField.geteps().getScalar(i, j) = 1e-6;
      //  pow(parameters_.turbulent.cmu, 0.75) * pow(flowField.getk().getScalar(i, j), 1.5)
      //  parameters_.geometry.lengthY;
    }
    void apply(TurbulentFlowFieldKE& flowField, int i, int j, int k) override {
      flowField.geteps().getScalar(
        i, j, k
      ) = pow(parameters_.turbulent.cmu, 0.75) * pow(flowField.getk().getScalar(i, j, k), 1.5)
          / parameters_.geometry.lengthY;
    }
  };

} // namespace Stencils