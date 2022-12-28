#include "StdAfx.hpp"

#include "nuTKEStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"
#include "TurbulentFlowFieldKE.hpp"

Stencils::nuTKEStencil::nuTKEStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowFieldKE>(parameters) {}

void Stencils::nuTKEStencil::apply(TurbulentFlowFieldKE& flowField, int i, int j) {

  flowField.getnuT().getScalar(
    i, j
  ) = parameters_.turbulent.cmu * fu(flowField, i, j) * flowField.getk().getScalar(i, j)
      * flowField.getk().getScalar(i, j) / flowField.geteps().getScalar(i, j);
}
