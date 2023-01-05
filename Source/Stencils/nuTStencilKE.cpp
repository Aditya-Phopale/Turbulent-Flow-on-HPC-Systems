#include "StdAfx.hpp"

#include "nuTStencilKE.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"
#include "TurbulentFlowFieldKE.hpp"

Stencils::nuTStencilKE::nuTStencilKE(const Parameters& parameters):
  FieldStencil<TurbulentFlowFieldKE>(parameters),
  parameters_(parameters) {}

void Stencils::nuTStencilKE::apply(TurbulentFlowFieldKE& flowField, int i, int j) {
  //

  flowField.getnuT().getScalar(
    i, j
  ) = parameters_.turbulent.cmu * fu(parameters_, flowField, i, j) * flowField.getk().getScalar(i, j)
      * flowField.getk().getScalar(i, j) / (flowField.geteps().getScalar(i, j));
}

void Stencils::nuTStencilKE::apply(
  TurbulentFlowFieldKE& flowField, int i, int j, int k
) /*yet to be defined properly*/ {

  flowField.getnuT().getScalar(
    i, j
  ) = parameters_.turbulent.cmu * fu(parameters_, flowField, i, j) * flowField.getk().getScalar(i, j)
      * flowField.getk().getScalar(i, j) / flowField.geteps().getScalar(i, j);
}
