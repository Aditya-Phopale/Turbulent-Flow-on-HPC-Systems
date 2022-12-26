#include "StdAfx.hpp"

#include "nuTKEStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"
#include "TurbulentFlowFieldKE.hpp"

Stencils::nuTStencil::nuTStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowFieldKE>(parameters) {}

void Stencils::nuTStencil::apply(TurbulentFlowFieldKE& flowField, int i, int j) {

  flowField.getnuT().getScalar(i, j) = 5;
}

void Stencils::nuTStencil::apply(TurbulentFlowFieldKE& flowField, int i, int j, int k) {

  flowField.getnuT().getScalar(i, j) = 5;
}