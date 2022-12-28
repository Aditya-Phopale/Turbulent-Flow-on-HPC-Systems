#include "StdAfx.hpp"

#include "nuTKEStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"
#include "TurbulentFlowFieldKE.hpp"

Stencils::nuTKEStencil::nuTKEStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowFieldKE>(parameters) {}

void Stencils::nuTKEStencil::apply(TurbulentFlowFieldKE& flowField, int i, int j) {
  float c_mu = 0.09;
  float f_mu, k, e;

  f_mu = computef_mu(flowField, i, j);
  k    = flowField.getk().getScalar(i, j);
  e    = flowField.geteps().getScalar(i, j);

  flowField.getnuT().getScalar(i, j) = c_mu * f_mu * k * k / e;
}
