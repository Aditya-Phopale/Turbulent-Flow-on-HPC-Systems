#include "StdAfx.hpp"

#include "nuTStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::nuTStencil::nuTStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) {}

void Stencils::nuTStencil::apply(TurbulentFlowField& flowField, int i, int j) {}

void Stencils::nuTStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {}