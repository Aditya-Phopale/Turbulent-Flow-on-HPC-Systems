#include "StdAfx.hpp"

#include "kStencil.hpp"

#include "Definitions.hpp"
#include "FGHStencil.hpp"
#include "StencilFunctions.hpp"

Stencils::kStencil::kStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) /*change flowfield here??*/ {}

void Stencils::TurbulentFGHStencil::apply(TurbulentFlowField& flowField, int i, int j) /*change flowfield here??*/ {
  // Load local velocities into the center layer of the local array
  loadLocalVelocity2D(flowField, localVelocity_, i, j);
  loadLocalViscosity2D(flowField, localViscosity_, i, j);
  loadLocalk2D(flowField, localk_, i, j);
  loadLocalMeshsize2D(parameters_, localMeshsize_, i, j);

  computek2D(localVelocity_, localViscosity_, localk_, localMeshsize_, parameters_, parameters_.timestep.dt);
}

void Stencils::kFGHStencil::apply(
  TurbulentFlowField& flowField, int i, int j, int k
) /*change flowfield here??*/ /*change flowfield here??*/ {}
