#include "StdAfx.hpp"

#include "kStencil.hpp"

#include "Definitions.hpp"
#include "FGHStencil.hpp"
#include "StencilFunctions.hpp"
#include "TurbulentFlowFieldKE.hpp"

Stencils::kStencil::kStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowFieldKE>(parameters) /*change flowfield here??*/ {}

void Stencils::kStencil::apply(TurbulentFlowFieldKE& flowField, int i, int j) /*change flowfield here??*/ {
  // Load local velocities into the center layer of the local array
  loadLocalVelocity2D(flowField, localVelocity_, i, j);
  loadLocalViscosity2D(flowField, localViscosity_, i, j);
  loadLocalK2D(flowField, localk_, i, j);
  loadLocalMeshsize2D(parameters_, localMeshsize_, i, j);

  flowField.getk().getScalar(i, j) = computek2D(
    flowField, localVelocity_, localViscosity_, localk_, localMeshsize_, parameters_, parameters_.timestep.dt, i, j
  );
}

void Stencils::kStencil::apply(
  TurbulentFlowFieldKE& flowField, int i, int j, int k
) /*change flowfield here??*/ /*change flowfield here??*/ {}
