#include "StdAfx.hpp"

#include "epsilonStencil.hpp"

#include "Definitions.hpp"
#include "FGHStencil.hpp"
#include "StencilFunctions.hpp"

Stencils::epsilonStencil::epsilonStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) /*change flowfield here??*/ {}

void Stencils::epsilonStencil::apply(TurbulentFlowField& flowField, int i, int j) /*change flowfield here??*/ {
  // Load local velocities into the center layer of the local array
  loadLocalVelocity2D(flowField, localVelocity_, i, j);
  loadLocalViscosity2D(flowField, localViscosity_, i, j);
  loadLocalEpsilon2D(flowField, localEpsilon_, i, j);
  loadLocalMeshsize2D(parameters_, localMeshsize_, i, j);
  computeEpsilon2D(
    flowField,
    localVelocity_,
    localViscosity_,
    localEpsilon_,
    localMeshsize_,
    parameters_,
    parameters_.timestep.dt,
    int i,
    int j
  );
}

void Stencils::kFGHStencil::apply(
  TurbulentFlowField& flowField, int i, int j, int k
) /*change flowfield here??*/ /*change flowfield here??*/ {}
