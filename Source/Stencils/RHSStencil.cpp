#include "StdAfx.hpp"
#include "RHSStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::RHSStencil::RHSStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {}

void Stencils::RHSStencil::apply(FlowField& flowField, int i, int j) {

  loadLocalMeshsize2D(parameters_, localMeshsize_, i, j);
  RealType& rhs = flowField.getRHS().getScalar(i, j);

  RealType* const value_p = flowField.getFGH().getVector(i, j);
  RealType* const value_f = flowField.getFGH().getVector(i - 1, j);
  RealType* const value_g = flowField.getFGH().getVector(i, j - 1);

  rhs = computeRHS2D(value_p, value_f, value_g, localMeshsize_, parameters_.timestep.dt);
}

void Stencils::RHSStencil::apply(FlowField& flowField, int i, int j, int k) {
  loadLocalMeshsize3D(parameters_, localMeshsize_, i, j, k);
  RealType& rhs = flowField.getRHS().getScalar(i, j, k);

  RealType* const value_p = flowField.getFGH().getVector(i, j, k);
  RealType* const value_f = flowField.getFGH().getVector(i - 1, j, k);
  RealType* const value_g = flowField.getFGH().getVector(i, j - 1, k);
  RealType* const value_h = flowField.getFGH().getVector(i, j, k - 1);

  rhs = computeRHS3D(value_p, value_f, value_g, value_h, localMeshsize_, parameters_.timestep.dt);
}