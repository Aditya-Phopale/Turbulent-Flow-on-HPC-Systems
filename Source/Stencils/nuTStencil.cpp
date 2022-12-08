#include "StdAfx.hpp"

#include "nuTStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::nuTStencil::nuTStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) {}

void Stencils::nuTStencil::apply(TurbulentFlowField& flowField, int i, int j) {
  RealType delta = 0.0;
  if (parameters_.turbulent.delta == 0) {
    delta = 0.0;
  } else if (parameters_.turbulent.delta == 1) {
    auto x = parameters_.meshsize->getPosX(i, j) + 0.5 * parameters_.meshsize->getDx(i, j);
    delta  = 4.91 * sqrt(x / ((parameters_.walls.scalarLeft) * parameters_.flow.Re));
  } else if (parameters_.turbulent.delta == 2) {
    auto x = parameters_.meshsize->getPosX(i, j) + 0.5 * parameters_.meshsize->getDx(i, j);
    delta  = 0.382 * pow(x, 0.8) / pow(parameters_.walls.vectorLeft[0] * parameters_.flow.Re, 0.2);
  }

  RealType lm = std::min(parameters_.turbulent.kappa * flowField.getheight().getScalar(i, j), 0.09 * delta);

  loadLocalVelocity2D(flowField, localVelocity_, i, j);
  loadLocalMeshsize2D(parameters_, localMeshsize_, i, j);

  RealType Sij = pow(dudx(localVelocity_, localMeshsize_), 2) + pow(dvdy(localVelocity_, localMeshsize_), 2)
                 + pow(dudy(localVelocity_, localMeshsize_) + dvdx(localVelocity_, localMeshsize_), 2);

  flowField.getnuT().getScalar(i, j) = lm * lm * sqrt(2 * Sij);
}

void Stencils::nuTStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
  RealType delta = 0.0;
  if (parameters_.turbulent.delta == 0) {
    delta = 0.0;
  } else if (parameters_.turbulent.delta == 1) {
    auto x = parameters_.meshsize->getPosX(i, j) + 0.5 * parameters_.meshsize->getDx(i, j);
    delta  = 4.91 * sqrt(x / ((parameters_.walls.scalarLeft) * parameters_.flow.Re));
  } else {
    auto x = parameters_.meshsize->getPosX(i, j) + 0.5 * parameters_.meshsize->getDx(i, j);
    delta  = 0.382 * pow(x, 0.8) / pow(parameters_.walls.vectorLeft[0] * parameters_.flow.Re, 0.2);
  }
  RealType lm = std::min(parameters_.turbulent.kappa * flowField.getheight().getScalar(i, j, k), 0.09 * delta);

  loadLocalVelocity3D(flowField, localVelocity_, i, j, k);
  loadLocalMeshsize3D(parameters_, localMeshsize_, i, j, k);

  RealType Sij = 0.0;

  // pow(dudx(localVelocity_, localMeshsize_), 2) + pow(dvdy(localVelocity_, localMeshsize_), 2)
  //                + pow(dwdz(localVelocity_, localMeshsize_), 2)
  //                + pow(dudy(localVelocity_, localMeshsize_) + dvdx(localVelocity_, localMeshsize_), 2)
  //                + pow(dudz(localVelocity_, localMeshsize_) + dwdx(localVelocity_, localMeshsize_), 2)
  //                + pow(dvdz(localVelocity_, localMeshsize_) + dwdy(localVelocity_, localMeshsize_),
  //                2);

  flowField.getnuT().getScalar(i, j, k) = lm * lm * sqrt(2 * Sij);
}