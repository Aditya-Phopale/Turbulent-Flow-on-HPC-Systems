#include "StdAfx.hpp"

#include "TurbulentFlowField.hpp"
#include "ViscosityBufferReadStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::ViscosityBufferReadStencil::ViscosityBufferReadStencil(
  const Parameters&     parameters,
  std::vector<RealType>& Left,
  std::vector<RealType>& Right,
  std::vector<RealType>& Top,
  std::vector<RealType>& Bottom,
  std::vector<RealType>& Front,
  std::vector<RealType>& Back
):
  BoundaryStencil<TurbulentFlowField>(parameters),
  Left_(Left),
  Right_(Right),
  Top_(Top),
  Bottom_(Bottom),
  Front_(Front),
  Back_(Back) {}

void Stencils::ViscosityBufferReadStencil::applyLeftWall(TurbulentFlowField& turbFlowField, int i, int j) {
  turbFlowField.getnuT().getScalar(i + 1, j) = Left_.at(j);
}
void Stencils::ViscosityBufferReadStencil::applyRightWall(TurbulentFlowField& turbFlowField, int i, int j) {
  turbFlowField.getnuT().getScalar(i, j) = Right_.at(j);
}
void Stencils::ViscosityBufferReadStencil::applyBottomWall(TurbulentFlowField& turbFlowField, int i, int j) {
  turbFlowField.getnuT().getScalar(i, j + 1) = Bottom_.at(i);
}
void Stencils::ViscosityBufferReadStencil::applyTopWall(TurbulentFlowField& turbFlowField, int i, int j) {
  turbFlowField.getnuT().getScalar(i, j) = Top_.at(i);
}

void Stencils::ViscosityBufferReadStencil::applyLeftWall(TurbulentFlowField& turbFlowField, int i, int j, int k) {
  turbFlowField.getnuT().getScalar(i + 1, j, k) = Left_.at(j * (parameters_.parallel.localSize[2] + 3) + k);
}
void Stencils::ViscosityBufferReadStencil::applyRightWall(TurbulentFlowField& turbFlowField, int i, int j, int k) {
  turbFlowField.getnuT().getScalar(i, j, k) = Right_.at(j * (parameters_.parallel.localSize[2] + 3) + k);
}
void Stencils::ViscosityBufferReadStencil::applyBottomWall(TurbulentFlowField& turbFlowField, int i, int j, int k) {
  (turbFlowField.getnuT().getScalar(i, j + 1, k)) = Bottom_.at(i * (parameters_.parallel.localSize[2] + 3) + k);
}
void Stencils::ViscosityBufferReadStencil::applyTopWall(TurbulentFlowField& turbFlowField, int i, int j, int k) {
  turbFlowField.getnuT().getScalar(i, j, k) = Top_.at(i * (parameters_.parallel.localSize[2] + 3) + k);
}
void Stencils::ViscosityBufferReadStencil::applyFrontWall(TurbulentFlowField& turbFlowField, int i, int j, int k) {
  turbFlowField.getnuT().getScalar(i, j, k + 1) = Front_.at(i * (parameters_.parallel.localSize[1] + 3) + j);
}
void Stencils::ViscosityBufferReadStencil::applyBackWall(TurbulentFlowField& turbFlowField, int i, int j, int k) {
  turbFlowField.getnuT().getScalar(i, j, k) = Back_.at(i * (parameters_.parallel.localSize[1] + 3) + j);
}
