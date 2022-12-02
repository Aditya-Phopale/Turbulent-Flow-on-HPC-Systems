#include "StdAfx.hpp"

#include "PressureBufferReadStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::PressureBufferReadStencil::PressureBufferReadStencil(
  const Parameters&     parameters,
  std::vector<RealType>& Left,
  std::vector<RealType>& Right,
  std::vector<RealType>& Top,
  std::vector<RealType>& Bottom,
  std::vector<RealType>& Front,
  std::vector<RealType>& Back
):
  BoundaryStencil<FlowField>(parameters),
  Left_(Left),
  Right_(Right),
  Top_(Top),
  Bottom_(Bottom),
  Front_(Front),
  Back_(Back) {}

void Stencils::PressureBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i + 1, j) = Left_.at(j);
}
void Stencils::PressureBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = Right_.at(j);
}
void Stencils::PressureBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j + 1) = Bottom_.at(i);
}
void Stencils::PressureBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = Top_.at(i);
}

void Stencils::PressureBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i + 1, j, k) = Left_.at(j * (parameters_.parallel.localSize[2] + 3) + k);
}
void Stencils::PressureBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = Right_.at(j * (parameters_.parallel.localSize[2] + 3) + k);
}
void Stencils::PressureBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  (flowField.getPressure().getScalar(i, j + 1, k)) = Bottom_.at(i * (parameters_.parallel.localSize[2] + 3) + k);
}
void Stencils::PressureBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = Top_.at(i * (parameters_.parallel.localSize[2] + 3) + k);
}
void Stencils::PressureBufferReadStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k + 1) = Front_.at(i * (parameters_.parallel.localSize[1] + 3) + j);
}
void Stencils::PressureBufferReadStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = Back_.at(i * (parameters_.parallel.localSize[1] + 3) + j);
}
