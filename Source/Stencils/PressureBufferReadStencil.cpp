#include "StdAfx.hpp"

#include "PressureBufferReadStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::PressureBufferReadStencil::PressureBufferReadStencil(
  const Parameters& parameters, RealType* bL, RealType* bR, RealType* bBo, RealType* bT, RealType* bF, RealType* bB
):
  BoundaryStencil<FlowField>(parameters),
  bLeft(bL),
  bRight(bR),
  bBottom(bB),
  bTop(bT),
  bFront(bF),
  bBack(bB) {}

void Stencils::PressureBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i + 1, j) = bLeft[j];
}
void Stencils::PressureBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = bRight[j];
}
void Stencils::PressureBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j + 1) = bBottom[i];
}
void Stencils::PressureBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = bTop[i];
}

void Stencils::PressureBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i + 2, j, k) = bLeft[j * parameters_.geometry.sizeZ + k];
}
void Stencils::PressureBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  bRight[j * parameters_.geometry.sizeZ + k] = flowField.getPressure().getScalar(i - 1, j, k);
}
void Stencils::PressureBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  bBottom[i * parameters_.geometry.sizeZ + k] = (flowField.getPressure().getScalar(i, j + 2, k));
}
void Stencils::PressureBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  bTop[i * parameters_.geometry.sizeZ + k] = flowField.getPressure().getScalar(i, j - 1, k);
}
void Stencils::PressureBufferReadStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  bFront[i * parameters_.geometry.sizeY + j] = flowField.getPressure().getScalar(i, j, k + 2);
}
void Stencils::PressureBufferReadStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  bBack[i * parameters_.geometry.sizeY + j] = flowField.getPressure().getScalar(i, j, k - 1);
}
