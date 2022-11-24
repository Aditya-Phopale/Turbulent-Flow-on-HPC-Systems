#include "StdAfx.hpp"

#include "PressureBufferFillStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::PressureBufferFillStencil::PressureBufferFillStencil(
  const Parameters& parameters, RealType* bL, RealType* bR, RealType* bBo, RealType* bT, RealType* bF, RealType* bB
):
  BoundaryStencil<FlowField>(parameters),
  bLeft(bL),
  bRight(bR),
  bBottom(bB),
  bTop(bT),
  bFront(bF),
  bBack(bB) {}

void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  bLeft[j] = flowField.getPressure().getScalar(i + 2, j);
}

void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  bRight[j] = flowField.getPressure().getScalar(i - 1, j);
}
void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  bBottom[i] = flowField.getPressure().getScalar(i, j + 2);
}
void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  bTop[i] = flowField.getPressure().getScalar(i, j - 1);
}

void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  bLeft[j * parameters_.geometry.sizeZ + k] = flowField.getPressure().getScalar(i + 2, j, k);
}
void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  bRight[j * parameters_.geometry.sizeZ + k] = flowField.getPressure().getScalar(i - 1, j, k);
}
void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  bBottom[i * parameters_.geometry.sizeZ + k] = (flowField.getPressure().getScalar(i, j + 2, k));
}
void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  bTop[i * parameters_.geometry.sizeZ + k] = flowField.getPressure().getScalar(i, j - 1, k);
}
void Stencils::PressureBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  bFront[i * parameters_.geometry.sizeY + j] = flowField.getPressure().getScalar(i, j, k + 2);
}
void Stencils::PressureBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  bBack[i * parameters_.geometry.sizeY + j] = flowField.getPressure().getScalar(i, j, k - 1);
}
