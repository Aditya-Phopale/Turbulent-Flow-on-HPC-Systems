#include "StdAfx.hpp"

#include "PressureBufferFillStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::PressureBufferFillStencil::
  PressureBufferFillStencil(const Parameters& parameters, std::vector<RealType>& left, std::vector<RealType>& right, std::vector<RealType>& top, std::vector<RealType>& bottom, std::vector<RealType>& front, std::vector<RealType>& back):
  BoundaryStencil<FlowField>(parameters), Left_(left), Right_(right), Top_(top), Bottom_(bottom), Front_(front), Back_(back) {
}

void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  Left_.at(j) = flowField.getPressure().getScalar(i + 2, j);
}

void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  Right_.at(j) = flowField.getPressure().getScalar(i - 1, j);
}
void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  Bottom_.at(i) = flowField.getPressure().getScalar(i, j + 2);
}
void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  Top_.at(i) = flowField.getPressure().getScalar(i, j - 1);
}

void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  Left_.at(j * (parameters_.parallel.localSize[2] + 3) + k) = flowField.getPressure().getScalar(i + 2, j, k);
}
void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  Right_.at(j * (parameters_.parallel.localSize[2] + 3) + k) = flowField.getPressure().getScalar(i - 1, j, k);
}
void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  Bottom_.at(i * (parameters_.parallel.localSize[2] + 3) + k) = (flowField.getPressure().getScalar(i, j + 2, k));
}
void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  Top_.at(i * (parameters_.parallel.localSize[2] + 3) + k) = flowField.getPressure().getScalar(i, j - 1, k);
}
void Stencils::PressureBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  Front_.at(i * (parameters_.parallel.localSize[1] + 3) + j) = flowField.getPressure().getScalar(i, j, k + 2);
}
void Stencils::PressureBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  Back_.at(i * (parameters_.parallel.localSize[1] + 3) + j) = flowField.getPressure().getScalar(i, j, k - 1);
}

