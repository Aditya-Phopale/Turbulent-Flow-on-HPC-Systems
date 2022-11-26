#include "StdAfx.hpp"

#include "VelocityBufferFillStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::VelocityBufferFillStencil::VelocityBufferFillStencil(
  const Parameters&                   parameters,
  std::vector<RealType>& left,
  std::vector<RealType>& right,
  std::vector<RealType>& top,
  std::vector<RealType>& bottom,
  std::vector<RealType>& front,
  std::vector<RealType>& back
):
  BoundaryStencil<FlowField>(parameters),
  Left_(left),
  Right_(right),
  Top_(top),
  Bottom_(bottom),
  Front_(front),
  Back_(back) {
}

void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  Left_.at(j) = flowField.getVelocity().getVector(i + 2, j)[0];
  Left_.at(flowField.getCellsY() + j) = flowField.getVelocity().getVector(i + 2, j)[1];
}

void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  Right_.at(j) = flowField.getVelocity().getVector(i - 2, j)[0];
  Right_.at(flowField.getCellsY() + j) = flowField.getVelocity().getVector(i - 1, j)[1];
}

void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  Bottom_.at(i) = flowField.getVelocity().getVector(i, j + 2)[0];
  Bottom_.at(flowField.getCellsX() + i) = flowField.getVelocity().getVector(i, j + 2)[1];
}

void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  Top_.at(i) = flowField.getVelocity().getVector(i, j - 1)[0];
  Top_.at(flowField.getCellsX() + i) = flowField.getVelocity().getVector(i, j - 2)[1];
}

void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  Left_.at(j * (parameters_.geometry.sizeZ + 3) + k) = flowField.getVelocity().getVector(i + 2, j, k)[0];
  Left_.at(j * (parameters_.geometry.sizeZ + 3) + k + (parameters_.geometry.sizeY + 3)*(parameters_.geometry.sizeZ + 3)) = flowField.getVelocity().getVector(i + 2, j, k)[1];
  Left_.at(j * (parameters_.geometry.sizeZ + 3) + k + 2*(parameters_.geometry.sizeY + 3)*(parameters_.geometry.sizeZ + 3)) = flowField.getVelocity().getVector(i + 2, j, k)[2];
}
void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  Right_.at(j * (parameters_.geometry.sizeZ + 3) + k) = flowField.getVelocity().getVector(i - 2, j, k)[0];
  Right_.at(j * (parameters_.geometry.sizeZ + 3) + k + (parameters_.geometry.sizeY + 3)*(parameters_.geometry.sizeZ + 3)) = flowField.getVelocity().getVector(i - 1, j, k)[1];
  Right_.at(j * (parameters_.geometry.sizeZ + 3) + k + 2*(parameters_.geometry.sizeY + 3)*(parameters_.geometry.sizeZ + 3)) = flowField.getVelocity().getVector(i - 1, j, k)[2];
}
void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  Bottom_.at(i * (parameters_.geometry.sizeZ + 3) + k) = flowField.getVelocity().getVector(i, j + 2, k)[0];
  Bottom_.at(i * (parameters_.geometry.sizeZ + 3) + k + (parameters_.geometry.sizeX + 3)*(parameters_.geometry.sizeZ + 3)) = flowField.getVelocity().getVector(i, j + 2, k)[1];
  Bottom_.at(i * (parameters_.geometry.sizeZ + 3) + k + 2*(parameters_.geometry.sizeX + 3)*(parameters_.geometry.sizeZ + 3)) = flowField.getVelocity().getVector(i, j + 2, k)[2];
}
void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  Top_.at(i * (parameters_.geometry.sizeZ + 3) + k) = flowField.getVelocity().getVector(i, j - 1, k)[0];
  Top_.at(i * (parameters_.geometry.sizeZ + 3) + k + (parameters_.geometry.sizeX + 3)*(parameters_.geometry.sizeZ + 3)) = flowField.getVelocity().getVector(i, j - 2, k)[1];
  Top_.at(i * (parameters_.geometry.sizeZ + 3) + k + 2*(parameters_.geometry.sizeX + 3)*(parameters_.geometry.sizeZ + 3)) = flowField.getVelocity().getVector(i, j - 1, k)[2];
}
void Stencils::VelocityBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  Front_.at(i * (parameters_.geometry.sizeY + 3) + j) = flowField.getVelocity().getVector(i, j, k + 2)[0];
  Front_.at(i * (parameters_.geometry.sizeY + 3) + j + (parameters_.geometry.sizeX + 3)*(parameters_.geometry.sizeY + 3)) = flowField.getVelocity().getVector(i, j, k + 2)[1];
  Front_.at(i * (parameters_.geometry.sizeY + 3) + j + 2*(parameters_.geometry.sizeX + 3)*(parameters_.geometry.sizeY + 3)) = flowField.getVelocity().getVector(i, j, k + 2)[2];
}
void Stencils::VelocityBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  Back_.at(i * (parameters_.geometry.sizeY + 3) + j) = flowField.getVelocity().getVector(i, j, k - 1)[0];
  Back_.at(i * (parameters_.geometry.sizeY + 3) + j + (parameters_.geometry.sizeX + 3)*(parameters_.geometry.sizeY + 3)) = flowField.getVelocity().getVector(i, j, k - 1)[1];
  Back_.at(i * (parameters_.geometry.sizeY + 3) + j + 2*(parameters_.geometry.sizeX + 3)*(parameters_.geometry.sizeY + 3)) = flowField.getVelocity().getVector(i, j, k - 2)[2];
}
