#include "StdAfx.hpp"

#include "VelocityBufferReadStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::VelocityBufferReadStencil::VelocityBufferReadStencil(
  const Parameters&                  parameters,
  std::vector<std::vector<RealType>> Left,
  std::vector<std::vector<RealType>> Right,
  std::vector<std::vector<RealType>> Top,
  std::vector<std::vector<RealType>> Bottom
):
  BoundaryStencil<FlowField>(parameters),
  Left_(Left),
  Right_(Right),
  Top_(Top),
  Bottom_(Bottom) {
}

Stencils::VelocityBufferReadStencil::VelocityBufferReadStencil(
  const Parameters&                  parameters,
  std::vector<std::vector<RealType>> Left,
  std::vector<std::vector<RealType>> Right,
  std::vector<std::vector<RealType>> Top,
  std::vector<std::vector<RealType>> Bottom,
  std::vector<std::vector<RealType>> Front,
  std::vector<std::vector<RealType>> Back
):
  BoundaryStencil<FlowField>(parameters),
  Left_(Left),
  Right_(Right),
  Top_(Top),
  Bottom_(Bottom),
  Front_(Front),
  Back_(Back) {}

void Stencils::VelocityBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0]     = Left_.at(j).at(0);
  flowField.getVelocity().getVector(i + 1, j)[1] = Left_.at(j).at(1);
}
void Stencils::VelocityBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = Right_.at(j).at(0);
  flowField.getVelocity().getVector(i, j)[1] = Right_.at(j).at(1);
}
void Stencils::VelocityBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j + 1)[0] = Bottom_.at(i).at(0);
  flowField.getVelocity().getVector(i, j)[1]     = Bottom_.at(i).at(1);
}
void Stencils::VelocityBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = Top_.at(i).at(0);
  flowField.getVelocity().getVector(i, j)[1] = Top_.at(i).at(1);
}

void Stencils::VelocityBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0]     = Left_.at(j * (parameters_.geometry.sizeZ + 3) + k).at(0);
  flowField.getVelocity().getVector(i + 1, j, k)[1] = Left_.at(j * (parameters_.geometry.sizeZ + 3) + k).at(1);
  flowField.getVelocity().getVector(i + 1, j, k)[2] = Left_.at(j * (parameters_.geometry.sizeZ + 3) + k).at(2);
}
void Stencils::VelocityBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = Right_.at(j * (parameters_.geometry.sizeZ + 3) + k).at(0);
  flowField.getVelocity().getVector(i, j, k)[1] = Right_.at(j * (parameters_.geometry.sizeZ + 3) + k).at(1);
  flowField.getVelocity().getVector(i, j, k)[2] = Right_.at(j * (parameters_.geometry.sizeZ + 3) + k).at(2);
}
void Stencils::VelocityBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j + 1, k)[0] = Bottom_.at(i * (parameters_.geometry.sizeZ + 3) + k).at(0);
  flowField.getVelocity().getVector(i, j, k)[1]     = Bottom_.at(i * (parameters_.geometry.sizeZ + 3) + k).at(1);
  flowField.getVelocity().getVector(i, j + 1, k)[2] = Bottom_.at(i * (parameters_.geometry.sizeZ + 3) + k).at(2);
}
void Stencils::VelocityBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[1] = Top_.at(i * (parameters_.geometry.sizeZ + 3) + k).at(0);
  flowField.getVelocity().getVector(i, j, k)[1] = Top_.at(i * (parameters_.geometry.sizeZ + 3) + k).at(1);
  flowField.getVelocity().getVector(i, j, k)[1] = Top_.at(i * (parameters_.geometry.sizeZ + 3) + k).at(2);
}
void Stencils::VelocityBufferReadStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k + 1)[0] = Front_.at(i * (parameters_.geometry.sizeY + 3) + j).at(0);
  flowField.getVelocity().getVector(i, j, k + 1)[1] = Front_.at(i * (parameters_.geometry.sizeY + 3) + j).at(1);
  flowField.getVelocity().getVector(i, j, k)[2]     = Front_.at(i * (parameters_.geometry.sizeY + 3) + j).at(2);
}
void Stencils::VelocityBufferReadStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = Back_.at(i * (parameters_.geometry.sizeY + 3) + j).at(0);
  flowField.getVelocity().getVector(i, j, k)[1] = Back_.at(i * (parameters_.geometry.sizeY + 3) + j).at(1);
  flowField.getVelocity().getVector(i, j, k)[2] = Back_.at(i * (parameters_.geometry.sizeY + 3) + j).at(2);
}
