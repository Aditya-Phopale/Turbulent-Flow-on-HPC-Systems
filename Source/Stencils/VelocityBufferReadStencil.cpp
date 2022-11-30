#include "StdAfx.hpp"

#include "VelocityBufferReadStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::VelocityBufferReadStencil::VelocityBufferReadStencil(
  const Parameters&                  parameters,
  std::vector<RealType> Left,
  std::vector<RealType> Right,
  std::vector<RealType> Top,
  std::vector<RealType> Bottom,
  std::vector<RealType> Front,
  std::vector<RealType> Back
):
  BoundaryStencil<FlowField>(parameters),
  Left_(Left),
  Right_(Right),
  Top_(Top),
  Bottom_(Bottom),
  Front_(Front),
  Back_(Back) {}

void Stencils::VelocityBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0]     = Left_.at(j);
  flowField.getVelocity().getVector(i + 1, j)[1] = Left_.at(flowField.getCellsY() + j);
}
void Stencils::VelocityBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = Right_.at(j);
  flowField.getVelocity().getVector(i, j)[1] = Right_.at(flowField.getCellsY() + j);
}
void Stencils::VelocityBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j + 1)[0] = Bottom_.at(i);
  flowField.getVelocity().getVector(i, j)[1]     = Bottom_.at(flowField.getCellsX() + i);
}
void Stencils::VelocityBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = Top_.at(i);
  flowField.getVelocity().getVector(i, j)[1] = Top_.at(flowField.getCellsX() + i);
}

void Stencils::VelocityBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0]     = Left_.at(j * (flowField.getCellsZ()) + k);
  flowField.getVelocity().getVector(i + 1, j, k)[1] = Left_.at(j * (flowField.getCellsZ()) + k + (flowField.getCellsY())*(flowField.getCellsZ()));
  flowField.getVelocity().getVector(i + 1, j, k)[2] = Left_.at(j * (flowField.getCellsZ()) + k + 2*(flowField.getCellsY())*(flowField.getCellsZ()));
}
void Stencils::VelocityBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = Right_.at(j * (flowField.getCellsZ()) + k);
  flowField.getVelocity().getVector(i, j, k)[1] = Right_.at(j * (flowField.getCellsZ()) + k + (flowField.getCellsY())*(flowField.getCellsZ()));
  flowField.getVelocity().getVector(i, j, k)[2] = Right_.at(j * (flowField.getCellsZ()) + k + 2*(flowField.getCellsY())*(flowField.getCellsZ()));
}
void Stencils::VelocityBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j + 1, k)[0] = Bottom_.at(i * (flowField.getCellsZ()) + k);
  flowField.getVelocity().getVector(i, j, k)[1]     = Bottom_.at(i * (flowField.getCellsZ()) + k + (flowField.getCellsX())*(flowField.getCellsZ()));
  flowField.getVelocity().getVector(i, j + 1, k)[2] = Bottom_.at(i * (flowField.getCellsZ()) + k + 2*(flowField.getCellsX())*(flowField.getCellsZ()));
}
void Stencils::VelocityBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[1] = Top_.at(i * (flowField.getCellsZ()) + k);
  flowField.getVelocity().getVector(i, j, k)[1] = Top_.at(i * (flowField.getCellsZ()) + k + (flowField.getCellsX())*(flowField.getCellsZ()));
  flowField.getVelocity().getVector(i, j, k)[1] = Top_.at(i * (flowField.getCellsZ()) + k + 2*(flowField.getCellsX())*(flowField.getCellsZ()));
}
void Stencils::VelocityBufferReadStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k + 1)[0] = Front_.at(i * (flowField.getCellsY()) + j);
  flowField.getVelocity().getVector(i, j, k + 1)[1] = Front_.at(i * (flowField.getCellsY()) + j + (flowField.getCellsX())*(flowField.getCellsY()));
  flowField.getVelocity().getVector(i, j, k)[2]     = Front_.at(i * (flowField.getCellsY()) + j + 2*(flowField.getCellsX())*(flowField.getCellsY()));
}
void Stencils::VelocityBufferReadStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = Back_.at(i * (flowField.getCellsY()) + j);
  flowField.getVelocity().getVector(i, j, k)[1] = Back_.at(i * (flowField.getCellsY()) + j + (flowField.getCellsX())*(flowField.getCellsY()));
  flowField.getVelocity().getVector(i, j, k)[2] = Back_.at(i * (flowField.getCellsY()) + j + 2*(flowField.getCellsX())*(flowField.getCellsY()));
}
