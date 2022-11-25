#include "StdAfx.hpp"

#include "VelocityBufferReadStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::VelocityBufferReadStencil::VelocityBufferReadStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters) {
  if (parameters_.geometry.dim == 2) {
    LeftRight.resize((parameters_.geometry.sizeY + 3), std::vector<RealType>(2));
    TopBottom.resize((parameters_.geometry.sizeX + 3), std::vector<RealType>(2));
  }

  if (parameters_.geometry.dim == 3) {
    LeftRight.resize((parameters_.geometry.sizeY + 3) * (parameters_.geometry.sizeZ + 3), std::vector<RealType>(3));
    TopBottom.resize((parameters_.geometry.sizeX + 3) * (parameters_.geometry.sizeZ + 3), std::vector<RealType>(3));
    FrontBack.resize((parameters_.geometry.sizeX + 3) * (parameters_.geometry.sizeY + 3), std::vector<RealType>(3));
  }
}

void Stencils::VelocityBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i + 1, j)[1] = LeftRight.at(j);
}
void Stencils::VelocityBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[1] = LeftRight[j];
}
void Stencils::VelocityBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j + 1)[1] = bBottom_v[i];
}
void Stencils::VelocityBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[1] = bTop_v[i];
}

void Stencils::VelocityBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i + 1, j, k)[1] = LeftRight[j * parameters_.geometry.sizeZ + k];
}
void Stencils::VelocityBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[1] = LeftRight_v[j * parameters_.geometry.sizeZ + k];
}
void Stencils::VelocityBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  (flowField.getVelocity().getVector(i, j + 1, k))[1] = bBottom_v[i * parameters_.geometry.sizeZ + k];
}
void Stencils::VelocityBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[1] = bTop_v[i * parameters_.geometry.sizeZ + k];
}
void Stencils::VelocityBufferReadStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k + 1)[1] = bFront_v[i * parameters_.geometry.sizeY + j];
}
void Stencils::VelocityBufferReadStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[1] = bBack_v[i * parameters_.geometry.sizeY + j];
}
