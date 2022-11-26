#include "StdAfx.hpp"

#include "VelocityBufferFillStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::VelocityBufferFillStencil::VelocityBufferFillStencil(
  const Parameters&                   parameters,
  std::vector<std::vector<RealType>>& left,
  std::vector<std::vector<RealType>>& right,
  std::vector<std::vector<RealType>>& top,
  std::vector<std::vector<RealType>>& bottom,
  std::vector<std::vector<RealType>>& front,
  std::vector<std::vector<RealType>>& back
):
  BoundaryStencil<FlowField>(parameters),
  Left_(left),
  Right_(right),
  Top_(top),
  Bottom_(bottom),
  Front_(front),
  Back_(back) {
  // if (parameters_.geometry.dim == 2) {
  //   Left_.resize(2,std::vector<RealType>(parameters_.geometry.sizeY + 3));
  //   Right_.resize(2,std::vector<RealType>(parameters_.geometry.sizeY + 3));
  //   Top_.resize(2, std::vector<RealType>(parameters_.geometry.sizeX + 3));
  //   Bottom_.resize(2, std::vector<RealType>(parameters_.geometry.sizeX + 3));
  // }

  // if (parameters_.geometry.dim == 3) {
  //   Left_.resize(3,  std::vector<RealType>((parameters_.geometry.sizeY + 3) * (parameters_.geometry.sizeZ + 3)));
  //   Right_.resize(3,  std::vector<RealType>((parameters_.geometry.sizeY + 3) * (parameters_.geometry.sizeZ + 3)));
  //   Top_.resize(3,  std::vector<RealType>((parameters_.geometry.sizeX + 3) * (parameters_.geometry.sizeZ + 3)));
  //   Bottom_.resize(3,  std::vector<RealType>((parameters_.geometry.sizeX + 3) * (parameters_.geometry.sizeZ + 3)));
  //   Front_.resize(3, std::vector<RealType>((parameters_.geometry.sizeX + 3) * (parameters_.geometry.sizeY + 3)));
  //   Back_.resize(3, std::vector<RealType>((parameters_.geometry.sizeX + 3) * (parameters_.geometry.sizeY + 3)));
  // }
}

void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  Left_.at(j).at(0) = flowField.getVelocity().getVector(i + 2, j)[0];
  Left_.at(j).at(1) = flowField.getVelocity().getVector(i + 2, j)[1];
}

void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  Right_.at(j).at(0) = flowField.getVelocity().getVector(i - 2, j)[0];
  Right_.at(j).at(1) = flowField.getVelocity().getVector(i - 1, j)[1];
}

void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  Bottom_.at(i).at(0) = flowField.getVelocity().getVector(i, j + 2)[0];
  Bottom_.at(i).at(1) = flowField.getVelocity().getVector(i, j + 2)[1];
}

void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  Top_.at(i).at(0) = flowField.getVelocity().getVector(i, j - 1)[0];
  Top_.at(i).at(1) = flowField.getVelocity().getVector(i, j - 2)[1];
}

void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  Left_.at(j * (parameters_.geometry.sizeZ + 3) + k).at(0) = flowField.getVelocity().getVector(i + 2, j, k)[0];
  Left_.at(j * (parameters_.geometry.sizeZ + 3) + k).at(1) = flowField.getVelocity().getVector(i + 2, j, k)[1];
  Left_.at(j * (parameters_.geometry.sizeZ + 3) + k).at(2) = flowField.getVelocity().getVector(i + 2, j, k)[2];
}
void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  Right_.at(j * (parameters_.geometry.sizeZ + 3) + k).at(0) = flowField.getVelocity().getVector(i - 2, j, k)[0];
  Right_.at(j * (parameters_.geometry.sizeZ + 3) + k).at(1) = flowField.getVelocity().getVector(i - 1, j, k)[1];
  Right_.at(j * (parameters_.geometry.sizeZ + 3) + k).at(2) = flowField.getVelocity().getVector(i - 1, j, k)[2];
}
void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  Bottom_.at(i * (parameters_.geometry.sizeZ + 3) + k).at(0) = flowField.getVelocity().getVector(i, j + 2, k)[0];
  Bottom_.at(i * (parameters_.geometry.sizeZ + 3) + k).at(1) = flowField.getVelocity().getVector(i, j + 2, k)[1];
  Bottom_.at(i * (parameters_.geometry.sizeZ + 3) + k).at(2) = flowField.getVelocity().getVector(i, j + 2, k)[2];
}
void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  Top_.at(i * (parameters_.geometry.sizeZ + 3) + k).at(0) = flowField.getVelocity().getVector(i, j - 1, k)[0];
  Top_.at(i * (parameters_.geometry.sizeZ + 3) + k).at(1) = flowField.getVelocity().getVector(i, j - 2, k)[1];
  Top_.at(i * (parameters_.geometry.sizeZ + 3) + k).at(2) = flowField.getVelocity().getVector(i, j - 1, k)[2];
}
void Stencils::VelocityBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  Front_.at(i * (parameters_.geometry.sizeY + 3) + j).at(0) = flowField.getVelocity().getVector(i, j, k + 2)[0];
  Front_.at(i * (parameters_.geometry.sizeY + 3) + j).at(1) = flowField.getVelocity().getVector(i, j, k + 2)[1];
  Front_.at(i * (parameters_.geometry.sizeY + 3) + j).at(2) = flowField.getVelocity().getVector(i, j, k + 2)[2];
}
void Stencils::VelocityBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  Back_.at(i * (parameters_.geometry.sizeY + 3) + j).at(0) = flowField.getVelocity().getVector(i, j, k - 1)[0];
  Back_.at(i * (parameters_.geometry.sizeY + 3) + j).at(1) = flowField.getVelocity().getVector(i, j, k - 1)[1];
  Back_.at(i * (parameters_.geometry.sizeY + 3) + j).at(2) = flowField.getVelocity().getVector(i, j, k - 2)[2];
}
