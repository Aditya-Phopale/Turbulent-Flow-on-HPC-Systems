#include "StdAfx.hpp"

#include "kBufferReadStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"
#include "TurbulentFlowFieldKE.hpp"

Stencils::kBufferReadStencil::kBufferReadStencil(
  const Parameters&      parameters,
  std::vector<RealType>& Left,
  std::vector<RealType>& Right,
  std::vector<RealType>& Top,
  std::vector<RealType>& Bottom,
  std::vector<RealType>& Front,
  std::vector<RealType>& Back
):
  BoundaryStencil<TurbulentFlowFieldKE>(parameters),
  Left_(Left),
  Right_(Right),
  Top_(Top),
  Bottom_(Bottom),
  Front_(Front),
  Back_(Back) {}

void Stencils::kBufferReadStencil::applyLeftWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j) {
  turbFlowFieldKE.getk().getScalar(i + 1, j) = Left_.at(j);
}
void Stencils::kBufferReadStencil::applyRightWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j) {
  turbFlowFieldKE.getk().getScalar(i, j) = Right_.at(j);
}
void Stencils::kBufferReadStencil::applyBottomWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j) {
  turbFlowFieldKE.getk().getScalar(i, j + 1) = Bottom_.at(i);
}
void Stencils::kBufferReadStencil::applyTopWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j) {
  turbFlowFieldKE.getk().getScalar(i, j) = Top_.at(i);
}

void Stencils::kBufferReadStencil::applyLeftWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) {
  turbFlowFieldKE.getk().getScalar(i + 1, j, k) = Left_.at(j * (parameters_.parallel.localSize[2] + 3) + k);
}
void Stencils::kBufferReadStencil::applyRightWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) {
  turbFlowFieldKE.getk().getScalar(i, j, k) = Right_.at(j * (parameters_.parallel.localSize[2] + 3) + k);
}
void Stencils::kBufferReadStencil::applyBottomWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) {
  (turbFlowFieldKE.getk().getScalar(i, j + 1, k)) = Bottom_.at(i * (parameters_.parallel.localSize[2] + 3) + k);
}
void Stencils::kBufferReadStencil::applyTopWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) {
  turbFlowFieldKE.getk().getScalar(i, j, k) = Top_.at(i * (parameters_.parallel.localSize[2] + 3) + k);
}
void Stencils::kBufferReadStencil::applyFrontWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) {
  turbFlowFieldKE.getk().getScalar(i, j, k + 1) = Front_.at(i * (parameters_.parallel.localSize[1] + 3) + j);
}
void Stencils::kBufferReadStencil::applyBackWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) {
  turbFlowFieldKE.getk().getScalar(i, j, k) = Back_.at(i * (parameters_.parallel.localSize[1] + 3) + j);
}
