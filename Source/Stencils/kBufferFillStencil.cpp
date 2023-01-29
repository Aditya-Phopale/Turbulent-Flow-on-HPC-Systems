#include "StdAfx.hpp"

#include "kBufferFillStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::kBufferFillStencil::kBufferFillStencil(
  const Parameters&      parameters,
  std::vector<RealType>& left,
  std::vector<RealType>& right,
  std::vector<RealType>& top,
  std::vector<RealType>& bottom,
  std::vector<RealType>& front,
  std::vector<RealType>& back
):
  BoundaryStencil<TurbulentFlowFieldKE>(parameters),
  Left_(left),
  Right_(right),
  Top_(top),
  Bottom_(bottom),
  Front_(front),
  Back_(back) {}

void Stencils::kBufferFillStencil::applyLeftWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j) {
  Left_.at(j) = turbFlowFieldKE.getk().getScalar(i + 2, j);
}

void Stencils::kBufferFillStencil::applyRightWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j) {
  Right_.at(j) = turbFlowFieldKE.getk().getScalar(i - 1, j);
}
void Stencils::kBufferFillStencil::applyBottomWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j) {
  Bottom_.at(i) = turbFlowFieldKE.getk().getScalar(i, j + 2);
}
void Stencils::kBufferFillStencil::applyTopWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j) {
  Top_.at(i) = turbFlowFieldKE.getk().getScalar(i, j - 1);
}

void Stencils::kBufferFillStencil::applyLeftWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) {
  Left_.at(j * (parameters_.parallel.localSize[2] + 3) + k) = turbFlowFieldKE.getk().getScalar(i + 2, j, k);
}
void Stencils::kBufferFillStencil::applyRightWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) {
  Right_.at(j * (parameters_.parallel.localSize[2] + 3) + k) = turbFlowFieldKE.getk().getScalar(i - 1, j, k);
}
void Stencils::kBufferFillStencil::applyBottomWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) {
  Bottom_.at(i * (parameters_.parallel.localSize[2] + 3) + k) = (turbFlowFieldKE.getk().getScalar(i, j + 2, k));
}
void Stencils::kBufferFillStencil::applyTopWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) {
  Top_.at(i * (parameters_.parallel.localSize[2] + 3) + k) = turbFlowFieldKE.getk().getScalar(i, j - 1, k);
}
void Stencils::kBufferFillStencil::applyFrontWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) {
  Front_.at(i * (parameters_.parallel.localSize[1] + 3) + j) = turbFlowFieldKE.getk().getScalar(i, j, k + 2);
}
void Stencils::kBufferFillStencil::applyBackWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) {
  Back_.at(i * (parameters_.parallel.localSize[1] + 3) + j) = turbFlowFieldKE.getk().getScalar(i, j, k - 1);
}
