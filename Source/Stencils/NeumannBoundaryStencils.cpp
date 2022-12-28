#include "StdAfx.hpp"

#include "NeumannBoundaryStencils.hpp"

Stencils::NeumannVelocityBoundaryStencil::NeumannVelocityBoundaryStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters) {}

void Stencils::NeumannVelocityBoundaryStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i - 1, j)[0] = flowField.getVelocity().getVector(i, j)[0];
  flowField.getVelocity().getVector(i, j)[1]     = flowField.getVelocity().getVector(i + 1, j)[1];
}

void Stencils::NeumannVelocityBoundaryStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = flowField.getVelocity().getVector(i - 1, j)[0];
  flowField.getVelocity().getVector(i, j)[1] = flowField.getVelocity().getVector(i - 1, j)[1];
}

void Stencils::NeumannVelocityBoundaryStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0]     = flowField.getVelocity().getVector(i, j + 1)[0];
  flowField.getVelocity().getVector(i, j - 1)[1] = flowField.getVelocity().getVector(i, j)[1];
}

void Stencils::NeumannVelocityBoundaryStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = flowField.getVelocity().getVector(i, j - 1)[0];
  flowField.getVelocity().getVector(i, j)[1] = flowField.getVelocity().getVector(i, j - 1)[1];
}

void Stencils::NeumannVelocityBoundaryStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i - 1, j, k)[0] = flowField.getVelocity().getVector(i, j, k)[0];
  flowField.getVelocity().getVector(i, j, k)[1]     = flowField.getVelocity().getVector(i + 1, j, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2]     = flowField.getVelocity().getVector(i + 1, j, k)[2];
}

void Stencils::NeumannVelocityBoundaryStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = flowField.getVelocity().getVector(i - 1, j, k)[0];
  flowField.getVelocity().getVector(i, j, k)[1] = flowField.getVelocity().getVector(i - 1, j, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2] = flowField.getVelocity().getVector(i - 1, j, k)[2];
}

void Stencils::NeumannVelocityBoundaryStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0]     = flowField.getVelocity().getVector(i, j + 1, k)[0];
  flowField.getVelocity().getVector(i, j - 1, k)[1] = flowField.getVelocity().getVector(i, j, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2]     = flowField.getVelocity().getVector(i, j + 1, k)[2];
}

void Stencils::NeumannVelocityBoundaryStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = flowField.getVelocity().getVector(i, j - 1, k)[0];
  flowField.getVelocity().getVector(i, j, k)[1] = flowField.getVelocity().getVector(i, j - 1, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2] = flowField.getVelocity().getVector(i, j - 1, k)[2];
}

void Stencils::NeumannVelocityBoundaryStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0]     = flowField.getVelocity().getVector(i, j, k + 1)[0];
  flowField.getVelocity().getVector(i, j, k)[1]     = flowField.getVelocity().getVector(i, j, k + 1)[1];
  flowField.getVelocity().getVector(i, j, k - 1)[2] = flowField.getVelocity().getVector(i, j, k)[2];
}

void Stencils::NeumannVelocityBoundaryStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = flowField.getVelocity().getVector(i, j, k - 1)[0];
  flowField.getVelocity().getVector(i, j, k)[1] = flowField.getVelocity().getVector(i, j, k - 1)[1];
  flowField.getVelocity().getVector(i, j, k)[2] = flowField.getVelocity().getVector(i, j, k - 1)[2];
}

Stencils::NeumannFGHBoundaryStencil::NeumannFGHBoundaryStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters) {}

// These are left empty. The right values should be computed by the FGH body stencil.
void Stencils::NeumannFGHBoundaryStencil::applyLeftWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::NeumannFGHBoundaryStencil::applyRightWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::NeumannFGHBoundaryStencil::applyBottomWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::NeumannFGHBoundaryStencil::applyTopWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}

void Stencils::NeumannFGHBoundaryStencil::applyLeftWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::NeumannFGHBoundaryStencil::applyRightWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::NeumannFGHBoundaryStencil::applyBottomWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::NeumannFGHBoundaryStencil::applyTopWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::NeumannFGHBoundaryStencil::applyFrontWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::NeumannFGHBoundaryStencil::applyBackWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}

Stencils::NeumannKBoundaryStencil::NeumannKBoundaryStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowFieldKE>(parameters) {}

void Stencils::NeumannKBoundaryStencil::applyLeftWall(TurbulentFlowFieldKE& flowField, int i, int j) {
  flowField.getk().getScalar(i, j) = flowField.getk().getScalar(i + 1, j);
  // flowField.getVelocity().getVector(i, j)[1]     = flowField.getVelocity().getVector(i + 1, j)[1];
}

void Stencils::NeumannKBoundaryStencil::applyRightWall(TurbulentFlowFieldKE& flowField, int i, int j) {
  flowField.getk().getScalar(i, j) = flowField.getk().getScalar(i - 1, j);
  // flowField.getVelocity().getVector(i, j)[0] = flowField.getVelocity().getVector(i - 1, j)[0];
  // flowField.getVelocity().getVector(i, j)[1] = flowField.getVelocity().getVector(i - 1, j)[1];
}

void Stencils::NeumannKBoundaryStencil::applyBottomWall(TurbulentFlowFieldKE& flowField, int i, int j) {
  flowField.getk().getScalar(i, j) = flowField.getk().getScalar(i, j + 1);
  // flowField.getVelocity().getVector(i, j)[0]     = flowField.getVelocity().getVector(i, j + 1)[0];
  // flowField.getVelocity().getVector(i, j - 1)[1] = flowField.getVelocity().getVector(i, j)[1];
}

void Stencils::NeumannKBoundaryStencil::applyTopWall(TurbulentFlowFieldKE& flowField, int i, int j) {
  flowField.getk().getScalar(i, j) = flowField.getk().getScalar(i, j - 1);
  // flowField.getVelocity().getVector(i, j)[0] = flowField.getVelocity().getVector(i, j - 1)[0];
  // flowField.getVelocity().getVector(i, j)[1] = flowField.getVelocity().getVector(i, j - 1)[1];
}

void Stencils::NeumannKBoundaryStencil::applyLeftWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
  flowField.getk().getScalar(i, j, k) = flowField.getk().getScalar(i + 1, j, k);
  // flowField.getVelocity().getVector(i - 1, j, k)[0] = flowField.getVelocity().getVector(i, j, k)[0];
  // flowField.getVelocity().getVector(i, j, k)[1]     = flowField.getVelocity().getVector(i + 1, j, k)[1];
  // flowField.getVelocity().getVector(i, j, k)[2]     = flowField.getVelocity().getVector(i + 1, j, k)[2];
}

void Stencils::NeumannKBoundaryStencil::applyRightWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
  flowField.getk().getScalar(i, j, k) = flowField.getk().getScalar(i - 1, j, k);
  // flowField.getVelocity().getVector(i, j, k)[0] = flowField.getVelocity().getVector(i - 1, j, k)[0];
  // flowField.getVelocity().getVector(i, j, k)[1] = flowField.getVelocity().getVector(i - 1, j, k)[1];
  // flowField.getVelocity().getVector(i, j, k)[2] = flowField.getVelocity().getVector(i - 1, j, k)[2];
}

void Stencils::NeumannKBoundaryStencil::applyBottomWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
  flowField.getk().getScalar(i, j, k) = flowField.getk().getScalar(i, j + 1, k);
  // flowField.getVelocity().getVector(i, j, k)[0]     = flowField.getVelocity().getVector(i, j + 1, k)[0];
  // flowField.getVelocity().getVector(i, j - 1, k)[1] = flowField.getVelocity().getVector(i, j, k)[1];
  // flowField.getVelocity().getVector(i, j, k)[2]     = flowField.getVelocity().getVector(i, j + 1, k)[2];
}

void Stencils::NeumannKBoundaryStencil::applyTopWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
  flowField.getk().getScalar(i, j, k) = flowField.getk().getScalar(i, j - 1, k);
  // flowField.getVelocity().getVector(i, j, k)[0] = flowField.getVelocity().getVector(i, j - 1, k)[0];
  // flowField.getVelocity().getVector(i, j, k)[1] = flowField.getVelocity().getVector(i, j - 1, k)[1];
  // flowField.getVelocity().getVector(i, j, k)[2] = flowField.getVelocity().getVector(i, j - 1, k)[2];
}

void Stencils::NeumannKBoundaryStencil::applyFrontWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
  flowField.getk().getScalar(i, j, k) = flowField.getk().getScalar(i, j, k + 1);
  // flowField.getVelocity().getVector(i, j, k)[0]     = flowField.getVelocity().getVector(i, j, k + 1)[0];
  // flowField.getVelocity().getVector(i, j, k)[1]     = flowField.getVelocity().getVector(i, j, k + 1)[1];
  // flowField.getVelocity().getVector(i, j, k - 1)[2] = flowField.getVelocity().getVector(i, j, k)[2];
}

void Stencils::NeumannKBoundaryStencil::applyBackWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
  flowField.getk().getScalar(i, j, k) = flowField.getk().getScalar(i, j, k - 1);
  // flowField.getVelocity().getVector(i, j, k)[0] = flowField.getVelocity().getVector(i, j, k - 1)[0];
  // flowField.getVelocity().getVector(i, j, k)[1] = flowField.getVelocity().getVector(i, j, k - 1)[1];
  // flowField.getVelocity().getVector(i, j, k)[2] = flowField.getVelocity().getVector(i, j, k - 1)[2];
}

Stencils::NeumannEpsilonBoundaryStencil::NeumannEpsilonBoundaryStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowFieldKE>(parameters) {}

void Stencils::NeumannEpsilonBoundaryStencil::applyLeftWall(TurbulentFlowFieldKE& flowField, int i, int j) {
  flowField.geteps().getScalar(i, j) = flowField.geteps().getScalar(i + 1, j);
  // flowField.getVelocity().getVector(i, j)[1]     = flowField.getVelocity().getVector(i + 1, j)[1];
}

void Stencils::NeumannEpsilonBoundaryStencil::applyRightWall(TurbulentFlowFieldKE& flowField, int i, int j) {
  flowField.geteps().getScalar(i, j) = flowField.geteps().getScalar(i - 1, j);
  // flowField.getVelocity().getVector(i, j)[0] = flowField.getVelocity().getVector(i - 1, j)[0];
  // flowField.getVelocity().getVector(i, j)[1] = flowField.getVelocity().getVector(i - 1, j)[1];
}

void Stencils::NeumannEpsilonBoundaryStencil::applyBottomWall(TurbulentFlowFieldKE& flowField, int i, int j) {
  flowField.geteps().getScalar(i, j) = flowField.geteps().getScalar(i, j + 1);
  // flowField.getVelocity().getVector(i, j)[0]     = flowField.getVelocity().getVector(i, j + 1)[0];
  // flowField.getVelocity().getVector(i, j - 1)[1] = flowField.getVelocity().getVector(i, j)[1];
}

void Stencils::NeumannEpsilonBoundaryStencil::applyTopWall(TurbulentFlowFieldKE& flowField, int i, int j) {
  flowField.geteps().getScalar(i, j) = flowField.geteps().getScalar(i, j - 1);
  // flowField.getVelocity().getVector(i, j)[0] = flowField.getVelocity().getVector(i, j - 1)[0];
  // flowField.getVelocity().getVector(i, j)[1] = flowField.getVelocity().getVector(i, j - 1)[1];
}

void Stencils::NeumannEpsilonBoundaryStencil::applyLeftWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
  flowField.geteps().getScalar(i, j, k) = flowField.geteps().getScalar(i + 1, j, k);
  // flowField.getVelocity().getVector(i - 1, j, k)[0] = flowField.getVelocity().getVector(i, j, k)[0];
  // flowField.getVelocity().getVector(i, j, k)[1]     = flowField.getVelocity().getVector(i + 1, j, k)[1];
  // flowField.getVelocity().getVector(i, j, k)[2]     = flowField.getVelocity().getVector(i + 1, j, k)[2];
}

void Stencils::NeumannEpsilonBoundaryStencil::applyRightWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
  flowField.geteps().getScalar(i, j, k) = flowField.geteps().getScalar(i - 1, j, k);
  // flowField.getVelocity().getVector(i, j, k)[0] = flowField.getVelocity().getVector(i - 1, j, k)[0];
  // flowField.getVelocity().getVector(i, j, k)[1] = flowField.getVelocity().getVector(i - 1, j, k)[1];
  // flowField.getVelocity().getVector(i, j, k)[2] = flowField.getVelocity().getVector(i - 1, j, k)[2];
}

void Stencils::NeumannEpsilonBoundaryStencil::applyBottomWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
  flowField.geteps().getScalar(i, j, k) = flowField.geteps().getScalar(i, j + 1, k);
  // flowField.getVelocity().getVector(i, j, k)[0]     = flowField.getVelocity().getVector(i, j + 1, k)[0];
  // flowField.getVelocity().getVector(i, j - 1, k)[1] = flowField.getVelocity().getVector(i, j, k)[1];
  // flowField.getVelocity().getVector(i, j, k)[2]     = flowField.getVelocity().getVector(i, j + 1, k)[2];
}

void Stencils::NeumannEpsilonBoundaryStencil::applyTopWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
  flowField.geteps().getScalar(i, j, k) = flowField.geteps().getScalar(i, j - 1, k);
  // flowField.getVelocity().getVector(i, j, k)[0] = flowField.getVelocity().getVector(i, j - 1, k)[0];
  // flowField.getVelocity().getVector(i, j, k)[1] = flowField.getVelocity().getVector(i, j - 1, k)[1];
  // flowField.getVelocity().getVector(i, j, k)[2] = flowField.getVelocity().getVector(i, j - 1, k)[2];
}

void Stencils::NeumannEpsilonBoundaryStencil::applyFrontWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
  flowField.geteps().getScalar(i, j, k) = flowField.geteps().getScalar(i, j, k + 1);
  // flowField.getVelocity().getVector(i, j, k)[0]     = flowField.getVelocity().getVector(i, j, k + 1)[0];
  // flowField.getVelocity().getVector(i, j, k)[1]     = flowField.getVelocity().getVector(i, j, k + 1)[1];
  // flowField.getVelocity().getVector(i, j, k - 1)[2] = flowField.getVelocity().getVector(i, j, k)[2];
}

void Stencils::NeumannEpsilonBoundaryStencil::applyBackWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
  flowField.geteps().getScalar(i, j, k) = flowField.geteps().getScalar(i, j, k - 1);
  // flowField.getVelocity().getVector(i, j, k)[0] = flowField.getVelocity().getVector(i, j, k - 1)[0];
  // flowField.getVelocity().getVector(i, j, k)[1] = flowField.getVelocity().getVector(i, j, k - 1)[1];
  // flowField.getVelocity().getVector(i, j, k)[2] = flowField.getVelocity().getVector(i, j, k - 1)[2];
}