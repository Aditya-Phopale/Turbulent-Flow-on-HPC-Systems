#include "StdAfx.hpp"

#include "MaxEpsStencil.hpp"

Stencils::MaxEpsStencil::MaxEpsStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowFieldKE>(parameters) {

  reset();
}

void Stencils::MaxEpsStencil::apply(TurbulentFlowFieldKE& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxEpsStencil::apply(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
  cellMaxValue(flowField, i, j, k);
}

// void Stencils::MaxUStencil::applyLeftWall(TurbulentFlowFieldKE& flowField, int i, int j) { cellMaxValue(flowField, i,
// j); }

// void Stencils::MaxUStencil::applyRightWall(TurbulentFlowFieldKE& flowField, int i, int j) { cellMaxValue(flowField,
// i, j); }

// void Stencils::MaxUStencil::applyBottomWall(TurbulentFlowFieldKE& flowField, int i, int j) { cellMaxValue(flowField,
// i, j); }

// void Stencils::MaxUStencil::applyTopWall(TurbulentFlowFieldKE& flowField, int i, int j) { cellMaxValue(flowField, i,
// j); }

// void Stencils::MaxUStencil::applyLeftWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
//   cellMaxValue(flowField, i, j, k);
// }

// void Stencils::MaxUStencil::applyRightWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
//   cellMaxValue(flowField, i, j, k);
// }

// void Stencils::MaxUStencil::applyBottomWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
//   cellMaxValue(flowField, i, j, k);
// }

// void Stencils::MaxUStencil::applyTopWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
//   cellMaxValue(flowField, i, j, k);
// }

// void Stencils::MaxUStencil::applyFrontWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
//   cellMaxValue(flowField, i, j, k);
// }

// void Stencils::MaxUStencil::applyBackWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
//   cellMaxValue(flowField, i, j, k);
// }

void Stencils::MaxEpsStencil::cellMaxValue(TurbulentFlowFieldKE& flowField, int i, int j) {
  RealType       epsvalue = flowField.geteps().getScalar(i, j);
  const RealType dx     = FieldStencil<TurbulentFlowFieldKE>::parameters_.meshsize->getDx(i, j);
  const RealType dy     = FieldStencil<TurbulentFlowFieldKE>::parameters_.meshsize->getDy(i, j);
  if (fabs(epsvalue * (1 / (dx * dx) + 1 / (dy * dy))) > maxValues_) {
    maxValues_ = fabs(epsvalue * (1 / (dx * dx) + 1 / (dy * dy)));
  }
}

void Stencils::MaxEpsStencil::cellMaxValue(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
  RealType       epsvalue = flowField.getk().getScalar(i, j, k);
  const RealType dx     = FieldStencil<TurbulentFlowFieldKE>::parameters_.meshsize->getDx(i, j, k);
  const RealType dy     = FieldStencil<TurbulentFlowFieldKE>::parameters_.meshsize->getDy(i, j, k);
  const RealType dz     = FieldStencil<TurbulentFlowFieldKE>::parameters_.meshsize->getDz(i, j, k);
  if (fabs(epsvalue * (1 / (dx * dx) + 1 / (dy * dy) + 1 / (dz * dz))) > maxValues_) {
    maxValues_ = fabs(epsvalue * (1 / (dx * dx) + 1 / (dy * dy) + 1 / (dz * dz)));
  }
}

void Stencils::MaxEpsStencil::reset() { maxValues_ = 0; }

const RealType Stencils::MaxEpsStencil::getMaxValues() const { return maxValues_; }
