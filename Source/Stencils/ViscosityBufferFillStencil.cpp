#include "StdAfx.hpp"

#include "ViscosityBufferFillStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"


Stencils::ViscosityBufferFillStencil::
  ViscosityBufferFillStencil(const Parameters& parameters, std::vector<RealType>& left, std::vector<RealType>& right, std::vector<RealType>& top, std::vector<RealType>& bottom, std::vector<RealType>& front, std::vector<RealType>& back):
  BoundaryStencil<TurbulentFlowField>(parameters), Left_(left), Right_(right), Top_(top), Bottom_(bottom), Front_(front), Back_(back) {
}

void Stencils::ViscosityBufferFillStencil::applyLeftWall(TurbulentFlowField& turbFlowField, int i, int j) {
  Left_.at(j) = turbFlowField.getnuT().getScalar(i + 2, j);
}

void Stencils::ViscosityBufferFillStencil::applyRightWall(TurbulentFlowField& turbFlowField, int i, int j) {
  Right_.at(j) = turbFlowField.getnuT().getScalar(i - 1, j);
}
void Stencils::ViscosityBufferFillStencil::applyBottomWall(TurbulentFlowField& turbFlowField, int i, int j) {
  Bottom_.at(i) = turbFlowField.getnuT().getScalar(i, j + 2);
}
void Stencils::ViscosityBufferFillStencil::applyTopWall(TurbulentFlowField& turbFlowField, int i, int j) {
  Top_.at(i) = turbFlowField.getnuT().getScalar(i, j - 1);
}

void Stencils::ViscosityBufferFillStencil::applyLeftWall(TurbulentFlowField& turbFlowField, int i, int j, int k) {
  Left_.at(j * (parameters_.parallel.localSize[2] + 3) + k) = turbFlowField.getnuT().getScalar(i + 2, j, k);
}
void Stencils::ViscosityBufferFillStencil::applyRightWall(TurbulentFlowField& turbFlowField, int i, int j, int k) {
  Right_.at(j * (parameters_.parallel.localSize[2] + 3) + k) = turbFlowField.getnuT().getScalar(i - 1, j, k);
}
void Stencils::ViscosityBufferFillStencil::applyBottomWall(TurbulentFlowField& turbFlowField, int i, int j, int k) {
  Bottom_.at(i * (parameters_.parallel.localSize[2] + 3) + k) = (turbFlowField.getnuT().getScalar(i, j + 2, k));
}
void Stencils::ViscosityBufferFillStencil::applyTopWall(TurbulentFlowField& turbFlowField, int i, int j, int k) {
  Top_.at(i * (parameters_.parallel.localSize[2] + 3) + k) = turbFlowField.getnuT().getScalar(i, j - 1, k);
}
void Stencils::ViscosityBufferFillStencil::applyFrontWall(TurbulentFlowField& turbFlowField, int i, int j, int k) {
  Front_.at(i * (parameters_.parallel.localSize[1] + 3) + j) = turbFlowField.getnuT().getScalar(i, j, k + 2);
}
void Stencils::ViscosityBufferFillStencil::applyBackWall(TurbulentFlowField& turbFlowField, int i, int j, int k) {
  Back_.at(i * (parameters_.parallel.localSize[1] + 3) + j) = turbFlowField.getnuT().getScalar(i, j, k - 1);
}
