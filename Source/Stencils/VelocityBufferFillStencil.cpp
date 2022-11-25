#include "StdAfx.hpp"

#include "VelocityBufferFillStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::VelocityBufferFillStencil::VelocityBufferFillStencil(const Parameters& parameters, RealType* bL_v, RealType* bR_v, RealType* bBo_v, RealType* bT_v, RealType* bF_v, RealType* bB_v):
  BoundaryStencil<FlowField>(parameters),
  bLeft_v(bL_v),
  bRight_v(bR_v),
  bTop_v(bT_v),
  bBottom_v(bBo_v),
  bFront_v(bF_v),
  bBack_v(bB_v) {
}

void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j){
  bLeft_v[j] = flowField.getVelocity().getVector(i+2,j)[1];
}

void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j){
  bRight_v[j] = flowField.getVelocity().getVector(i-1,j)[1];
}
void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j){
  bBottom_v[i] = flowField.getVelocity().getVector(i,j+2)[1];
}
void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j){
  bTop_v[i] = flowField.getVelocity().getVector(i,j-1)[1];
}


void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k){
  bLeft_v[j * parameters_.geometry.sizeZ + k] = flowField.getVelocity().getVector(i+2,j, k)[1];
}
void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k){
  bRight_v[j * parameters_.geometry.sizeZ + k] = flowField.getVelocity().getVector(i-1,j, k)[1];
}
void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k){
  bBottom_v[i * parameters_.geometry.sizeZ + k] = flowField.getVelocity().getVector(i,j+2)[1];
}
void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k){
  bTop_v[i * parameters_.geometry.sizeZ + k] = flowField.getVelocity().getVector(i,j-1, k)[1];
}
void Stencils::VelocityBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k){
  bFront_v[i * parameters_.geometry.sizeY + j] = flowField.getVelocity().getVector(i,j, k+2)[1];
}
void Stencils::VelocityBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k){
  bBack_v[i * parameters_.geometry.sizeY + j] = flowField.getVelocity().getVector(i,j, k-1)[1];
}

