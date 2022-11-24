#include "StdAfx.hpp"

#include "VelocityBufferFillStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::VelocityBufferFillStencil::VelocityBufferFillStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters) {
}

void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j){
  bLeft[j] = flowField.getVelocity().getVector(i+2,j)[1];
}

void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j){
  bRight[j] = flowField.getVelocity().getVector(i-1,j)[1];
}
void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j){
  bBottom[i] = flowField.getVelocity().getVector(i,j+2)[1];
}
void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j){
  bTop[i] = flowField.getVelocity().getVector(i,j-1)[1];
}


void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k){
  bLeft[j * parameters_.geometry.sizeZ + k] = flowField.getVelocity().getVector(i+2,j, k)[1];
}
void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k){
  bRight[j * parameters_.geometry.sizeZ + k] = flowField.getVelocity().getVector(i-1,j, k)[1];
}
void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k){
  bBottom[i * parameters_.geometry.sizeZ + k] = flowField.getVelocity().getVector(i,j+2)[1];
}
void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k){
  bTop[i * parameters_.geometry.sizeZ + k] = flowField.getVelocity().getVector(i,j-1, k)[1];
}
void Stencils::VelocityBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k){
  bFront[i * parameters_.geometry.sizeY + j] = flowField.getVelocity().getVector(i,j, k+2)[1];
}
void Stencils::VelocityBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k){
  bBack[i * parameters_.geometry.sizeY + j] = flowField.getVelocity().getVector(i,j, k-1)[1];
}

