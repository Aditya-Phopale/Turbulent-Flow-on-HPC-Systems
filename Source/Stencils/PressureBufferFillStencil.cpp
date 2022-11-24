#include "StdAfx.hpp"

#include "PressureBufferFillStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::PressureBufferFillStencil::PressureBufferFillStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters) {
}

void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j){
  bLeft.push_back(flowField.getPressure().getScalar(i+2,j));
}

void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j){
  bRight1.push_back(flowField.getPressure().getScalar(i-1,j));
  bRight2.push_back(flowField.getPressure().getScalar(i-2,j));
}
void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j){
  bBottom.push_back(flowField.getPressure().getScalar(i,j));
}
void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j){
  bTop.push_back(flowField.getPressure().getScalar(i,j));
}


void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k){
  bLeft.push_back(flowField.getPressure().getScalar(i,j));
}
void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k){
  bRight.push_back(flowField.getPressure().getScalar(i,j));
}
void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k){
  bBottom.push_back(flowField.getPressure().getScalar(i,j));
}
void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k){
  bTop.push_back(flowField.getPressure().getScalar(i,j));
}
void Stencils::PressureBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k){
  bFront.push_back(flowField.getPressure().getScalar(i,j));
}
void Stencils::PressureBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k){
  bBack.push_back(flowField.getPressure().getScalar(i,j));
}

