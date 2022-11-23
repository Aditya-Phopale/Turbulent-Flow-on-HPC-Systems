#include "StdAfx.hpp"

#include "PressureBufferFillStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::PressureBufferFillStencil::PressureBufferFillStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {

  pressurefill_iteratorx(parameters)
}

void Stencils::PressureBufferFillStencil::fill(FlowField& flowField) {
  int sizeX       = parameters_.parallel.localSize[0];
  int sizeY       = parameters_.parallel.localSize[1];
  int sizeZ       = parameters_.parallel.localSize[2];
  int rnk         = parameters_.parallel.rank;
  int lowOffsetX  = rnk * sizeX;
  int highOffsetX = lowOffsetX + sizeX;

  int lowOffsetY  = rnk * sizeY;
  int highOffsetY = lowOffsetY + sizeY;

  int lowOffsetZ  = rnk * sizeZ;
  int highOffsetZ = lowOffsetZ + sizeZ;

  auto pressure = flowField.getPressure();
  ParallelBoundaryIterator<FlowField> pressurefill_iteratorx(pressure, parameters_, );
  ParallelBoundaryIterator<FlowField> pressurefill_iteratory;
  ParallelBoundaryIterator<FlowField> pressurefill_iteratorz;
  if (parameters_.geometry.dim == 2)
    RealType* buffer[sizeX + sizeY];
  else
    RealType* buffer[sizeX + sizeY + sizeZ];
}