#include "StdAfx.hpp"

#include "timeStepStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::timeStepStencil::timeStepStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {
  Mindt = parameters.timestep.dt;
}

void Stencils::timeStepStencil::apply(FlowField& flowField, int i, int j) {
  RealType factor
    = 1.0
        / (FieldStencil<FlowField>::parameters_.meshsize->getDxMin() * FieldStencil<FlowField>::parameters_.meshsize->getDxMin())
      + 1.0
          / (FieldStencil<FlowField>::parameters_.meshsize->getDyMin() * FieldStencil<FlowField>::parameters_.meshsize->getDyMin());

  RealType vtotal = 1 / FieldStencil<FlowField>::parameters_.flow.Re + flowField.getnuT().getScalar(i, j);
  Mindt           = std::min(1 / (2 * vtotal * factor), Mindt);
}

void Stencils::timeStepStencil::apply(FlowField& flowField, int i, int j, int k) {
  RealType factor
    = 1.0
        / (FieldStencil<FlowField>::parameters_.meshsize->getDxMin() * FieldStencil<FlowField>::parameters_.meshsize->getDxMin())
      + 1.0
          / (FieldStencil<FlowField>::parameters_.meshsize->getDyMin() * FieldStencil<FlowField>::parameters_.meshsize->getDyMin())
      + 1.0
          / (FieldStencil<FlowField>::parameters_.meshsize->getDzMin() * FieldStencil<FlowField>::parameters_.meshsize->getDzMin());

  Mindt = std::min(
    (FieldStencil<FlowField>::parameters_.flow.Re + 1 / flowField.getnuT().getScalar(i, j, k)) / (2 * factor), Mindt
  );
}

void Stencils::timeStepStencil::reset() { Mindt = 100; } // temporary solution

RealType Stencils::timeStepStencil::getDt() { return Mindt; }
