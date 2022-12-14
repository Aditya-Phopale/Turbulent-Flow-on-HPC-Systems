#include "StdAfx.hpp"

#include "timeStepStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::timeStepStencil::timeStepStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) {
  Mindt = parameters.timestep.dt;
}

void Stencils::timeStepStencil::apply(TurbulentFlowField& flowField, int i, int j) {
  RealType factor
    = 1.0
        / (FieldStencil<TurbulentFlowField>::parameters_.meshsize->getDxMin() * FieldStencil<TurbulentFlowField>::parameters_.meshsize->getDxMin())
      + 1.0
          / (FieldStencil<TurbulentFlowField>::parameters_.meshsize->getDyMin() * FieldStencil<TurbulentFlowField>::parameters_.meshsize->getDyMin());

  RealType vtotal = 1 / FieldStencil<TurbulentFlowField>::parameters_.flow.Re + flowField.getnuT().getScalar(i, j);
  Mindt           = std::min(1 / (2 * vtotal * factor), Mindt);
}

void Stencils::timeStepStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
  RealType factor
    = 1.0
        / (FieldStencil<TurbulentFlowField>::parameters_.meshsize->getDxMin() * FieldStencil<TurbulentFlowField>::parameters_.meshsize->getDxMin())
      + 1.0
          / (FieldStencil<TurbulentFlowField>::parameters_.meshsize->getDyMin() * FieldStencil<TurbulentFlowField>::parameters_.meshsize->getDyMin())
      + 1.0
          / (FieldStencil<TurbulentFlowField>::parameters_.meshsize->getDzMin() * FieldStencil<TurbulentFlowField>::parameters_.meshsize->getDzMin());

  Mindt = std::min(
    (FieldStencil<TurbulentFlowField>::parameters_.flow.Re + 1 / flowField.getnuT().getScalar(i, j, k)) / (2 * factor),
    Mindt
  );
}

void Stencils::timeStepStencil::reset() { Mindt = MY_FLOAT_MAX; } 

RealType Stencils::timeStepStencil::getDt() { return Mindt; }
