#include "StdAfx.hpp"

#include "timeStepStencilKE.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::timeStepStencilKE::timeStepStencilKE(const Parameters& parameters):
  FieldStencil<TurbulentFlowFieldKE>(parameters) {
  Mindt = parameters.timestep.dt;
}

void Stencils::timeStepStencilKE::apply(TurbulentFlowFieldKE& flowField, int i, int j) {
  RealType factor
    = 1.0
        / (FieldStencil<TurbulentFlowFieldKE>::parameters_.meshsize->getDxMin() * FieldStencil<TurbulentFlowFieldKE>::parameters_.meshsize->getDxMin())
      + 1.0
          / (FieldStencil<TurbulentFlowFieldKE>::parameters_.meshsize->getDyMin() * FieldStencil<TurbulentFlowFieldKE>::parameters_.meshsize->getDyMin());

  RealType vtotal = 1 / FieldStencil<TurbulentFlowFieldKE>::parameters_.flow.Re + flowField.getnuT().getScalar(i, j);
  Mindt           = std::min(1 / (2 * vtotal * factor), Mindt);
}

void Stencils::timeStepStencilKE::apply(TurbulentFlowFieldKE& flowField, int i, int j, int k) {
  RealType factor
    = 1.0
        / (FieldStencil<TurbulentFlowFieldKE>::parameters_.meshsize->getDxMin() * FieldStencil<TurbulentFlowFieldKE>::parameters_.meshsize->getDxMin())
      + 1.0
          / (FieldStencil<TurbulentFlowFieldKE>::parameters_.meshsize->getDyMin() * FieldStencil<TurbulentFlowFieldKE>::parameters_.meshsize->getDyMin())
      + 1.0
          / (FieldStencil<TurbulentFlowFieldKE>::parameters_.meshsize->getDzMin() * FieldStencil<TurbulentFlowFieldKE>::parameters_.meshsize->getDzMin());

  Mindt = std::min(
    (FieldStencil<TurbulentFlowFieldKE>::parameters_.flow.Re + 1 / flowField.getnuT().getScalar(i, j, k)
    ) / (2 * factor),
    Mindt
  );
}

void Stencils::timeStepStencilKE::reset() { Mindt = MY_FLOAT_MAX; }

RealType Stencils::timeStepStencilKE::getDt() { return Mindt; }
