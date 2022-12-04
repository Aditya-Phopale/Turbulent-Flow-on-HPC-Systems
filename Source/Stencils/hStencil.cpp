#include "StdAfx.hpp"

#include "hStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::hStencil::hStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {}

void Stencils::hStencil::apply(FlowField& flowField, int i, int j) {
    auto xPos = parameters_.meshsize->getPosX(i,j);
    auto yPos = parameters_.meshsize->getPosY(i,j);

    if (parameters_.simulation.scenario == "channel")
    flowField.geth(i,j) = std::min(xPos,(parameters_.geometry.lengthX - xPos), yPos, (parameters_.geometry.lengthY));
}

void Stencils::hStencil::apply(FlowField& flowField, int i, int j, int k) {
    auto xPos = parameters_.meshsize->getPosX(i,j,k);
    auto yPos = parameters_.meshsize->getPosY(i,j,k);
    auto zPos = parameters_.meshsize->getPosZ(i,j,k);
}   


