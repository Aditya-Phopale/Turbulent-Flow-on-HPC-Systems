#include "StdAfx.hpp"

#include "hStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::hStencil::hStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {}

void Stencils::hStencil::apply(FlowField& flowField, int i, int j) {
  auto xPos = parameters_.meshsize->getPosX(i, j);
  auto yPos = parameters_.meshsize->getPosY(i, j);

  if (parameters_.bfStep.xRatio < 0 || parameters_.bfStep.yRatio < 0) {
    flowField.geth(i, j) = std::min(
      xPos, (parameters_.geometry.lengthX - xPos), yPos, (parameters_.geometry.lengthY - yPos)
    );
    // ****Remove x pos for cells near to inlet****
  } else {
    if (xPos <= parameters_.bfStep.xRatio * parameters_.geometry.lengthX) {
      flowField.geth(i, j) = std::min(
        yPos - (parameters_.bfStep.yRatio * parameters_.geometry.lengthY), (parameters_.geometry.lengthY - yPos)
      );
    } else {
      flowField.geth(i, j) = std::min(yPos, (parameters_.geometry.lengthY - yPos));
    }
  }
}

void Stencils::hStencil::apply(FlowField& flowField, int i, int j, int k) {
  auto xPos = parameters_.meshsize->getPosX(i, j, k);
  auto yPos = parameters_.meshsize->getPosY(i, j, k);
  auto zPos = parameters_.meshsize->getPosZ(i, j, k);

  if (parameters_.bfStep.xRatio < 0 || parameters_.bfStep.yRatio < 0) {
    flowField.geth(i, j) = std::min(
      xPos, (parameters_.geometry.lengthX - xPos), yPos, (parameters_.geometry.lengthY - yPos), zPos, parameters_.geometry.lengthZ - zPos)
    );
    // ****Remove x pos for cells near to inlet****
  } else {
    if (xPos <= parameters_.bfStep.xRatio * parameters_.geometry.lengthX) {
      flowField.geth(i, j) = std::min(
        yPos - (parameters_.bfStep.yRatio * parameters_.geometry.lengthY), (parameters_.geometry.lengthY - yPos), zPos, parameters_.geometry.lengthZ - zPos)
      );
    } else {
      flowField.geth(i, j) = std::min(yPos, (parameters_.geometry.lengthY - yPos), zPos, parameters_.geometry.lengthZ - zPos));
    }
  }
}
