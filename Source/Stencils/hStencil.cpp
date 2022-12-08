#include "StdAfx.hpp"

#include "hStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

Stencils::hStencil::hStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) {}

void Stencils::hStencil::apply(TurbulentFlowField& flowField, int i, int j) {
  auto xPos = parameters_.meshsize->getPosX(i, j) + 0.5 * parameters_.meshsize->getDx(i, j);
  auto yPos = parameters_.meshsize->getPosY(i, j) + 0.5 * parameters_.meshsize->getDy(i, j);

  // For CHannel Flow
  if (parameters_.bfStep.xRatio < 0 || parameters_.bfStep.yRatio < 0) {
    flowField.getheight().getScalar(i, j) = std::min(yPos, (parameters_.geometry.lengthY - yPos));

    // ****Remove x pos for cells near to inlet****
  }
  // For BFS
  else {
    if (xPos <= parameters_.bfStep.xRatio * parameters_.geometry.lengthX) {
      flowField.getheight().getScalar(i, j) = std::min(
        yPos - (parameters_.bfStep.yRatio * parameters_.geometry.lengthY), (parameters_.geometry.lengthY - yPos)
      );
    } else if (yPos < parameters_.bfStep.yRatio * parameters_.geometry.lengthY) {
      flowField.getheight().getScalar(i, j) = std::min(
        xPos - (parameters_.bfStep.xRatio * parameters_.geometry.lengthX),
        std::min((parameters_.geometry.lengthX - xPos), std::min(yPos, (parameters_.geometry.lengthY - yPos)))
      );
    } else {
      flowField.getheight().getScalar(i, j) = std::min(yPos, (parameters_.geometry.lengthY - yPos));
    }
  }
}

void Stencils::hStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
  auto xPos = parameters_.meshsize->getPosX(i, j, k) + 0.5 * parameters_.meshsize->getDx(i, j, k);
  auto yPos = parameters_.meshsize->getPosY(i, j, k) + 0.5 * parameters_.meshsize->getDy(i, j, k);
  auto zPos = parameters_.meshsize->getPosZ(i, j, k) + 0.5 * parameters_.meshsize->getDz(i, j, k);

  if (parameters_.bfStep.xRatio < 0 || parameters_.bfStep.yRatio < 0) {
    flowField.getheight().getScalar(i, j, k) = std::min(
      yPos, std::min((parameters_.geometry.lengthY - yPos), std::min(zPos, parameters_.geometry.lengthZ - zPos))
    );
    // ****Remove x pos for cells near to inlet****
  } else {
    if (xPos <= parameters_.bfStep.xRatio * parameters_.geometry.lengthX) {
      flowField.getheight().getScalar(i, j, k) = std::min(
        yPos - (parameters_.bfStep.yRatio * parameters_.geometry.lengthY),
        std::min((parameters_.geometry.lengthY - yPos), std::min(zPos, (parameters_.geometry.lengthZ - zPos)))
      );
    } else {
      flowField.getheight().getScalar(i, j, k) = std::min(
        yPos, std::min((parameters_.geometry.lengthY - yPos), std::min(zPos, (parameters_.geometry.lengthZ - zPos)))
      );
    }
  }
}
