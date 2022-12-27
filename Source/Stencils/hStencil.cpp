#include "StdAfx.hpp"

#include "hStencil.hpp"

#include "Definitions.hpp"
#include "StencilFunctions.hpp"

template <class FlowFieldType>
Stencils::hStencil<FlowFieldType>::hStencil(const Parameters& parameters):
  FieldStencil<FlowFieldType>(parameters) {}

template <class FlowFieldType>
void Stencils::hStencil<FlowFieldType>::apply(FlowFieldType& flowField, int i, int j) {
  auto xPos = FieldStencil<FlowFieldType>::parameters_.meshsize->getPosX(i, j)
              + 0.5 * FieldStencil<FlowFieldType>::parameters_.meshsize->getDx(i, j);
  auto yPos = FieldStencil<FlowFieldType>::parameters_.meshsize->getPosY(i, j)
              + 0.5 * FieldStencil<FlowFieldType>::parameters_.meshsize->getDy(i, j);

  // For CHannel Flow
  if (FieldStencil<FlowFieldType>::parameters_.bfStep.xRatio < 0 || FieldStencil<FlowFieldType>::parameters_.bfStep.yRatio < 0) {
    flowField.getheight().getScalar(i, j) = std::min(
      yPos, (FieldStencil<FlowFieldType>::parameters_.geometry.lengthY - yPos)
    );
  }
  // For BFS
  else {
    if (xPos <= FieldStencil<FlowFieldType>::parameters_.bfStep.xRatio * FieldStencil<FlowFieldType>::parameters_.geometry.lengthX) {
      flowField.getheight().getScalar(i, j) = std::min(
        yPos
          - (FieldStencil<FlowFieldType>::parameters_.bfStep.yRatio * FieldStencil<FlowFieldType>::parameters_.geometry.lengthY),
        (FieldStencil<FlowFieldType>::parameters_.geometry.lengthY - yPos)
      );
    } else if (yPos < FieldStencil<FlowFieldType>::parameters_.bfStep.yRatio * FieldStencil<FlowFieldType>::parameters_.geometry.lengthY) {
      flowField.getheight().getScalar(i, j) = std::min(
        xPos
          - (FieldStencil<FlowFieldType>::parameters_.bfStep.xRatio * FieldStencil<FlowFieldType>::parameters_.geometry.lengthX),
        std::min(yPos, (FieldStencil<FlowFieldType>::parameters_.geometry.lengthY - yPos))
      );
    } else {
      flowField.getheight().getScalar(i, j) = std::min(
        yPos, (FieldStencil<FlowFieldType>::parameters_.geometry.lengthY - yPos)
      );
    }
  }
}

template <class FlowFieldType>
void Stencils::hStencil<FlowFieldType>::apply(FlowFieldType& flowField, int i, int j, int k) {
  auto xPos = FieldStencil<FlowFieldType>::parameters_.meshsize->getPosX(i, j, k)
              + 0.5 * FieldStencil<FlowFieldType>::parameters_.meshsize->getDx(i, j, k);
  auto yPos = FieldStencil<FlowFieldType>::parameters_.meshsize->getPosY(i, j, k)
              + 0.5 * FieldStencil<FlowFieldType>::parameters_.meshsize->getDy(i, j, k);
  auto zPos = FieldStencil<FlowFieldType>::parameters_.meshsize->getPosZ(i, j, k)
              + 0.5 * FieldStencil<FlowFieldType>::parameters_.meshsize->getDz(i, j, k);

  if (FieldStencil<FlowFieldType>::parameters_.bfStep.xRatio < 0 || FieldStencil<FlowFieldType>::parameters_.bfStep.yRatio < 0) {
    flowField.getheight().getScalar(i, j, k) = std::min(
      yPos,
      std::min(
        (FieldStencil<FlowFieldType>::parameters_.geometry.lengthY - yPos),
        std::min(zPos, FieldStencil<FlowFieldType>::parameters_.geometry.lengthZ - zPos)
      )
    );
    // ****Remove x pos for cells near to inlet****
  } else {
    if (xPos <= FieldStencil<FlowFieldType>::parameters_.bfStep.xRatio * FieldStencil<FlowFieldType>::parameters_.geometry.lengthX) {
      flowField.getheight().getScalar(i, j, k) = std::min(
        yPos
          - (FieldStencil<FlowFieldType>::parameters_.bfStep.yRatio * FieldStencil<FlowFieldType>::parameters_.geometry.lengthY),
        std::min(
          (FieldStencil<FlowFieldType>::parameters_.geometry.lengthY - yPos),
          std::min(zPos, (FieldStencil<FlowFieldType>::parameters_.geometry.lengthZ - zPos))
        )
      );
    } else if (yPos < FieldStencil<FlowFieldType>::parameters_.bfStep.yRatio * FieldStencil<FlowFieldType>::parameters_.geometry.lengthY) {
      flowField.getheight().getScalar(i, j, k) = std::min(
        xPos
          - (FieldStencil<FlowFieldType>::parameters_.bfStep.xRatio * FieldStencil<FlowFieldType>::parameters_.geometry.lengthX),
        std::min(
          yPos,
          std::min(
            (FieldStencil<FlowFieldType>::parameters_.geometry.lengthY - yPos),
            std::min(zPos, (FieldStencil<FlowFieldType>::parameters_.geometry.lengthZ - zPos))
          )
        )
      );
    } else {
      flowField.getheight().getScalar(i, j) = std::min(
        yPos,
        std::min(
          (FieldStencil<FlowFieldType>::parameters_.geometry.lengthY - yPos),
          std::min(zPos, FieldStencil<FlowFieldType>::parameters_.geometry.lengthZ - zPos)
        )
      );
    }
  }
}
