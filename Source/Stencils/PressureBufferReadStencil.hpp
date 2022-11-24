#pragma once

#include "BoundaryStencil.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Iterators.hpp"
#include "Parameters.hpp"

namespace Stencils {

  class PressureBufferReadStencil: public BoundaryStencil<FlowField> {
  private:
    RealType* bLeft;
    RealType* bRight;
    RealType* bTop;
    RealType* bBottom;
    RealType* bFront;
    RealType* bBack;

  public:
    PressureBufferReadStencil(const Parameters&, RealType*, RealType*, RealType*, RealType*, RealType*, RealType*);
    ~PressureBufferReadStencil() override = default;

    void applyLeftWall(FlowField& flowField, int i, int j) override;
    void applyRightWall(FlowField& flowField, int i, int j) override;
    void applyBottomWall(FlowField& flowField, int i, int j) override;
    void applyTopWall(FlowField& flowField, int i, int j) override;

    void applyLeftWall(FlowField& flowField, int i, int j, int k) override;
    void applyRightWall(FlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(FlowField& flowField, int i, int j, int k) override;
    void applyTopWall(FlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(FlowField& flowField, int i, int j, int k) override;
    void applyBackWall(FlowField& flowField, int i, int j, int k) override;
  };
} // namespace Stencils