#pragma once

#include <vector>
#include "BoundaryStencil.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "Iterators.hpp"

namespace Stencils {

  class VelocityBufferFillStencil: public BoundaryStencil<FlowField> {
  private:
    RealType* bLeft;
    RealType* bRight;
    RealType* bTop;
    RealType* bBottom;
    RealType* bFront;
    RealType* bBack;
  public:
    VelocityBufferFillStencil(const Parameters& parameters);
    ~VelocityBufferFillStencil() override = default;

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