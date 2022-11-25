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
    RealType* bLeft_v;
    RealType* bRight_v;
    RealType* bTop_v;
    RealType* bBottom_v;
    RealType* bFront_v;
    RealType* bBack_v;
  public:
    VelocityBufferFillStencil(const Parameters& parameters, RealType*, RealType*, RealType*, RealType*, RealType*, RealType*);
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