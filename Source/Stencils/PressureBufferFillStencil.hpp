#pragma once

#include "BoundaryStencil.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "Iterators.hpp"

namespace Stencils {

  class PressureBufferFillStencil: public BoundaryStencil<FlowField> {
  private:
    std::vector<RealType> bLeft;
    std::vector<RealType> bRight;
    std::vector<RealType> bTop;
    std::vector<RealType> bBottom;
    std::vector<RealType> bFront;
    std::vector<RealType> bBack;
  public:
    PressureBufferFillStencil(const Parameters& parameters);
    ~PressureBufferFillStencil() override = default;

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