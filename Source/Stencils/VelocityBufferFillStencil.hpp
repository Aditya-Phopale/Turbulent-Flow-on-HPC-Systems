#pragma once

#include <vector>

#include "BoundaryStencil.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Iterators.hpp"
#include "Parameters.hpp"

namespace Stencils {

  class VelocityBufferFillStencil: public BoundaryStencil<FlowField> {
  private:
    std::vector<RealType>& Left_;
    std::vector<RealType>& Right_;
    std::vector<RealType>& Top_;
    std::vector<RealType>& Bottom_;
    std::vector<RealType>& Front_;
    std::vector<RealType>& Back_;

  public:
    VelocityBufferFillStencil(const Parameters& parameters, std::vector<RealType>&, std::vector<RealType>&, std::vector<RealType>&, std::vector<RealType>&, std::vector<RealType>&, std::vector<RealType>&);
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