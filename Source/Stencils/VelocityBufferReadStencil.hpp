#pragma once

#include <vector>

#include "BoundaryStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  template <class FlowFieldType>
  class VelocityBufferReadStencil: public BoundaryStencil<FlowFieldType> {
  private:
    std::vector<RealType> Left_;
    std::vector<RealType> Right_;
    std::vector<RealType> Top_;
    std::vector<RealType> Bottom_;
    std::vector<RealType> Front_;
    std::vector<RealType> Back_;

  public:
    VelocityBufferReadStencil(const Parameters&, std::vector<RealType>, std::vector<RealType>, std::vector<RealType>, std::vector<RealType>);
    VelocityBufferReadStencil(const Parameters&, std::vector<RealType>, std::vector<RealType>, std::vector<RealType>, std::vector<RealType>, std::vector<RealType>, std::vector<RealType>);

    ~VelocityBufferReadStencil() override = default;

    void applyLeftWall(FlowFieldType& flowField, int i, int j) override;
    void applyRightWall(FlowFieldType& flowField, int i, int j) override;
    void applyBottomWall(FlowFieldType& flowField, int i, int j) override;
    void applyTopWall(FlowFieldType& flowField, int i, int j) override;

    void applyLeftWall(FlowFieldType& flowField, int i, int j, int k) override;
    void applyRightWall(FlowFieldType& flowField, int i, int j, int k) override;
    void applyBottomWall(FlowFieldType& flowField, int i, int j, int k) override;
    void applyTopWall(FlowFieldType& flowField, int i, int j, int k) override;
    void applyFrontWall(FlowFieldType& flowField, int i, int j, int k) override;
    void applyBackWall(FlowFieldType& flowField, int i, int j, int k) override;
  };
} // namespace Stencils