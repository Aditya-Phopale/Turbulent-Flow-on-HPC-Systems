#pragma once

#include <vector>

#include "BoundaryStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {
  template <class FlowFieldType>
  class ViscosityBufferReadStencil: public BoundaryStencil<FlowFieldType> {
  private:
    std::vector<RealType> Left_;
    std::vector<RealType> Right_;
    std::vector<RealType> Top_;
    std::vector<RealType> Bottom_;
    std::vector<RealType> Front_;
    std::vector<RealType> Back_;

  public:
    ViscosityBufferReadStencil(const Parameters&, std::vector<RealType>&, std::vector<RealType>&, std::vector<RealType>&, std::vector<RealType>&, std::vector<RealType>&, std::vector<RealType>&);
    ~ViscosityBufferReadStencil() override = default;

    void applyLeftWall(FlowFieldType& FlowField, int i, int j) override;
    void applyRightWall(FlowFieldType& FlowField, int i, int j) override;
    void applyBottomWall(FlowFieldType& FlowField, int i, int j) override;
    void applyTopWall(FlowFieldType& FlowField, int i, int j) override;

    void applyLeftWall(FlowFieldType& FlowField, int i, int j, int k) override;
    void applyRightWall(FlowFieldType& FlowField, int i, int j, int k) override;
    void applyBottomWall(FlowFieldType& FlowField, int i, int j, int k) override;
    void applyTopWall(FlowFieldType& FlowField, int i, int j, int k) override;
    void applyFrontWall(FlowFieldType& FlowField, int i, int j, int k) override;
    void applyBackWall(FlowFieldType& FlowField, int i, int j, int k) override;
  };
} // namespace Stencils