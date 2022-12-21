#pragma once

#include <vector>

#include "BoundaryStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  class ViscosityBufferFillStencil: public BoundaryStencil<TurbulentFlowField> {
  private:
    std::vector<RealType>& Left_;
    std::vector<RealType>& Right_;
    std::vector<RealType>& Top_;
    std::vector<RealType>& Bottom_;
    std::vector<RealType>& Front_;
    std::vector<RealType>& Back_;

  public:
    ViscosityBufferFillStencil(const Parameters&, std::vector<RealType>&, std::vector<RealType>&, std::vector<RealType>&, std::vector<RealType>&, std::vector<RealType>&, std::vector<RealType>&);
    ~ViscosityBufferFillStencil() override = default;

    void applyLeftWall(TurbulentFlowField& turbFlowField, int i, int j) override;
    void applyRightWall(TurbulentFlowField& turbFlowField, int i, int j) override;
    void applyBottomWall(TurbulentFlowField& turbFlowField, int i, int j) override;
    void applyTopWall(TurbulentFlowField& turbFlowField, int i, int j) override;

    void applyLeftWall(TurbulentFlowField& turbFlowField, int i, int j, int k) override;
    void applyRightWall(TurbulentFlowField& turbFlowField, int i, int j, int k) override;
    void applyBottomWall(TurbulentFlowField& turbFlowField, int i, int j, int k) override;
    void applyTopWall(TurbulentFlowField& turbFlowField, int i, int j, int k) override;
    void applyFrontWall(TurbulentFlowField& turbFlowField, int i, int j, int k) override;
    void applyBackWall(TurbulentFlowField& turbFlowField, int i, int j, int k) override;

  };
} // namespace Stencils