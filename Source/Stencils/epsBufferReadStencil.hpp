#pragma once

#include <vector>

#include "BoundaryStencil.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowFieldKE.hpp"

namespace Stencils {

  class epsBufferReadStencil: public BoundaryStencil<TurbulentFlowFieldKE> {
  private:
    std::vector<RealType> Left_;
    std::vector<RealType> Right_;
    std::vector<RealType> Top_;
    std::vector<RealType> Bottom_;
    std::vector<RealType> Front_;
    std::vector<RealType> Back_;

  public:
    epsBufferReadStencil(const Parameters&, std::vector<RealType>&, std::vector<RealType>&, std::vector<RealType>&, std::vector<RealType>&, std::vector<RealType>&, std::vector<RealType>&);
    ~epsBufferReadStencil() override = default;

    void applyLeftWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j) override;
    void applyRightWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j) override;
    void applyBottomWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j) override;
    void applyTopWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j) override;

    void applyLeftWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) override;
    void applyRightWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) override;
    void applyBottomWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) override;
    void applyTopWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) override;
    void applyFrontWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) override;
    void applyBackWall(TurbulentFlowFieldKE& turbFlowFieldKE, int i, int j, int k) override;
  };
} // namespace Stencils