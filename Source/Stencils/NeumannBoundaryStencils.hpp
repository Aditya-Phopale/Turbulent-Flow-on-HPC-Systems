#pragma once

#include "BoundaryStencil.hpp"
#include "DataStructures.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowFieldKE.hpp"

namespace Stencils {

  class NeumannVelocityBoundaryStencil: public BoundaryStencil<FlowField> {
  public:
    NeumannVelocityBoundaryStencil(const Parameters& parameters);
    ~NeumannVelocityBoundaryStencil() override = default;

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

  class NeumannFGHBoundaryStencil: public BoundaryStencil<FlowField> {
  public:
    NeumannFGHBoundaryStencil(const Parameters& parameters);
    ~NeumannFGHBoundaryStencil() override = default;

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

  class NeumannKBoundaryStencil: public BoundaryStencil<TurbulentFlowFieldKE> {
  public:
    NeumannKBoundaryStencil(const Parameters& parameters);
    ~NeumannKBoundaryStencil() override = default;

    void applyLeftWall(TurbulentFlowFieldKE& flowField, int i, int j) override;
    void applyRightWall(TurbulentFlowFieldKE& flowField, int i, int j) override;
    void applyBottomWall(TurbulentFlowFieldKE& flowField, int i, int j) override;
    void applyTopWall(TurbulentFlowFieldKE& flowField, int i, int j) override;

    void applyLeftWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
    void applyRightWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
    void applyBottomWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
    void applyTopWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
    void applyFrontWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
    void applyBackWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
  };

  class NeumannEpsilonBoundaryStencil: public BoundaryStencil<TurbulentFlowFieldKE> {
  public:
    NeumannEpsilonBoundaryStencil(const Parameters& parameters);
    ~NeumannEpsilonBoundaryStencil() override = default;

    void applyLeftWall(TurbulentFlowFieldKE& flowField, int i, int j) override;
    void applyRightWall(TurbulentFlowFieldKE& flowField, int i, int j) override;
    void applyBottomWall(TurbulentFlowFieldKE& flowField, int i, int j) override;
    void applyTopWall(TurbulentFlowFieldKE& flowField, int i, int j) override;

    void applyLeftWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
    void applyRightWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
    void applyBottomWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
    void applyTopWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
    void applyFrontWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
    void applyBackWall(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
