#pragma once
#include "Simulation.hpp"
#include "TurbulentFlowFieldKE.hpp"

#include "ParallelManagers/PetscParallelManager.hpp"
#include "Stencils/hStencilKE.hpp"
#include "Stencils/MaxUStencil.hpp"
#include "Stencils/nuTStencilKE.hpp"
#include "Stencils/timeStepStencilKE.hpp"
#include "Stencils/TurbulentVTKStencilKE.hpp"

class TurbulentSimulationKE: public Simulation {
protected:
  TurbulentFlowFieldKE& turbflowFieldKE_;

  // Stencils::MaxUStencil             maxUStencil_;
  // FieldIterator<FlowField>          maxUFieldIterator_;
  // GlobalBoundaryIterator<TurbulentFlowFieldKE> maxUBoundaryIterator_;

  // // Set up the boundary conditions
  // GlobalBoundaryFactory             globalBoundaryFactory_;
  // GlobalBoundaryIterator<FlowField> wallVelocityIterator_;
  // GlobalBoundaryIterator<FlowField> wallFGHIterator_;

  Stencils::TurbulentFGHStencilKE     TurbulentFGHStencilKE_;
  FieldIterator<TurbulentFlowFieldKE> TurbulentFGHIteratorKE_;

  Stencils::nuTStencilKE              nuTStencilKE_;
  FieldIterator<TurbulentFlowFieldKE> nuTIteratorKE_;

  Stencils::hStencilKE                hStencilKE_;
  FieldIterator<TurbulentFlowFieldKE> hIteratorKE_;

  // Stencils::kStencil                  kStencil_;
  // FieldIterator<TurbulentFlowFieldKE> kIterator_;

  // Stencils::eStencil                  eStencil_;
  // FieldIterator<TurbulentFlowFieldKE> eIterator_;

  Stencils::timeStepStencilKE         dtStencil_;
  FieldIterator<TurbulentFlowFieldKE> dtIterator_;

  //   Stencils::VelocityStencil velocityStencil_;
  //   Stencils::ObstacleStencil obstacleStencil_;
  //   FieldIterator<FlowField>  velocityIterator_;
  //   FieldIterator<FlowField>  obstacleIterator_;

  //   Stencils::RHSStencil     rhsStencil_;
  //   FieldIterator<FlowField> rhsIterator_;

  //   std::unique_ptr<Solvers::LinearSolver> solver_;
  ParallelManagers::PetscParallelManager<TurbulentFlowFieldKE> ppmTurbulentKE_;

  virtual void setTimeStep() override;

public:
  TurbulentSimulationKE(Parameters&, TurbulentFlowFieldKE&);
  virtual ~TurbulentSimulationKE() override = default;

  /** Initialises the flow field according to the scenario */
  virtual void initializeFlowField() override;

  virtual void solveTimestep() override;

  virtual void hUpdate() override;

  virtual void nuTUpdate() override;

  // virtual void keUpdate() override;

  // virtual void eUpdate() override;

  void plotVTK(int, RealType) override;

  /** Plots the flow field */
  // virtual void plotVTK(int timeStep, RealType simulationTime);
};