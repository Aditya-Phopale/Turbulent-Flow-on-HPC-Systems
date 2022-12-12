#pragma once
#include "Simulation.hpp"
#include "TurbulentFlowField.hpp"

#include "ParallelManagers/PetscParallelManager.hpp"
#include "Stencils/hStencil.hpp"
#include "Stencils/nuTStencil.hpp"
#include "Stencils/timeStepStencil.hpp"
#include "Stencils/TurbulentVTKStencil.hpp"

class TurbulentSimulation: public Simulation {
protected:
  //   Stencils::MaxUStencil             maxUStencil_;
  //   FieldIterator<FlowField>          maxUFieldIterator_;
  //   GlobalBoundaryIterator<FlowField> maxUBoundaryIterator_;

  //   // Set up the boundary conditions
  //   GlobalBoundaryFactory             globalBoundaryFactory_;
  //   GlobalBoundaryIterator<FlowField> wallVelocityIterator_;
  //   GlobalBoundaryIterator<FlowField> wallFGHIterator_;

  Stencils::TurbulentFGHStencil TurbulentFGHStencil_;
  FieldIterator<FlowField>      TurbulentFGHIterator_;

  Stencils::nuTStencil     nuTStencil_;
  FieldIterator<FlowField> nuTIterator_;

  Stencils::hStencil       hStencil_;
  FieldIterator<FlowField> hIterator_;

  Stencils::timeStepStencil dtStencil_;
  FieldIterator<FlowField>  dtIterator_;

  //   Stencils::VelocityStencil velocityStencil_;
  //   Stencils::ObstacleStencil obstacleStencil_;
  //   FieldIterator<FlowField>  velocityIterator_;
  //   FieldIterator<FlowField>  obstacleIterator_;

  //   Stencils::RHSStencil     rhsStencil_;
  //   FieldIterator<FlowField> rhsIterator_;

  //   std::unique_ptr<Solvers::LinearSolver> solver_;
  // /ParallelManagers::PetscParallelManager<TurbulentFlowField> ppmTurbulent_;

  virtual void setTimeStep() override;

public:
  TurbulentSimulation(Parameters&, FlowField&);
  virtual ~TurbulentSimulation() override = default;

  /** Initialises the flow field according to the scenario */
  virtual void initializeFlowField() override;

  virtual void solveTimestep() override;

  virtual void hUpdate() override;

  virtual void nuTUpdate() override;

  void plotVTK(int, RealType) override;

  /** Plots the flow field */
  // virtual void plotVTK(int timeStep, RealType simulationTime);
};