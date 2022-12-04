#include "Simulation.hpp"
#include "ViscStencil.hpp"

class TurbulentSimulation: public Simulation {
protected:
  //   FlowField& flowField_;

  //   Stencils::MaxUStencil             maxUStencil_;
  //   FieldIterator<FlowField>          maxUFieldIterator_;
  //   GlobalBoundaryIterator<FlowField> maxUBoundaryIterator_;

  //   // Set up the boundary conditions
  //   GlobalBoundaryFactory             globalBoundaryFactory_;
  //   GlobalBoundaryIterator<FlowField> wallVelocityIterator_;
  //   GlobalBoundaryIterator<FlowField> wallFGHIterator_;

  Stencils::FGHStencil     fghStencil_;
  FieldIterator<FlowField> fghIterator_;

  Stencils::ViscStencil    viscStencil_;
  FieldIterator<FlowField> viscIterator_;

  //   Stencils::VelocityStencil velocityStencil_;
  //   Stencils::ObstacleStencil obstacleStencil_;
  //   FieldIterator<FlowField>  velocityIterator_;
  //   FieldIterator<FlowField>  obstacleIterator_;

  //   Stencils::RHSStencil     rhsStencil_;
  //   FieldIterator<FlowField> rhsIterator_;

  //   std::unique_ptr<Solvers::LinearSolver> solver_;

  virtual void setTimeStep();

public:
  TurbulentSimulation(Parameters&, FlowField&);
  virtual ~TurbulentSimulation() = default;

  /** Initialises the flow field according to the scenario */
  virtual void initializeFlowField();

  virtual void solveTimestep();

  /** Plots the flow field */
  virtual void plotVTK(int timeStep, RealType simulationTime);
};