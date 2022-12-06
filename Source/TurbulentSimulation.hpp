#include "hStencil.hpp"
#include "nuTStencil.hpp"
#include "Simulation.hpp"

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

  Stencils::FGHStencil     TurbulentFGHStencil_;
  FieldIterator<FlowField> TurbulentFGHIterator_;

  Stencils::nuTStencil     nuTStencil_;
  FieldIterator<FlowField> nuTIterator_;

  Stencils::hStencil       hStecnil_;
  FieldIterator<FlowField> hIterator_;

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

  void hUpdate();

  /** Plots the flow field */
  virtual void plotVTK(int timeStep, RealType simulationTime);
};