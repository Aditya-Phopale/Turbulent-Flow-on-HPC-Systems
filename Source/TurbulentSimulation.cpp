#include "TurbulentSimulation.hpp"

TurbulentSimulation::TurbulentSimulation(Parameters& parameters, TurbulentFlowField& flowField):
  Simulation(parameters, flowField),
  turbflowField_(flowField),
  TurbulentFGHStencil_(parameters),
  TurbulentFGHIterator_(flowField, parameters, TurbulentFGHStencil_),
  nuTStencil_(parameters),
  nuTIterator_(turbflowField_, parameters, nuTStencil_, 1, 0),
  hStencil_(parameters),
  hIterator_(turbflowField_, parameters, hStencil_, 1, 0) {}

void TurbulentSimulation::initializeFlowField() {
  Simulation::initializeFlowField();
  hUpdate();
  nuTUpdate();
}

void TurbulentSimulation::solveTimestep() {

  // Communicate viscosity

  // Determine and set max. timestep which is allowed in this simulation
  setTimeStep();
  // Compute FGH
  TurbulentFGHIterator_.iterate();
  // Set global boundary values
  wallFGHIterator_.iterate();
  // TODO WS1: compute the right hand side (RHS)
  rhsIterator_.iterate();
  // Solve for pressure
  solver_->solve();
  // TODO WS2: communicate pressure values
  // Compute velocity
  velocityIterator_.iterate();
  obstacleIterator_.iterate();
  // TODO WS2: communicate velocity values
  // Iterate for velocities on the boundary
  wallVelocityIterator_.iterate();

  // Calculate viscosity
}

void TurbulentSimulation::setTimeStep() {} //{ tIterator_.iterate(); }

void TurbulentSimulation::hUpdate() { hIterator_.iterate(); }

void TurbulentSimulation::nuTUpdate() { nuTIterator_.iterate(); }