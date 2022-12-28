#include "TurbulentSimulation.hpp"

TurbulentSimulation::TurbulentSimulation(Parameters& parameters, TurbulentFlowField& flowField):
  Simulation(parameters, flowField),
  turbflowField_(flowField),
  TurbulentFGHStencil_(parameters),
  TurbulentFGHIterator_(turbflowField_, parameters, TurbulentFGHStencil_),
  nuTStencil_(parameters),
  nuTKEStencil_(parameters),
  nuTIterator_(turbflowField_, parameters, nuTStencil_, 0, 0),
  hStencil_(parameters),
  hIterator_(turbflowField_, parameters, hStencil_, 0, 0),
  dtStencil_(parameters),
  dtIterator_(turbflowField_, parameters, dtStencil_),
  ppmTurbulent_(parameters, turbflowField_) {}

void TurbulentSimulation::initializeFlowField() {
  Simulation::initializeFlowField();
  hUpdate();
  // nuTUpdate();
}

void TurbulentSimulation::solveTimestep() {
  nuTUpdate();
  // std::cout << "***************************************************************************\n";
  // turbflowField_.getnuT().show();
  // Communicate viscosity

  ppmTurbulent_.communicateViscosity();
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
  ppmTurbulent_.communicatePressure();
  // Compute velocity
  velocityIterator_.iterate();
  obstacleIterator_.iterate();
  // TODO WS2: communicate velocity values
  ppmTurbulent_.communicateVelocities();
  // Iterate for velocities on the boundary
  wallVelocityIterator_.iterate();

  // Calculate viscosity
}

void TurbulentSimulation::setTimeStep() {
  RealType localMin, globalMin;
  ASSERTION(parameters_.geometry.dim == 2 || parameters_.geometry.dim == 3);

  // Determine maximum velocity
  maxUStencil_.reset();
  maxUFieldIterator_.iterate();
  maxUBoundaryIterator_.iterate();
  dtStencil_.reset();
  dtIterator_.iterate();

  localMin = dtStencil_.getDt();

  if (parameters_.geometry.dim == 3) {
    parameters_.timestep.dt = 1.0 / (maxUStencil_.getMaxValues()[2] + std::numeric_limits<double>::min());
  } else {
    parameters_.timestep.dt = 1.0 / (maxUStencil_.getMaxValues()[0] + std::numeric_limits<double>::min());
  }

  localMin = std::min(
    localMin,
    std::min(
      parameters_.timestep.dt,
      std::min(
        1 / (maxUStencil_.getMaxValues()[1] + std::numeric_limits<double>::min()),
        1 / (maxUStencil_.getMaxValues()[0] + std::numeric_limits<double>::min())
      )
    )
  );
  // if (fetestexcept(FE_DIVBYZERO))
  //     std::cout <<"Exception occured\n";

  // Here, we select the type of operation before compiling. This allows to use the correct
  // data type for MPI. Not a concern for small simulations, but useful if using heterogeneous
  // machines.

  globalMin = MY_FLOAT_MAX;
  MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN, PETSC_COMM_WORLD);

  parameters_.timestep.dt = globalMin;
  parameters_.timestep.dt *= parameters_.timestep.tau;
}

void TurbulentSimulation::hUpdate() { hIterator_.iterate(); }

void TurbulentSimulation::nuTUpdate() { nuTIterator_.iterate(); }

void TurbulentSimulation::plotVTK(int timeStep, RealType simulationTime) {
  Stencils::TurbulentVTKStencil     TurbulentvtkStencil(parameters_);
  FieldIterator<TurbulentFlowField> TurbulentvtkIterator(turbflowField_, parameters_, TurbulentvtkStencil, 1, 0);

  TurbulentvtkIterator.iterate();
  TurbulentvtkStencil.write(turbflowField_, timeStep, simulationTime);
}