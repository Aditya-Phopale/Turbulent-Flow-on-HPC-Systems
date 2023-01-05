#include "TurbulentSimulationKE.hpp"

TurbulentSimulationKE::TurbulentSimulationKE(Parameters& parameters, TurbulentFlowFieldKE& flowField):
  Simulation(parameters, flowField),
  turbflowFieldKE_(flowField),
  TurbulentFGHStencilKE_(parameters),
  TurbulentFGHIteratorKE_(turbflowFieldKE_, parameters, TurbulentFGHStencilKE_),
  nuTStencilKE_(parameters),
  nuTIteratorKE_(turbflowFieldKE_, parameters, nuTStencilKE_, 1, 0),
  globalTurbulentBoundaryFactory_(parameters),
  wallkIterator_(globalTurbulentBoundaryFactory_.getGlobalBoundaryKIterator(turbflowFieldKE_)),
  wallEpsilonIterator_(globalTurbulentBoundaryFactory_.getGlobalBoundaryEpsilonIterator(turbflowFieldKE_)),
  hStencil_(parameters),
  hIterator_(turbflowFieldKE_, parameters, hStencil_, 0, 0),
  dtStencil_(parameters),
  dtIterator_(turbflowFieldKE_, parameters, dtStencil_),
  kStencil_(parameters),
  kIterator_(flowField, parameters, kStencil_, 1, 0),
  epsilonStencil_(parameters),
  epsilonIterator_(flowField, parameters, epsilonStencil_, 1, 0),
  ppmTurbulentKE_(parameters, turbflowFieldKE_) {}

void TurbulentSimulationKE::initializeFlowField() {
  Simulation::initializeFlowField();
  hUpdate();
  if (parameters_.simulation.scenario == "channel") {
    // Currently, a particular initialisation is only required for the taylor-green vortex.
    Stencils::InitkFlowFieldStencil     kInitStencil(parameters_);
    FieldIterator<TurbulentFlowFieldKE> kInitIterator(turbflowFieldKE_, parameters_, kInitStencil);
    kInitIterator.iterate();

    Stencils::InitEpsilonFlowFieldStencil epsInitStencil(parameters_);
    FieldIterator<TurbulentFlowFieldKE>   epsInitIterator(turbflowFieldKE_, parameters_, epsInitStencil);
    epsInitIterator.iterate();
  }
  // nuTUpdate();
}

void TurbulentSimulationKE::solveTimestep() {
  // nuTUpdate();
  // std::cout << "***************************************************************************\n";
  // turbflowFieldKE_.getnuT().show();
  // Communicate viscosity

  // ppmTurbulentKE_.communicateViscosity();
  // Determine and set max. timestep which is allowed in this simulation
  setTimeStep();
  std::cout << parameters_.timestep.dt << "\n";
  // parameters_.timestep.dt = 1e-4;
  // std::cout << "***************************************************************************\n";
  // turbflowFieldKE_.getk().show();
  wallkIterator_.iterate();
  wallEpsilonIterator_.iterate();

  kIterator_.iterate();
  std::cout << "******************************k*********************************************\n";
  turbflowFieldKE_.getk().show();
  epsilonIterator_.iterate();
  std::cout << "*****************************epsilon**********************************************\n";
  turbflowFieldKE_.geteps().show();
  nuTUpdate();
  std::cout << "*****************************nuT**********************************************\n";
  turbflowFieldKE_.getnuT().show();

  // Compute FGH
  TurbulentFGHIteratorKE_.iterate();
  // Set global boundary values
  wallFGHIterator_.iterate();
  // TODO WS1: compute the right hand side (RHS)
  rhsIterator_.iterate();
  // Solve for pressure
  solver_->solve();
  // TODO WS2: communicate pressure values
  ppmTurbulentKE_.communicatePressure();
  // Compute velocity
  velocityIterator_.iterate();
  obstacleIterator_.iterate();
  // TODO WS2: communicate velocity values
  ppmTurbulentKE_.communicateVelocities();
  // Iterate for velocities on the boundary
  wallVelocityIterator_.iterate();

  // Calculate viscosity
}

void TurbulentSimulationKE::setTimeStep() {
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

void TurbulentSimulationKE::hUpdate() { hIterator_.iterate(); }

void TurbulentSimulationKE::nuTUpdate() { nuTIteratorKE_.iterate(); }

void TurbulentSimulationKE::plotVTK(int timeStep, RealType simulationTime) {
  Stencils::TurbulentVTKStencilKE     TurbulentvtkStencilKE(parameters_);
  FieldIterator<TurbulentFlowFieldKE> TurbulentvtkIteratorKE(
    turbflowFieldKE_, parameters_, TurbulentvtkStencilKE, 1, 0
  );

  TurbulentvtkIteratorKE.iterate();
  TurbulentvtkStencilKE.write(turbflowFieldKE_, timeStep, simulationTime);
}