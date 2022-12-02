#include "PetscParallelManager.hpp"

ParallelManagers::PetscParallelManager::PetscParallelManager(Parameters& parameters, FlowField& flowfield):
  parameters_(parameters),
  flowfield_(flowfield) {}

void ParallelManagers::PetscParallelManager::communicatePressure() {
  std::vector<RealType> leftSend{0};
  std::vector<RealType> rightSend{0};
  std::vector<RealType> topSend{0};
  std::vector<RealType> bottomSend{0};
  std::vector<RealType> frontSend{0};
  std::vector<RealType> backSend{0};

  std::vector<RealType> leftReceive{0};
  std::vector<RealType> rightReceive{0};
  std::vector<RealType> topReceive{0};
  std::vector<RealType> bottomReceive{0};
  std::vector<RealType> frontReceive{0};
  std::vector<RealType> backReceive{0};

  if (parameters_.geometry.dim == 2) {
    leftSend.resize(parameters_.parallel.localSize[1] + 3);
    rightSend.resize(parameters_.parallel.localSize[1] + 3);
    topSend.resize(parameters_.parallel.localSize[0] + 3);
    bottomSend.resize(parameters_.parallel.localSize[0] + 3);
    frontSend = {0};
    backSend  = {0};
    leftReceive.resize(parameters_.parallel.localSize[1] + 3);
    rightReceive.resize(parameters_.parallel.localSize[1] + 3);
    topReceive.resize(parameters_.parallel.localSize[0] + 3);
    bottomReceive.resize(parameters_.parallel.localSize[0] + 3);
    frontReceive = {0};
    backReceive  = {0};
  } else if (parameters_.geometry.dim == 3) {
    leftSend.resize((parameters_.parallel.localSize[1] + 3) * (parameters_.parallel.localSize[2] + 3));
    rightSend.resize((parameters_.parallel.localSize[1] + 3) * (parameters_.parallel.localSize[2] + 3));
    topSend.resize((parameters_.parallel.localSize[0] + 3) * (parameters_.parallel.localSize[2] + 3));
    bottomSend.resize((parameters_.parallel.localSize[0] + 3) * (parameters_.parallel.localSize[2] + 3));
    frontSend.resize((parameters_.parallel.localSize[0] + 3) * (parameters_.parallel.localSize[1] + 3));
    backSend.resize((parameters_.parallel.localSize[0] + 3) * (parameters_.parallel.localSize[1] + 3));
    leftReceive.resize((parameters_.parallel.localSize[1] + 3) * (parameters_.parallel.localSize[2] + 3));
    rightReceive.resize((parameters_.parallel.localSize[1] + 3) * (parameters_.parallel.localSize[2] + 3));
    topReceive.resize((parameters_.parallel.localSize[0] + 3) * (parameters_.parallel.localSize[2] + 3));
    bottomReceive.resize((parameters_.parallel.localSize[0] + 3) * (parameters_.parallel.localSize[2] + 3));
    frontReceive.resize((parameters_.parallel.localSize[0] + 3) * (parameters_.parallel.localSize[1] + 3));
    backReceive.resize((parameters_.parallel.localSize[0] + 3) * (parameters_.parallel.localSize[1] + 3));
  }

  Stencils::PressureBufferFillStencil pfill_(
    parameters_, leftSend, rightSend, topSend, bottomSend, frontSend, backSend
  );
  ParallelBoundaryIterator<FlowField> pressureFillIterator(flowfield_, parameters_, pfill_, 0, 0);
  pressureFillIterator.iterate();
  MPI_Sendrecv(
    leftSend.data(),
    leftSend.size(),
    MPI_DOUBLE,
    parameters_.parallel.leftNb,
    0,
    rightReceive.data(),
    rightReceive.size(),
    MPI_DOUBLE,
    parameters_.parallel.rightNb,
    0,
    PETSC_COMM_WORLD,
    MPI_STATUS_IGNORE
  );

  MPI_Sendrecv(
    rightSend.data(),
    rightSend.size(),
    MPI_DOUBLE,
    parameters_.parallel.rightNb,
    0,
    leftReceive.data(),
    leftReceive.size(),
    MPI_DOUBLE,
    parameters_.parallel.leftNb,
    0,
    PETSC_COMM_WORLD,
    MPI_STATUS_IGNORE
  );

  MPI_Sendrecv(
    topSend.data(),
    topSend.size(),
    MPI_DOUBLE,
    parameters_.parallel.topNb,
    0,
    bottomReceive.data(),
    bottomReceive.size(),
    MPI_DOUBLE,
    parameters_.parallel.bottomNb,
    0,
    PETSC_COMM_WORLD,
    MPI_STATUS_IGNORE
  );

  MPI_Sendrecv(
    bottomSend.data(),
    bottomSend.size(),
    MPI_DOUBLE,
    parameters_.parallel.bottomNb,
    0,
    topReceive.data(),
    topReceive.size(),
    MPI_DOUBLE,
    parameters_.parallel.topNb,
    0,
    PETSC_COMM_WORLD,
    MPI_STATUS_IGNORE
  );

  MPI_Sendrecv(
    frontSend.data(),
    frontSend.size(),
    MPI_DOUBLE,
    parameters_.parallel.frontNb,
    0,
    backReceive.data(),
    backReceive.size(),
    MPI_DOUBLE,
    parameters_.parallel.backNb,
    0,
    PETSC_COMM_WORLD,
    MPI_STATUS_IGNORE
  );

  MPI_Sendrecv(
    backSend.data(),
    backSend.size(),
    MPI_DOUBLE,
    parameters_.parallel.backNb,
    0,
    frontReceive.data(),
    frontReceive.size(),
    MPI_DOUBLE,
    parameters_.parallel.frontNb,
    0,
    PETSC_COMM_WORLD,
    MPI_STATUS_IGNORE
  );

  Stencils::PressureBufferReadStencil pread_(
    parameters_, leftReceive, rightReceive, topReceive, bottomReceive, frontReceive, backReceive
  );
  ParallelBoundaryIterator<FlowField> pressureReadIterator(flowfield_, parameters_, pread_, 0, 0);
  pressureReadIterator.iterate();
}

void ParallelManagers::PetscParallelManager::communicateVelocities() {
  std::vector<RealType> leftSend;
  std::vector<RealType> rightSend;
  std::vector<RealType> topSend;
  std::vector<RealType> bottomSend;
  std::vector<RealType> frontSend;
  std::vector<RealType> backSend;

  std::vector<RealType> leftReceive;
  std::vector<RealType> rightReceive;
  std::vector<RealType> topReceive;
  std::vector<RealType> bottomReceive;
  std::vector<RealType> frontReceive;
  std::vector<RealType> backReceive;

  if (parameters_.geometry.dim == 2) {
    leftSend.resize(2 *(parameters_.parallel.localSize[1] + 3));
    rightSend.resize(2 *(parameters_.parallel.localSize[1] + 3));
    topSend.resize(2 *(parameters_.parallel.localSize[0] + 3));
    bottomSend.resize(2 *(parameters_.parallel.localSize[0] + 3));
    frontSend = {0};
    backSend  = {0};
    leftReceive.resize(2 *(parameters_.parallel.localSize[1] + 3));
    rightReceive.resize(2 *(parameters_.parallel.localSize[1] + 3));
    topReceive.resize(2 *(parameters_.parallel.localSize[0] + 3));
    bottomReceive.resize(2 *(parameters_.parallel.localSize[0] + 3));
    frontReceive = {0};
    backReceive  = {0};
  } else if (parameters_.geometry.dim == 3) {
    leftSend.resize(3 * ((parameters_.parallel.localSize[1] + 3) * (parameters_.parallel.localSize[2] + 3)));
    rightSend.resize(3 * ((parameters_.parallel.localSize[1] + 3) * (parameters_.parallel.localSize[2] + 3)));
    topSend.resize(3 * ((parameters_.parallel.localSize[0] + 3) * (parameters_.parallel.localSize[2] + 3)));
    bottomSend.resize(3 * ((parameters_.parallel.localSize[0] + 3) * (parameters_.parallel.localSize[2] + 3)));
    frontSend.resize(3 * ((parameters_.parallel.localSize[0] + 3) * (parameters_.parallel.localSize[1] + 3)));
    backSend.resize(3 * ((parameters_.parallel.localSize[0] + 3) * (parameters_.parallel.localSize[1] + 3)));
    leftReceive.resize(3 * ((parameters_.parallel.localSize[1] + 3) * (parameters_.parallel.localSize[2] + 3)));
    rightReceive.resize(3 * ((parameters_.parallel.localSize[1] + 3) * (parameters_.parallel.localSize[2] + 3)));
    topReceive.resize(3 * ((parameters_.parallel.localSize[0] + 3) * (parameters_.parallel.localSize[2] + 3)));
    bottomReceive.resize(3 * (parameters_.parallel.localSize[0] + 3) * (parameters_.parallel.localSize[2] + 3));
    frontReceive.resize(3 * (parameters_.parallel.localSize[0] + 3) * (parameters_.parallel.localSize[1] + 3));
    backReceive.resize(3 * (parameters_.parallel.localSize[0] + 3) * (parameters_.parallel.localSize[1] + 3));
  }

  Stencils::VelocityBufferFillStencil vfill_(
    parameters_, leftSend, rightSend, topSend, bottomSend, frontSend, backSend
  );
  ParallelBoundaryIterator<FlowField> velocityFillIterator(flowfield_, parameters_, vfill_, 0, 0);
  velocityFillIterator.iterate();
  MPI_Sendrecv(
    leftSend.data(),
    leftSend.size(),
    MPI_DOUBLE,
    parameters_.parallel.leftNb,
    0,
    rightReceive.data(),
    rightReceive.size(),
    MPI_DOUBLE,
    parameters_.parallel.rightNb,
    0,
    PETSC_COMM_WORLD,
    MPI_STATUS_IGNORE
  );

  MPI_Sendrecv(
    rightSend.data(),
    rightSend.size(),
    MPI_DOUBLE,
    parameters_.parallel.rightNb,
    0,
    leftReceive.data(),
    leftReceive.size(),
    MPI_DOUBLE,
    parameters_.parallel.leftNb,
    0,
    PETSC_COMM_WORLD,
    MPI_STATUS_IGNORE
  );

  MPI_Sendrecv(
    topSend.data(),
    topSend.size(),
    MPI_DOUBLE,
    parameters_.parallel.topNb,
    0,
    bottomReceive.data(),
    bottomReceive.size(),
    MPI_DOUBLE,
    parameters_.parallel.bottomNb,
    0,
    PETSC_COMM_WORLD,
    MPI_STATUS_IGNORE
  );

  MPI_Sendrecv(
    bottomSend.data(),
    bottomSend.size(),
    MPI_DOUBLE,
    parameters_.parallel.bottomNb,
    0,
    topReceive.data(),
    topReceive.size(),
    MPI_DOUBLE,
    parameters_.parallel.topNb,
    0,
    PETSC_COMM_WORLD,
    MPI_STATUS_IGNORE
  );

  MPI_Sendrecv(
    frontSend.data(),
    frontSend.size(),
    MPI_DOUBLE,
    parameters_.parallel.frontNb,
    0,
    backReceive.data(),
    backReceive.size(),
    MPI_DOUBLE,
    parameters_.parallel.backNb,
    0,
    PETSC_COMM_WORLD,
    MPI_STATUS_IGNORE
  );

  MPI_Sendrecv(
    backSend.data(),
    backSend.size(),
    MPI_DOUBLE,
    parameters_.parallel.backNb,
    0,
    frontReceive.data(),
    frontReceive.size(),
    MPI_DOUBLE,
    parameters_.parallel.frontNb,
    0,
    PETSC_COMM_WORLD,
    MPI_STATUS_IGNORE
  );

  Stencils::VelocityBufferReadStencil vread_(
    parameters_, leftReceive, rightReceive, topReceive, bottomReceive, frontReceive, backReceive
  );
  ParallelBoundaryIterator<FlowField> velocityReadIterator(flowfield_, parameters_, vread_, 0, 0);
  velocityReadIterator.iterate();
}
