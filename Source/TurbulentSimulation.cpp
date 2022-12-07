#include "TurbulentSimulation.hpp"

TurbulentSimulation::TurbulentSimulation(Parameters& parameters, TurbulentFlowField& flowField):
  Simulation(parameters, flowField),
  turbflowField_(flowField),
  TurbulentFGHStencil_(parameters),
  TurbulentFGHIterator_(turbflowField_, parameters, TurbulentFGHStencil_),
  nuTStencil_(parameters),
  nuTIterator_(turbflowField_, parameters, nuTStencil_),
  hStencil_(parameters),
  hIterator_(turbflowField_, parameters, hStencil_, 1, 0) {}

void TurbulentSimulation::solveTimestep() {}

void TurbulentSimulation::hUpdate() { hIterator_.iterate(); }