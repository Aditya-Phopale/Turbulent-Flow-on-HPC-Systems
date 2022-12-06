#include "TurbulentSimulation.hpp"

TurbulentSimulation::TurbulentSimulation(Parameters& parameters, FlowField& flowField):
  Simulation(parameters, flowField),
  TurbulentFGHStencil_(parameters_),
  TurbulentFGHIterator_(flowField_, parameters_, TurbulentFGHStencil_),
  nuTStencil_(parameters_),
  nuTIterator_(flowField_, parameters_, nuTStencil_) {}

void TurbulentSimulation::solveTimestep() {}

void TurbulentSimulation::hUpdate() { hIterator_.iterate(); }