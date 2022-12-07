#include "TurbulentSimulation.hpp"

TurbulentSimulation::TurbulentSimulation(Parameters& parameters, FlowField& flowField):
  Simulation(parameters, flowField),
  TurbulentFGHStencil_(parameters),
  TurbulentFGHIterator_(flowField_, parameters, TurbulentFGHStencil_),
  nuTStencil_(parameters),
  nuTIterator_(flowField_, parameters, nuTStencil_),
  hStencil_(parameters),
  hIterator_(flowField_, parameters, hStencil_) {}

void TurbulentSimulation::solveTimestep() {}

void TurbulentSimulation::hUpdate() { hIterator_.iterate(); }