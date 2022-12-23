#pragma once
#include "FlowField.hpp"

class TurbulentFlowFieldKE: public FlowField {
private:
  ScalarField k_;
  ScalarField eps_;

public:
  TurbulentFlowFieldKE(const Parameters& parameters):
    FlowField(parameters),
    k_(ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3)),
    eps_(ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3)) {}
  ~TurbulentFlowFieldKE() = default;

  ScalarField& getk() { return k_; }
  ScalarField& geteps() { return eps_; }
};
