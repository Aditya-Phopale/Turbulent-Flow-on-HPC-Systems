#pragma once
#include "FlowField.hpp"

class TurbulentFlowField: public FlowField {
private:
  ScalarField h_;
  ScalarField nuT_;

public:
  TurbulentFlowField(const Parameters& parameters):
    FlowField(parameters),
    h_(ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3)),
    nuT_(ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3)) {}
  ~TurbulentFlowField() = default;

  ScalarField& getheight() { return h_; }
  ScalarField& getnuT() { return nuT_; }

  void geth(RealType& height, int i, int j) { height = getheight().getScalar(i, j); }
  void getviscosity(RealType& viscosity, int i, int j) { viscosity = getnuT().getScalar(i, j); }
};