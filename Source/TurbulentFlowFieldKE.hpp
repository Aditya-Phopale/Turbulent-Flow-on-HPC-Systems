#pragma once
#include "FlowField.hpp"

class TurbulentFlowFieldKE: public FlowField {
private:
  const Parameters& parameters_;
  ScalarField       k_old_;
  ScalarField       k_;
  ScalarField       eps_old_;
  ScalarField       eps_;
  ScalarField       h_;
  ScalarField       nuT_;

public:
  TurbulentFlowFieldKE(const Parameters& parameters):
    FlowField(parameters),
    k_old_(ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3)),
    k_(ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3)),
    eps_old_(ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3)),
    eps_(ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3)),
    h_(ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3)),
    nuT_(ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3)),
    parameters_(parameters) {}
  ~TurbulentFlowFieldKE() = default;

  ScalarField& getk() { return k_; }
  ScalarField& geteps() { return eps_; }
  ScalarField& getkold() { return k_old_; }
  ScalarField& getepsold() { return eps_old_; }
  ScalarField& getheight() { return h_; }
  ScalarField& getnuT() { return nuT_; }

  void getviscosity(RealType& viscosity, int i, int j) { viscosity = getnuT().getScalar(i, j); }
  void getviscosity(RealType& viscosity, int i, int j, int k) { viscosity = getnuT().getScalar(i, j, k); }
  void updatekold() {
    for (int i = 0; i < parameters_.parallel.localSize[0] + 3; i++) {
      for (int j = 0; j < parameters_.parallel.localSize[1] + 3; j++) {
        k_old_.getScalar(i, j) = k_.getScalar(i, j);
      }
    }
  }
  void updateepsold() {
    for (int i = 0; i < parameters_.parallel.localSize[0] + 3; i++) {
      for (int j = 0; j < parameters_.parallel.localSize[1] + 3; j++) {
        eps_old_.getScalar(i, j) = eps_.getScalar(i, j);
      }
    }
  }
};