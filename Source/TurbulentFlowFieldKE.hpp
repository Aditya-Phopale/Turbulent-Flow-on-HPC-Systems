#pragma once
#include "FlowField.hpp"

class TurbulentFlowFieldKE: public FlowField {
private:
  ScalarField h_;
  ScalarField nuT_;
  ScalarField ke_;
  ScalarField e_;

public:
  TurbulentFlowFieldKE(const Parameters& parameters):
    FlowField(parameters),
    h_(ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3)),
    nuT_(ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3)),
    ke_(ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3)),
    e_(ScalarField(parameters.parallel.localSize[0] + 3, parameters.parallel.localSize[1] + 3)) {}
  ~TurbulentFlowFieldKE() = default;

  ScalarField& getheight() { return h_; }
  ScalarField& getnuT() { return nuT_; }
  ScalarField& getke() { return ke_; }
  ScalarField& gete() { return e_; }

  void geth(RealType& height, int i, int j) { height = getheight().getScalar(i, j); }
  void getviscosity(RealType& viscosity, int i, int j) { viscosity = getnuT().getScalar(i, j); }
  void getviscosity(RealType& viscosity, int i, int j, int k) { viscosity = getnuT().getScalar(i, j, k); }
  void getenergy(RealType& ke, int i, int j) { ke = getke().getScalar(i, j); }
  void getenergy(RealType& ke, int i, int j, int k) { ke = getke().getScalar(i, j, k); }
  void getdiss(RealType& e, int i, int j) { e = gete().getScalar(i, j); }
  void getdiss(RealType& e, int i, int j, int k) { e = gete().getScalar(i, j, k); }
};
