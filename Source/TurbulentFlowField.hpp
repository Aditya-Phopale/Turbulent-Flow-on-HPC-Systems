#pragma once
#include "FlowField.hpp"

class TurbulentFlowField: public FlowField {
private:
  ScalarField h_;
  ScalarField nuT_;

public:
  TurbulentFlowField(const Parameters& parameters):
    FlowField(parameters),
    h_(ScalarField(getNx() + 3, getNy() + 3)),
    nuT_(ScalarField(getNx() + 3, getNy() + 3)) {}
  virtual ~TurbulentFlowField() override = default;

  virtual ScalarField& getheight() override {
    std::cout << "Hello123\n";
    return h_;
  }
  virtual ScalarField& getnuT() override { return nuT_; }

  virtual void geth(RealType& height, int i, int j) override { height = getheight().getScalar(i, j); }
  virtual void getviscosity(RealType& viscosity, int i, int j) override { viscosity = getnuT().getScalar(i, j); }
};