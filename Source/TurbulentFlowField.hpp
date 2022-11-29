#include "FlowField.hpp"

class TurbulentFlowField: public FlowField {
private:
  ScalarField h_;
  ScalarField nuT_;

  public:
  TurbulentFlowField(const Parameters& parameters) : FlowField(parameters), h_(ScalarField(getNx() + 3, getNy() + 3)), nuT_(ScalarField(getNx() + 3, getNy() + 3)){

  }
  virtual ~TurbulentFlowField() = default;

  ScalarField& geth(){
    return h_;
  }
  ScalarField& getnuT(){
    return nuT_;
  }
};