#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"
#include "TurbulentFlowFieldKE.hpp"

namespace Stencils {

  class FGHStencil: public FieldStencil<FlowField> {
  private:
    // A local velocity variable that will be used to approximate derivatives. Size matches 3D
    // case, but can be used for 2D as well.
    RealType localVelocity_[27 * 3];
    RealType localMeshsize_[27 * 3];

  public:
    FGHStencil(const Parameters& parameters);
    ~FGHStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;
  };

  class TurbulentFGHStencil: public FieldStencil<TurbulentFlowField> {
  private:
    // A local velocity variable that will be used to approximate derivatives. Size matches 3D
    // case, but can be used for 2D as well.
    RealType localVelocity_[27 * 3];
    RealType localViscosity_[27 * 3];
    RealType localMeshsize_[27 * 3];

  public:
    TurbulentFGHStencil(const Parameters& parameters);
    ~TurbulentFGHStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

  class TurbulentFGHStencilKE: public FieldStencil<TurbulentFlowFieldKE> {
  private:
    // A local velocity variable that will be used to approximate derivatives. Size matches 3D
    // case, but can be used for 2D as well.
    RealType localVelocity_[27 * 3];
    RealType localViscosity_[27];
    RealType localKineticEnergy_[27];
    RealType localMeshsize_[27 * 3];

  public:
    TurbulentFGHStencilKE(const Parameters& parameters);
    ~TurbulentFGHStencilKE() override = default;

    void apply(TurbulentFlowFieldKE& flowField, int i, int j) override;
  };

} // namespace Stencils
