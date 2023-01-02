#pragma once

#include "Definitions.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"
#include "TurbulentFlowFieldKE.hpp"

namespace Stencils {

  // Load the local velocity cube with relevant velocities of the 2D plane
  inline void loadLocalVelocity2D(FlowField& flowField, RealType* const localVelocity, int i, int j) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        const RealType* const point                  = flowField.getVelocity().getVector(i + column, j + row);
        localVelocity[39 + 9 * row + 3 * column]     = point[0]; // x-component
        localVelocity[39 + 9 * row + 3 * column + 1] = point[1]; // y-component
      }
    }
  }

  inline void loadLocalViscosity2D(TurbulentFlowField& flowField, RealType* const localViscosity, int i, int j) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        localViscosity[39 + 9 * row + 3 * column] = flowField.getnuT().getScalar(i + column, j + row); // x-component
      }
    }
  }

  inline void loadLocalViscosity2D(TurbulentFlowFieldKE& flowField, RealType* const localViscosity, int i, int j) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        localViscosity[39 + 9 * row + 3 * column] = flowField.getnuT().getScalar(i + column, j + row); // x-component
      }
    }
  }

  inline void loadLocalK2D(TurbulentFlowFieldKE& flowField, RealType* const localK, int i, int j) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        localK[39 + 9 * row + 3 * column] = flowField.getk().getScalar(i + column, j + row); // x-component
      }
    }
  }

  inline void loadLocalEpsilon2D(TurbulentFlowFieldKE& flowField, RealType* const localEpsilon, int i, int j) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        localEpsilon[39 + 9 * row + 3 * column] = flowField.geteps().getScalar(i + column, j + row); // x-component
      }
    }
  }

  // Load the local velocity cube with surrounding velocities
  inline void loadLocalVelocity3D(FlowField& flowField, RealType* const localVelocity, int i, int j, int k) {
    for (int layer = -1; layer <= 1; layer++) {
      for (int row = -1; row <= 1; row++) {
        for (int column = -1; column <= 1; column++) {
          const RealType* const point = flowField.getVelocity().getVector(i + column, j + row, k + layer);
          localVelocity[39 + 27 * layer + 9 * row + 3 * column]     = point[0]; // x-component
          localVelocity[39 + 27 * layer + 9 * row + 3 * column + 1] = point[1]; // y-component
          localVelocity[39 + 27 * layer + 9 * row + 3 * column + 2] = point[2]; // z-component
        }
      }
    }
  }

  inline void loadLocalViscosity3D(TurbulentFlowField& flowField, RealType* const localViscosity, int i, int j, int k) {
    for (int layer = -1; layer <= 1; layer++) {
      for (int row = -1; row <= 1; row++) {
        for (int column = -1; column <= 1; column++) {
          localViscosity[39 + 27 * layer + 9 * row + 3 * column] = flowField.getnuT().getScalar(
            i + column, j + row, k + layer
          ); // x-component
        }
      }
    }
  }

  // Load local meshsize for 2D -> same as loadLocalVelocity2D, but invoking call to meshsize-ptr
  inline void loadLocalMeshsize2D(const Parameters& parameters, RealType* const localMeshsize, int i, int j) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        localMeshsize[39 + 9 * row + 3 * column]     = parameters.meshsize->getDx(i + column, j + row);
        localMeshsize[39 + 9 * row + 3 * column + 1] = parameters.meshsize->getDy(i + column, j + row);
      }
    }
  }

  // Load local meshsize for 3D
  inline void loadLocalMeshsize3D(const Parameters& parameters, RealType* const localMeshsize, int i, int j, int k) {
    for (int layer = -1; layer <= 1; layer++) {
      for (int row = -1; row <= 1; row++) {
        for (int column = -1; column <= 1; column++) {
          localMeshsize[39 + 27 * layer + 9 * row + 3 * column] = parameters.meshsize->getDx(
            i + column, j + row, k + layer
          );
          localMeshsize[39 + 27 * layer + 9 * row + 3 * column + 1] = parameters.meshsize->getDy(
            i + column, j + row, k + layer
          );
          localMeshsize[39 + 27 * layer + 9 * row + 3 * column + 2] = parameters.meshsize->getDz(
            i + column, j + row, k + layer
          );
        }
      }
    }
  }

  // Maps an index and a component to the corresponding value in the cube.
  inline int mapd(int i, int j, int k, int component) { return 39 + 27 * k + 9 * j + 3 * i + component; }

  // Derivative functions. They are applied to a cube of 3x3x3 cells. lv stands for the local velocity, lm represents
  // the local mesh sizes dudx <-> first derivative of u-component of velocity field w.r.t. x-direction.
  inline RealType dudx(const RealType* const lv, const RealType* const lm) {
    // Evaluate dudx in the cell center by a central difference
    const int index0 = mapd(0, 0, 0, 0);
    const int index1 = mapd(-1, 0, 0, 0);
    return (lv[index0] - lv[index1]) / lm[index0];
  }

  inline RealType dudy(const RealType* const lv, const RealType* const lm) {
    // Evaluate dudx in the cell center by a central difference
    const int indexc = mapd(0, 0, 0, 0);
    const int indexl = mapd(-1, 0, 0, 0);
    RealType  cv     = 0.5 * (lv[indexc] + lv[indexl]);

    const int indexlt = mapd(-1, 1, 0, 0);
    const int indext  = mapd(0, 1, 0, 0);
    RealType  tv      = 0.5 * (lv[indexlt] + lv[indext]);

    const int indexbl = mapd(-1, -1, 0, 0);
    const int indexb  = mapd(0, -1, 0, 0);
    RealType  bv      = 0.5 * (lv[indexbl] + lv[indexb]);

    return 0.5
           * ((tv - cv) / (0.5 * (lm[mapd(0, 1, 0, 1)] + lm[mapd(0, 0, 0, 1)])) + (cv - bv) / (0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, -1, 0, 1)])));
  }

  inline RealType dvdx(const RealType* const lv, const RealType* const lm) {
    // Evaluate dudx in the cell center by a central difference
    const int indexc = mapd(0, 0, 0, 1);
    const int indexb = mapd(0, -1, 0, 1);
    RealType  cv     = 0.5 * (lv[indexc] + lv[indexb]);

    const int indexl  = mapd(-1, 0, 0, 1);
    const int indexbl = mapd(-1, -1, 0, 1);
    RealType  lev     = 0.5 * (lv[indexl] + lv[indexbl]);

    const int indexr  = mapd(1, 0, 0, 1);
    const int indexbr = mapd(1, -1, 0, 1);
    RealType  rv      = 0.5 * (lv[indexr] + lv[indexbr]);

    return 0.5
           * ((rv - cv) / (0.5 * (lm[mapd(1, 0, 0, 0)] + lm[mapd(0, 0, 0, 0)])) + (cv - lev) / (0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(-1, 0, 0, 0)])));
  }

  inline RealType dudz(const RealType* const lv, const RealType* const lm) {
    // Evaluate dudx in the cell center by a central difference
    const int indexc   = mapd(0, 0, 0, 0);
    const int indexba  = mapd(0, 0, -1, 0);
    const int indexf   = mapd(0, 0, 1, 0);
    const int indexbal = mapd(-1, 0, -1, 0);
    const int indexfl  = mapd(-1, 0, 1, 0);
    const int indexl   = mapd(-1, 0, 0, 0);

    RealType topR     = (lv[indexf] - lv[indexc]) / (0.5 * (lm[mapd(0, 0, 1, 2)] + lm[mapd(0, 0, 0, 2)]));
    RealType topL     = (lv[indexfl] - lv[indexl]) / (0.5 * (lm[mapd(-1, 0, 1, 2)] + lm[mapd(-1, 0, 0, 2)]));
    RealType bottonmL = (lv[indexl] - lv[indexbal]) / (0.5 * (lm[mapd(-1, 0, 0, 2)] + lm[mapd(-1, 0, -1, 2)]));
    RealType bottomR  = (lv[indexc] - lv[indexba]) / (0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, -1, 2)]));
    return 0.25 * (topR + topL + bottomR + bottonmL);
  }
  // Implemented derivaive functions till now are correct

  inline RealType dwdx(const RealType* const lv, const RealType* const lm) {
    // Evaluate dudx in the cell center by a central difference
    const int indexc   = mapd(0, 0, 0, 2);
    const int indexl   = mapd(-1, 0, 0, 2);
    const int indexr   = mapd(1, 0, 0, 2);
    const int indexbal = mapd(-1, 0, -1, 2);
    const int indexba  = mapd(0, 0, -1, 2);
    const int indexbar = mapd(1, 0, -1, 2);

    RealType frontR = (lv[indexr] - lv[indexc]) / (0.5 * (lm[mapd(1, 0, 0, 0)] + lm[mapd(0, 0, 0, 0)]));
    RealType frontL = (lv[indexc] - lv[indexl]) / (0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(-1, 0, 0, 0)]));
    RealType backL  = (lv[indexba] - lv[indexbal]) / (0.5 * (lm[mapd(0, 0, -1, 0)] + lm[mapd(-1, 0, -1, 0)]));
    RealType backR  = (lv[indexbar] - lv[indexba]) / (0.5 * (lm[mapd(1, 0, -1, 0)] + lm[mapd(0, 0, -1, 0)]));
    return 0.25 * (frontR + frontL + backR + backL);
  }

  inline RealType dvdz(const RealType* const lv, const RealType* const lm) {
    // Evaluate dudx in the cell center by a central difference
    const int indexc   = mapd(0, 0, 0, 1);
    const int indexf   = mapd(0, 0, 1, 1);
    const int indexba  = mapd(0, 0, -1, 1);
    const int indexb   = mapd(0, -1, 0, 1);
    const int indexbf  = mapd(0, -1, 1, 1);
    const int indexbba = mapd(0, -1, -1, 1);

    RealType topF     = (lv[indexf] - lv[indexc]) / (0.5 * (lm[mapd(0, 0, 1, 2)] + lm[mapd(0, 0, 0, 2)]));
    RealType topBa    = (lv[indexc] - lv[indexba]) / (0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, -1, 2)]));
    RealType bottomF  = (lv[indexbf] - lv[indexb]) / (0.5 * (lm[mapd(0, -1, 1, 2)] + lm[mapd(0, -1, 0, 2)]));
    RealType bottomBa = (lv[indexb] - lv[indexbba]) / (0.5 * (lm[mapd(0, -1, 0, 2)] + lm[mapd(0, -1, -1, 2)]));
    return 0.25 * (topF + topBa + bottomF + bottomBa);
  }

  inline RealType dwdy(const RealType* const lv, const RealType* const lm) {
    // Evaluate dudx in the cell center by a central difference
    const int indexc   = mapd(0, 0, 0, 2);
    const int indext   = mapd(0, 1, 0, 2);
    const int indexb   = mapd(0, -1, 0, 2);
    const int indexba  = mapd(0, 0, -1, 2);
    const int indexbat = mapd(0, 1, -1, 2);
    const int indexbab = mapd(0, -1, -1, 2);

    RealType topF     = (lv[indext] - lv[indexc]) / (0.5 * (lm[mapd(0, 1, 0, 1)] + lm[mapd(0, 0, 0, 1)]));
    RealType topBa    = (lv[indexbat] - lv[indexba]) / (0.5 * (lm[mapd(0, 1, -1, 1)] + lm[mapd(0, 0, -1, 1)]));
    RealType BottomF  = (lv[indexc] - lv[indexb]) / (0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 0, -1, 1)]));
    RealType BottomBa = (lv[indexba] - lv[indexbab]) / (0.5 * (lm[mapd(0, 0, -1, 1)] + lm[mapd(0, -1, -1, 1)]));
    return 0.25 * (topF + topBa + BottomF + BottomBa);
  }

  inline RealType dvdy(const RealType* const lv, const RealType* const lm) {
    const int index0 = mapd(0, 0, 0, 1);
    const int index1 = mapd(0, -1, 0, 1);
    return (lv[index0] - lv[index1]) / lm[index0];
  }

  inline RealType dwdz(const RealType* const lv, const RealType* const lm) {
    const int index0 = mapd(0, 0, 0, 2);
    const int index1 = mapd(0, 0, -1, 2);
    return (lv[index0] - lv[index1]) / lm[index0];
  }

  inline RealType dkdx(const RealType* const lk, const RealType* const lm) {
    // given in pg 164 (followed same as pressure)
    const int index0 = mapd(1, 0, 0, 0);
    const int index1 = mapd(0, 0, 0, 0);
    return (lk[index0] - lk[index1]) / lm[index0];
  }

  inline RealType dkdy(const RealType* const lk, const RealType* const lm) {
    // given in pg 164 (followed same as pressure)
    const int index0 = mapd(0, 1, 0, 0);
    const int index1 = mapd(0, 0, 0, 0);
    return (lk[index0] - lk[index1]) / lm[index0];
  }

  inline RealType d2udx2(const RealType* const lv, const RealType* const lm) {
    // Evaluate the second derivative at the location of the u-component of the velocity field;
    // we therefore use the two neighbouring u-components and assume arbitrary mesh sizes in both
    // directions -> the formula arises from a straight-forward taylor expansion
    // -> for equal meshsizes, we obtain the usual [1 -2 1]-like stencil.
    const int indexM1 = mapd(-1, 0, 0, 0);
    const int index0  = mapd(0, 0, 0, 0);
    const int indexP1 = mapd(1, 0, 0, 0);

    const RealType dx0   = lm[index0];
    const RealType dx1   = lm[indexP1];
    const RealType dxSum = dx0 + dx1;
    return 2.0 * (lv[indexP1] / (dx1 * dxSum) - lv[index0] / (dx1 * dx0) + lv[indexM1] / (dx0 * dxSum));
  }

  // Second derivative of turbulent u-component w.r.t. x-direction, evaluated at the location of the u-component.
  inline RealType d2udx2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    // Evaluate the second derivative at the location of the u-component of the velocity field;
    // we therefore use the two neighbouring u-components and assume arbitrary mesh sizes in both
    // directions -> the formula arises from a straight-forward taylor expansion
    // -> for equal meshsizes, we obtain the usual [1 -2 1]-like stencil.
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    RealType dudxR = (lv[mapd(1, 0, 0, 0)] - lv[mapd(0, 0, 0, 0)]) / dx_P1;
    RealType dudxL = (lv[mapd(0, 0, 0, 0)] - lv[mapd(-1, 0, 0, 0)]) / dx_0;

    RealType ViscR = lvis[mapd(1, 0, 0, 0)] + 1 / parameters.flow.Re;
    RealType ViscL = lvis[mapd(0, 0, 0, 0)] + 1 / parameters.flow.Re;

    return (ViscR * dudxR - ViscL * dudxL) / dx_0;
  }

  inline RealType d2udy2(const RealType* const lv, const RealType* const lm) {
    // Average mesh sizes, since the component u is located in the middle of the cell's face.
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);
    const RealType dySum = dy0 + dy1;
    return 2.0
           * (lv[mapd(0, 1, 0, 0)] / (dy1 * dySum) - lv[mapd(0, 0, 0, 0)] / (dy1 * dy0) + lv[mapd(0, -1, 0, 0)] / (dy0 * dySum));
  }

  inline RealType d2udy2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    RealType dudyT = (lv[mapd(0, 1, 0, 0)] - lv[mapd(0, 0, 0, 0)]) / dy1;
    RealType dudyB = (lv[mapd(0, 0, 0, 0)] - lv[mapd(0, -1, 0, 0)]) / dy0;

    RealType ViscT1 = lvis[mapd(0, 1, 0, 0)];
    RealType ViscT2 = lvis[mapd(1, 1, 0, 0)];
    RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    RealType ViscM2 = lvis[mapd(1, 0, 0, 0)];
    RealType ViscB1 = lvis[mapd(0, -1, 0, 0)];
    RealType ViscB2 = lvis[mapd(1, -1, 0, 0)];

    RealType ViscT = (0.5 * (0.5 * (ViscT1 * dx_P1 + ViscT2 * dx_0) / dx1) * dy_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dy_P1
                     ) / dy1
                     + 1 / parameters.flow.Re;

    RealType ViscB = (0.5 * (0.5 * (ViscB1 * dx_P1 + ViscB2 * dx_0) / dx1) * dy_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dy_M1
                     ) / dy0
                     + 1 / parameters.flow.Re;

    return (ViscT * dudyT - ViscB * dudyB) / dy_0;
  }

  inline RealType d2udz2(const RealType* const lv, const RealType* const lm) {
    const RealType dz_M1 = lm[mapd(0, 0, -1, 2)];
    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(0, 0, 1, 2)];
    const RealType dz0   = 0.5 * (dz_0 + dz_M1);
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);
    const RealType dzSum = dz0 + dz1;
    return 2.0
           * (lv[mapd(0, 0, 1, 0)] / (dz1 * dzSum) - lv[mapd(0, 0, 0, 0)] / (dz1 * dz0) + lv[mapd(0, 0, -1, 0)] / (dz0 * dzSum));
  }

  inline RealType d2udz2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dz_M1 = lm[mapd(0, 0, -1, 2)];
    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(0, 0, 1, 2)];
    const RealType dz0   = 0.5 * (dz_0 + dz_M1);
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    RealType dudzT = (lv[mapd(0, 0, 1, 0)] - lv[mapd(0, 0, 0, 0)]) / dz1;
    RealType dudzB = (lv[mapd(0, 0, 0, 0)] - lv[mapd(0, 0, -1, 0)]) / dz0;

    RealType ViscT1 = lvis[mapd(0, 0, 1, 0)];
    RealType ViscT2 = lvis[mapd(1, 0, 1, 0)];
    RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    RealType ViscM2 = lvis[mapd(1, 0, 0, 0)];
    RealType ViscB1 = lvis[mapd(0, 0, -1, 0)];
    RealType ViscB2 = lvis[mapd(1, 0, -1, 0)];

    RealType ViscT = (0.5 * (0.5 * (ViscT1 * dx_P1 + ViscT2 * dx_0) / dx1) * dz_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dz_P1
                     ) / dz1
                     + 1 / parameters.flow.Re;

    RealType ViscB = (0.5 * (0.5 * (ViscB1 * dx_P1 + ViscB2 * dx_0) / dx1) * dz_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dz_M1
                     ) / dz0
                     + 1 / parameters.flow.Re;

    return (ViscT * dudzT - ViscB * dudzB) / dz_0;
  }

  // Second derivative of the v-component, evaluated at the location of the v-component.
  inline RealType d2vdx2(const RealType* const lv, const RealType* const lm) {
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);
    const RealType dxSum = dx0 + dx1;
    return 2.0
           * (lv[mapd(1, 0, 0, 1)] / (dx1 * dxSum) - lv[mapd(0, 0, 0, 1)] / (dx1 * dx0) + lv[mapd(-1, 0, 0, 1)] / (dx0 * dxSum));
  }

  inline RealType d2vdx2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType ViscR1 = lvis[mapd(1, 0, 0, 0)];
    const RealType ViscR2 = lvis[mapd(1, 1, 0, 0)];
    const RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    const RealType ViscM2 = lvis[mapd(0, 1, 0, 0)];
    const RealType ViscL1 = lvis[mapd(-1, 0, 0, 0)];
    const RealType ViscL2 = lvis[mapd(-1, 1, 0, 0)];

    RealType dvdxR = (lv[mapd(1, 0, 0, 1)] - lv[mapd(0, 0, 0, 1)]) / dx1;
    RealType dvdxL = (lv[mapd(0, 0, 0, 1)] - lv[mapd(-1, 0, 0, 1)]) / dx0;

    RealType ViscR = (0.5 * (0.5 * (ViscR1 * dx_0 + ViscM1 * dx_P1) / (dx1)) * (dy_P1)
                      + 0.5 * (0.5 * (ViscR2 * dx_0 + ViscM2 * dx_P1) / (dx1)) * (dy_0)
                     ) / dy1
                     + 1 / parameters.flow.Re;

    RealType ViscL = (0.5 * (0.5 * (ViscL1 * dx_0 + ViscM1 * dx_M1) / (dx0)) * (dy_P1)
                      + 0.5 * (0.5 * (ViscL2 * dx_0 + ViscM2 * dx_M1) / (dx0)) * (dy_0)
                     ) / dy1
                     + 1 / parameters.flow.Re;

    return (ViscR * dvdxR - ViscL * dvdxL) / dx_0;
  }

  inline RealType d2vdydx(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    RealType dvdxT = lv[mapd(1, 0, 0, 1)] - lv[mapd(0, 0, 0, 1)] / dx1;
    RealType dvdxB = lv[mapd(1, -1, 0, 1)] - lv[mapd(0, -1, 0, 1)] / dx1;

    RealType ViscT1 = lvis[mapd(0, 1, 0, 0)];
    RealType ViscT2 = lvis[mapd(1, 1, 0, 0)];
    RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    RealType ViscM2 = lvis[mapd(1, 0, 0, 0)];
    RealType ViscB1 = lvis[mapd(0, -1, 0, 0)];
    RealType ViscB2 = lvis[mapd(1, -1, 0, 0)];

    RealType ViscT = (0.5 * (0.5 * (ViscT1 * dx_P1 + ViscT2 * dx_0) / dx1) * dy_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dy_P1
                     ) / dy1
                     + 1 / parameters.flow.Re;

    RealType ViscB = (0.5 * (0.5 * (ViscB1 * dx_P1 + ViscB2 * dx_0) / dx1) * dy_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dy_M1
                     ) / dy0
                     + 1 / parameters.flow.Re;

    return (ViscT * dvdxT - ViscB * dvdxB) / dy_0;
  }

  inline RealType d2udxdy(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType ViscR1 = lvis[mapd(1, 0, 0, 0)];
    const RealType ViscR2 = lvis[mapd(1, 1, 0, 0)];
    const RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    const RealType ViscM2 = lvis[mapd(0, 1, 0, 0)];
    const RealType ViscL1 = lvis[mapd(-1, 0, 0, 0)];
    const RealType ViscL2 = lvis[mapd(-1, 1, 0, 0)];

    RealType dudyR = (lv[mapd(0, 1, 0, 0)] - lv[mapd(0, 0, 0, 0)]) / dy1;
    RealType dudyL = (lv[mapd(-1, 1, 0, 0)] - lv[mapd(-1, 0, 0, 0)]) / dy1;

    RealType ViscR = (0.5 * (0.5 * (ViscR1 * dx_0 + ViscM1 * dx_P1) / (dx1)) * (dy_P1)
                      + 0.5 * (0.5 * (ViscR2 * dx_0 + ViscM2 * dx_P1) / (dx1)) * (dy_0)
                     ) / dy1
                     + 1 / parameters.flow.Re;

    RealType ViscL = (0.5 * (0.5 * (ViscL1 * dx_0 + ViscM1 * dx_M1) / (dx0)) * (dy_P1)
                      + 0.5 * (0.5 * (ViscL2 * dx_0 + ViscM2 * dx_M1) / (dx0)) * (dy_0)
                     ) / dy1
                     + 1 / parameters.flow.Re;

    return ((ViscR * dudyR) - (ViscL * dudyL)) / dx_0;
  }

  inline RealType d2wdxdz(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(0, 0, 1, 2)];
    const RealType dz_M1 = lm[mapd(0, 0, -1, 2)];
    const RealType dz0   = 0.5 * (dz_0 + dz_M1);
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);

    RealType dwdxT = lv[mapd(1, 0, 0, 2)] - lv[mapd(0, 0, 0, 2)] / dx1;
    RealType dwdxB = lv[mapd(1, 0, -1, 2)] - lv[mapd(0, 0, -1, 2)] / dx1;

    RealType ViscT1 = lvis[mapd(0, 0, 1, 0)];
    RealType ViscT2 = lvis[mapd(1, 0, 1, 0)];
    RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    RealType ViscM2 = lvis[mapd(1, 0, 0, 0)];
    RealType ViscB1 = lvis[mapd(0, 0, -1, 0)];
    RealType ViscB2 = lvis[mapd(1, 0, -1, 0)];

    RealType ViscT = (0.5 * (0.5 * (ViscT1 * dx_P1 + ViscT2 * dx_0) / dx1) * dz_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dz_P1
                     ) / dz1
                     + 1 / parameters.flow.Re;

    RealType ViscB = (0.5 * (0.5 * (ViscB1 * dx_P1 + ViscB2 * dx_0) / dx1) * dz_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dz_M1
                     ) / dz0
                     + 1 / parameters.flow.Re;

    return (ViscT * dwdxT - ViscB * dwdxB) / dz_0;
  }

  inline RealType d2wdydz(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(1, 0, 0, 1)];
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(0, 0, 1, 2)];
    const RealType dz_M1 = lm[mapd(0, 0, -1, 2)];
    const RealType dz0   = 0.5 * (dz_0 + dz_M1);
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);

    RealType dwdyT = lv[mapd(0, 1, 0, 2)] - lv[mapd(0, 0, 0, 2)] / dy1;
    RealType dwdyB = lv[mapd(0, 1, -1, 2)] - lv[mapd(0, 0, -1, 2)] / dy1;

    RealType ViscT1 = lvis[mapd(0, 0, 1, 0)];
    RealType ViscT2 = lvis[mapd(1, 0, 1, 0)];
    RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    RealType ViscM2 = lvis[mapd(1, 0, 0, 0)];
    RealType ViscB1 = lvis[mapd(0, 0, -1, 0)];
    RealType ViscB2 = lvis[mapd(1, 0, -1, 0)];

    RealType ViscT = (0.5 * (0.5 * (ViscT1 * dy_P1 + ViscT2 * dy_0) / dy1) * dz_0
                      + 0.5 * (0.5 * (ViscM1 * dy_P1 + ViscM2 * dy_0) / dy1) * dz_P1
                     ) / dz1
                     + 1 / parameters.flow.Re;

    RealType ViscB = (0.5 * (0.5 * (ViscB1 * dy_P1 + ViscB2 * dy_0) / dy1) * dz_0
                      + 0.5 * (0.5 * (ViscM1 * dy_P1 + ViscM2 * dy_0) / dy1) * dz_M1
                     ) / dz0
                     + 1 / parameters.flow.Re;

    return (ViscT * dwdyT - ViscB * dwdyB) / dz_0;
  }

  inline RealType d2udzdx(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    RealType dvdxT = lv[mapd(1, 0, 0, 1)] - lv[mapd(0, 0, 0, 1)] / dx1;
    RealType dvdxB = lv[mapd(1, -1, 0, 1)] - lv[mapd(0, -1, 0, 1)] / dx1;

    RealType ViscT1 = lvis[mapd(0, 1, 0, 0)];
    RealType ViscT2 = lvis[mapd(1, 1, 0, 0)];
    RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    RealType ViscM2 = lvis[mapd(1, 0, 0, 0)];
    RealType ViscB1 = lvis[mapd(0, -1, 0, 0)];
    RealType ViscB2 = lvis[mapd(1, -1, 0, 0)];

    RealType ViscT = (0.5 * (0.5 * (ViscT1 * dx_P1 + ViscT2 * dx_0) / dx1) * dy_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dy_P1
                     ) / dy1
                     + 1 / parameters.flow.Re;

    RealType ViscB = (0.5 * (0.5 * (ViscB1 * dx_P1 + ViscB2 * dx_0) / dx1) * dy_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dy_M1
                     ) / dy0
                     + 1 / parameters.flow.Re;

    return (ViscT * dvdxT - ViscB * dvdxB) / dy_0;
  }

  inline RealType d2vdzdy(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    RealType dvdxT = lv[mapd(1, 0, 0, 1)] - lv[mapd(0, 0, 0, 1)] / dx1;
    RealType dvdxB = lv[mapd(1, -1, 0, 1)] - lv[mapd(0, -1, 0, 1)] / dx1;

    RealType ViscT1 = lvis[mapd(0, 1, 0, 0)];
    RealType ViscT2 = lvis[mapd(1, 1, 0, 0)];
    RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    RealType ViscM2 = lvis[mapd(1, 0, 0, 0)];
    RealType ViscB1 = lvis[mapd(0, -1, 0, 0)];
    RealType ViscB2 = lvis[mapd(1, -1, 0, 0)];

    RealType ViscT = (0.5 * (0.5 * (ViscT1 * dx_P1 + ViscT2 * dx_0) / dx1) * dy_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dy_P1
                     ) / dy1
                     + 1 / parameters.flow.Re;

    RealType ViscB = (0.5 * (0.5 * (ViscB1 * dx_P1 + ViscB2 * dx_0) / dx1) * dy_0
                      + 0.5 * (0.5 * (ViscM1 * dx_P1 + ViscM2 * dx_0) / dx1) * dy_M1
                     ) / dy0
                     + 1 / parameters.flow.Re;

    return (ViscT * dvdxT - ViscB * dvdxB) / dy_0;
  }
  inline RealType d2vdy2(const RealType* const lv, const RealType* const lm) {
    const int indexM1 = mapd(0, -1, 0, 1);
    const int index0  = mapd(0, 0, 0, 1);
    const int indexP1 = mapd(0, 1, 0, 1);

    const RealType dy0   = lm[index0];
    const RealType dy1   = lm[indexP1];
    const RealType dySum = dy0 + dy1;
    return 2.0 * (lv[indexP1] / (dy1 * dySum) - lv[index0] / (dy1 * dy0) + lv[indexM1] / (dy0 * dySum));
  }

  inline RealType d2vdy2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    RealType dvdyT = (lv[mapd(0, 1, 0, 1)] - lv[mapd(0, 0, 0, 1)]) / dy_P1;
    RealType dvdyB = (lv[mapd(0, 0, 0, 1)] - lv[mapd(0, -1, 0, 1)]) / dy_0;

    RealType ViscT = lvis[mapd(0, 1, 0, 0)] + 1 / parameters.flow.Re;
    RealType ViscB = lvis[mapd(0, 0, 0, 0)] + 1 / parameters.flow.Re;

    return (ViscT * dvdyT - ViscB * dvdyB) / dy1;
  }

  inline RealType d2vdz2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(1, 0, 0, 2)];
    const RealType dz_M1 = lm[mapd(-1, 0, 0, 2)];
    const RealType dz0   = 0.5 * (dz_0 + dz_M1);
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);

    const RealType ViscR1 = lvis[mapd(0, 0, 1, 0)];
    const RealType ViscR2 = lvis[mapd(0, 1, 1, 0)];
    const RealType ViscM1 = lvis[mapd(0, 0, 0, 0)];
    const RealType ViscM2 = lvis[mapd(0, 1, 0, 0)];
    const RealType ViscL1 = lvis[mapd(0, 0, -1, 0)];
    const RealType ViscL2 = lvis[mapd(0, 1, -1, 0)];

    RealType dvdzR = (lv[mapd(1, 0, 0, 1)] - lv[mapd(0, 0, 0, 1)]) / dz1;
    RealType dvdzL = (lv[mapd(0, 0, 0, 1)] - lv[mapd(-1, 0, 0, 1)]) / dz0;

    RealType ViscR = (0.5 * (0.5 * (ViscR1 * dz_0 + ViscM1 * dz_P1) / (dz1)) * (dy_P1)
                      + 0.5 * (0.5 * (ViscR2 * dz_0 + ViscM2 * dz_P1) / (dz1)) * (dy_0)
                     ) / dy1
                     + 1 / parameters.flow.Re;

    RealType ViscL = (0.5 * (0.5 * (ViscL1 * dz_0 + ViscM1 * dz_M1) / (dz0)) * (dy_P1)
                      + 0.5 * (0.5 * (ViscL2 * dz_0 + ViscM2 * dz_M1) / (dz0)) * (dy_0)
                     ) / dy1
                     + 1 / parameters.flow.Re;

    return (ViscR * dvdzR - ViscL * dvdzL) / dz_0;
  }

  inline RealType d2vdz2(const RealType* const lv, const RealType* const lm) {
    const RealType dz_M1 = lm[mapd(0, 0, -1, 2)];
    const RealType dz_0  = lm[mapd(0, 0, 0, 2)];
    const RealType dz_P1 = lm[mapd(0, 0, 1, 2)];
    const RealType dz0   = 0.5 * (dz_0 + dz_M1);
    const RealType dz1   = 0.5 * (dz_0 + dz_P1);
    const RealType dzSum = dz0 + dz1;
    return 2.0
           * (lv[mapd(0, 0, 1, 1)] / (dz1 * dzSum) - lv[mapd(0, 0, 0, 1)] / (dz1 * dz0) + lv[mapd(0, 0, -1, 1)] / (dz0 * dzSum));
  }

  // Second derivative of the w-component, evaluated at the location of the w-component.
  inline RealType d2wdx2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);
    const RealType dxSum = dx0 + dx1;
    return 2.0
           * (lv[mapd(1, 0, 0, 2)] / (dx1 * dxSum) - lv[mapd(0, 0, 0, 2)] / (dx1 * dx0) + lv[mapd(-1, 0, 0, 2)] / (dx0 * dxSum));
  }

  inline RealType d2wdx2(const RealType* const lv, const RealType* const lm) {
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);
    const RealType dxSum = dx0 + dx1;
    return 2.0
           * (lv[mapd(1, 0, 0, 2)] / (dx1 * dxSum) - lv[mapd(0, 0, 0, 2)] / (dx1 * dx0) + lv[mapd(-1, 0, 0, 2)] / (dx0 * dxSum));
  }

  inline RealType d2wdy2(const RealType* const lv, const RealType* const lm) {
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);
    const RealType dySum = dy0 + dy1;
    return 2.0
           * (lv[mapd(0, 1, 0, 2)] / (dy1 * dySum) - lv[mapd(0, 0, 0, 2)] / (dy1 * dy0) + lv[mapd(0, -1, 0, 2)] / (dy0 * dySum));
  }

  inline RealType d2wdy2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);
    const RealType dySum = dy0 + dy1;
    return 2.0
           * (lv[mapd(0, 1, 0, 2)] / (dy1 * dySum) - lv[mapd(0, 0, 0, 2)] / (dy1 * dy0) + lv[mapd(0, -1, 0, 2)] / (dy0 * dySum));
  }

  inline RealType d2wdz2(const RealType* const lv, const RealType* const lm) {
    const int index_M1 = mapd(0, 0, -1, 2);
    const int index_0  = mapd(0, 0, 0, 2);
    const int index_P1 = mapd(0, 0, 1, 2);

    const RealType dz0   = lm[index_0];
    const RealType dz1   = lm[index_P1];
    const RealType dzSum = dz0 + dz1;
    return 2.0 * (lv[index_P1] / (dz1 * dzSum) - lv[index_0] / (dz1 * dz0) + lv[index_M1] / (dz0 * dzSum));
  }

  inline RealType d2wdz2(
    const RealType* const lv, const RealType* const lvis, const Parameters& parameters, const RealType* const lm
  ) {
    const int index_M1 = mapd(0, 0, -1, 2);
    const int index_0  = mapd(0, 0, 0, 2);
    const int index_P1 = mapd(0, 0, 1, 2);

    const RealType dz0   = lm[index_0];
    const RealType dz1   = lm[index_P1];
    const RealType dzSum = dz0 + dz1;
    return 2.0 * (lv[index_P1] / (dz1 * dzSum) - lv[index_0] / (dz1 * dz0) + lv[index_M1] / (dz0 * dzSum));
  }

  // First derivative of product (u*v), evaluated at the location of the v-component.
  inline RealType duvdx(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 1, 0, 0)]) *
        (lv[mapd(0, 0, 0, 1)] + lv[mapd(1, 0, 0, 1)])) -
        ((lv[mapd(-1, 0, 0, 0)] + lv[mapd(-1, 1, 0, 0)]) *
            (lv[mapd(-1, 0, 0, 1)] + lv[mapd(0, 0, 0, 1)])))
        + parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 1, 0, 0)]) *
            (lv[mapd(0, 0, 0, 1)] - lv[mapd(1, 0, 0, 1)])) -
            (fabs(lv[mapd(-1, 0, 0, 0)] + lv[mapd(-1, 1, 0, 0)]) *
                (lv[mapd(-1, 0, 0, 1)] - lv[mapd(0, 0, 0, 1)])))
        ) / lm[mapd(0, 0, 0, 0)];
#endif

    const RealType hxShort = 0.5 * lm[mapd(0, 0, 0, 0)]; // Distance of corner points in x-direction from center v-value
    const RealType hxLong0 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(-1, 0, 0, 0)]); // Distance between center and west
                                                                                   // v-value
    const RealType hxLong1 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hyShort = 0.5 * lm[mapd(0, 0, 0, 1)]; // Distance of center u-value from upper edge of cell
    const RealType hyLong = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]); // Distance of north and center u-value

    const RealType u00 = lv[mapd(0, 0, 0, 0)];
    const RealType u01 = lv[mapd(0, 1, 0, 0)];
    const RealType v00 = lv[mapd(0, 0, 0, 1)];
    const RealType v10 = lv[mapd(1, 0, 0, 1)];

    const RealType uM10 = lv[mapd(-1, 0, 0, 0)];
    const RealType uM11 = lv[mapd(-1, 1, 0, 0)];
    const RealType vM10 = lv[mapd(-1, 0, 0, 1)];

    // This a central difference expression for the first-derivative. We therefore linearly interpolate u*v onto the
    // surface of the current cell (in 2D: upper left and upper right corner) and then take the central difference.
    const RealType secondOrder = (((hyLong - hyShort) / hyLong * u00 + hyShort / hyLong * u01
                                  ) * ((hxLong1 - hxShort) / hxLong1 * v00 + hxShort / hxLong1 * v10)
                                  - ((hyLong - hyShort) / hyLong * uM10 + hyShort / hyLong * uM11
                                    ) * ((hxLong0 - hxShort) / hxLong0 * v00 + hxShort / hxLong0 * vM10))
                                 / (2.0 * hxShort);

    // This is a forward-difference in donor-cell style. We apply donor cell and again interpolate the velocity values
    // (u-comp.) onto the surface of the cell. We then apply the standard donor cell scheme. This will, however, result
    // in non-equal mesh spacing evaluations (in case of stretched meshes).
    const RealType kr = (hyLong - hyShort) / hyLong * u00 + hyShort / hyLong * u01;
    const RealType kl = (hyLong - hyShort) / hyLong * uM10 + hyShort / hyLong * uM11;

    const RealType firstOrder
      = 1.0 / (4.0 * hxShort)
        * (kr * (v00 + v10) - kl * (vM10 + v00) + fabs(kr) * (v00 - v10) - fabs(kl) * (vM10 - v00));

    // Return linear combination of central and donor-cell difference
    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in duvdx");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. y for u*v at location of u-component. For details on implementation, see duvdx.
  inline RealType duvdy(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 1)] + lv[mapd(1, 0, 0, 1)]) *
        (lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 1, 0, 0)])) -
        ((lv[mapd(0, -1, 0, 1)] + lv[mapd(1, -1, 0, 1)]) *
            (lv[mapd(0, -1, 0, 0)] + lv[mapd(0, 0, 0, 0)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 1)] + lv[mapd(1, 0, 0, 1)]) *
            (lv[mapd(0, 0, 0, 0)] - lv[mapd(0, 1, 0, 0)])) -
            (fabs(lv[mapd(0, -1, 0, 1)] + lv[mapd(1, -1, 0, 1)]) *
                (lv[mapd(0, -1, 0, 0)] - lv[mapd(0, 0, 0, 0)])))) /
        lm[mapd(0, 0, 0, 1)];
#endif

    const RealType hyShort = 0.5 * lm[mapd(0, 0, 0, 1)]; // Distance of corner points in x-direction from center v-value
    const RealType hyLong0 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, -1, 0, 1)]); // Distance between center and west
                                                                                   // v-value
    const RealType hyLong1 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hxShort = 0.5 * lm[mapd(0, 0, 0, 0)]; // Distance of center u-value from upper edge of cell
    const RealType hxLong = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]); // Distance of north and center u-value

    const RealType v00 = lv[mapd(0, 0, 0, 1)];
    const RealType v10 = lv[mapd(1, 0, 0, 1)];
    const RealType u00 = lv[mapd(0, 0, 0, 0)];
    const RealType u01 = lv[mapd(0, 1, 0, 0)];

    const RealType v0M1 = lv[mapd(0, -1, 0, 1)];
    const RealType v1M1 = lv[mapd(1, -1, 0, 1)];
    const RealType u0M1 = lv[mapd(0, -1, 0, 0)];

    const RealType secondOrder = (((hxLong - hxShort) / hxLong * v00 + hxShort / hxLong * v10
                                  ) * ((hyLong1 - hyShort) / hyLong1 * u00 + hyShort / hyLong1 * u01)
                                  - ((hxLong - hxShort) / hxLong * v0M1 + hxShort / hxLong * v1M1
                                    ) * ((hyLong0 - hyShort) / hyLong0 * u00 + hyShort / hyLong0 * u0M1))
                                 / (2.0 * hyShort);

    const RealType kr = (hxLong - hxShort) / hxLong * v00 + hxShort / hxLong * v10;
    const RealType kl = (hxLong - hxShort) / hxLong * v0M1 + hxShort / hxLong * v1M1;

    const RealType firstOrder
      = 1.0 / (4.0 * hyShort)
        * (kr * (u00 + u01) - kl * (u0M1 + u00) + fabs(kr) * (u00 - u01) - fabs(kl) * (u0M1 - u00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in duvdy");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. x for u*w at location of w-component. For details on implementation, see duvdx.
  inline RealType duwdx(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 0, 1, 0)]) *
        (lv[mapd(0, 0, 0, 2)] + lv[mapd(1, 0, 0, 2)])) -
        ((lv[mapd(-1, 0, 0, 0)] + lv[mapd(-1, 0, 1, 0)]) *
            (lv[mapd(-1, 0, 0, 2)] + lv[mapd(0, 0, 0, 2)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 0, 1, 0)]) *
            (lv[mapd(0, 0, 0, 2)] - lv[mapd(1, 0, 0, 2)])) -
            (fabs(lv[mapd(-1, 0, 0, 0)] + lv[mapd(-1, 0, 1, 0)]) *
                (lv[mapd(-1, 0, 0, 2)] - lv[mapd(0, 0, 0, 2)])))) /
        lm[mapd(0, 0, 0, 0)];
#endif

    const RealType hxShort = 0.5 * lm[mapd(0, 0, 0, 0)]; // Distance of corner points in x-direction from center v-value
    const RealType hxLong0 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(-1, 0, 0, 0)]); // Distance between center and west
                                                                                   // v-value
    const RealType hxLong1 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hzShort = 0.5 * lm[mapd(0, 0, 0, 2)]; // Distance of center u-value from upper edge of cell
    const RealType hzLong = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]); // Distance of north and center u-value

    const RealType u00 = lv[mapd(0, 0, 0, 0)];
    const RealType u01 = lv[mapd(0, 0, 1, 0)];
    const RealType w00 = lv[mapd(0, 0, 0, 2)];
    const RealType w10 = lv[mapd(1, 0, 0, 2)];

    const RealType uM10 = lv[mapd(-1, 0, 0, 0)];
    const RealType uM11 = lv[mapd(-1, 0, 1, 0)];
    const RealType wM10 = lv[mapd(-1, 0, 0, 2)];

    const RealType secondOrder = (((hzLong - hzShort) / hzLong * u00 + hzShort / hzLong * u01
                                  ) * ((hxLong1 - hxShort) / hxLong1 * w00 + hxShort / hxLong1 * w10)
                                  - ((hzLong - hzShort) / hzLong * uM10 + hzShort / hzLong * uM11
                                    ) * ((hxLong0 - hxShort) / hxLong0 * w00 + hxShort / hxLong0 * wM10))
                                 / (2.0 * hxShort);

    const RealType kr = (hzLong - hzShort) / hzLong * u00 + hzShort / hzLong * u01;
    const RealType kl = (hzLong - hzShort) / hzLong * uM10 + hzShort / hzLong * uM11;

    const RealType firstOrder
      = 1.0 / (4.0 * hxShort)
        * (kr * (w00 + w10) - kl * (wM10 + w00) + fabs(kr) * (w00 - w10) - fabs(kl) * (wM10 - w00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in duwdx");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. z for u*w at location of u-component. For details on implementation, see duvdx.
  inline RealType duwdz(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 2)] + lv[mapd(1, 0, 0, 2)]) *
        (lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 0, 1, 0)])) -
        ((lv[mapd(0, 0, -1, 2)] + lv[mapd(1, 0, -1, 2)]) *
            (lv[mapd(0, 0, -1, 0)] + lv[mapd(0, 0, 0, 0)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 2)] + lv[mapd(1, 0, 0, 2)]) *
            (lv[mapd(0, 0, 0, 0)] - lv[mapd(0, 0, 1, 0)])) -
            (fabs(lv[mapd(0, 0, -1, 2)] + lv[mapd(1, 0, -1, 2)]) *
                (lv[mapd(0, 0, -1, 0)] - lv[mapd(0, 0, 0, 0)])))) /
        lm[mapd(0, 0, 0, 2)];
#endif

    const RealType hzShort = 0.5 * lm[mapd(0, 0, 0, 2)]; // Distance of corner points in x-direction from center v-value
    const RealType hzLong0 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, -1, 2)]); // Distance between center and west
                                                                                   // v-value
    const RealType hzLong1 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hxShort = 0.5 * lm[mapd(0, 0, 0, 0)]; // Distance of center u-value from upper edge of cell
    const RealType hxLong = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]); // Distance of north and center u-value

    const RealType w00 = lv[mapd(0, 0, 0, 2)];
    const RealType w10 = lv[mapd(1, 0, 0, 2)];
    const RealType u00 = lv[mapd(0, 0, 0, 0)];
    const RealType u01 = lv[mapd(0, 0, 1, 0)];

    const RealType w0M1 = lv[mapd(0, 0, -1, 2)];
    const RealType w1M1 = lv[mapd(1, 0, -1, 2)];
    const RealType u0M1 = lv[mapd(0, 0, -1, 0)];

    const RealType secondOrder = (((hxLong - hxShort) / hxLong * w00 + hxShort / hxLong * w10
                                  ) * ((hzLong1 - hzShort) / hzLong1 * u00 + hzShort / hzLong1 * u01)
                                  - ((hxLong - hxShort) / hxLong * w0M1 + hxShort / hxLong * w1M1
                                    ) * ((hzLong0 - hzShort) / hzLong0 * u00 + hzShort / hzLong0 * u0M1))
                                 / (2.0 * hzShort);

    const RealType kr = (hxLong - hxShort) / hxLong * w00 + hxShort / hxLong * w10;
    const RealType kl = (hxLong - hxShort) / hxLong * w0M1 + hxShort / hxLong * w1M1;

    const RealType firstOrder
      = 1.0 / (4.0 * hzShort)
        * (kr * (u00 + u01) - kl * (u0M1 + u00) + fabs(kr) * (u00 - u01) - fabs(kl) * (u0M1 - u00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in duwdz");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. y for v*w at location of w-component. For details on implementation, see duvdx.
  inline RealType dvwdy(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 0, 1, 1)]) *
        (lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 1, 0, 2)])) -
        ((lv[mapd(0, -1, 0, 1)] + lv[mapd(0, -1, 1, 1)]) *
            (lv[mapd(0, -1, 0, 2)] + lv[mapd(0, 0, 0, 2)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 0, 1, 1)]) *
            (lv[mapd(0, 0, 0, 2)] - lv[mapd(0, 1, 0, 2)])) -
            (fabs(lv[mapd(0, -1, 0, 1)] + lv[mapd(0, -1, 1, 1)]) *
                (lv[mapd(0, -1, 0, 2)] - lv[mapd(0, 0, 0, 2)])))) /
        lm[mapd(0, 0, 0, 1)];
#endif

    const RealType hyShort = 0.5 * lm[mapd(0, 0, 0, 1)]; // Distance of corner points in x-direction from center v-value
    const RealType hyLong0 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, -1, 0, 1)]); // Distance between center and west
                                                                                   // v-value
    const RealType hyLong1 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hzShort = 0.5 * lm[mapd(0, 0, 0, 2)]; // Distance of center u-value from upper edge of cell
    const RealType hzLong = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]); // Distance of north and center u-value

    const RealType v00 = lv[mapd(0, 0, 0, 1)];
    const RealType v01 = lv[mapd(0, 0, 1, 1)];
    const RealType w00 = lv[mapd(0, 0, 0, 2)];
    const RealType w10 = lv[mapd(0, 1, 0, 2)];

    const RealType vM10 = lv[mapd(0, -1, 0, 1)];
    const RealType vM11 = lv[mapd(0, -1, 1, 1)];
    const RealType wM10 = lv[mapd(0, -1, 0, 2)];

    const RealType secondOrder = (((hzLong - hzShort) / hzLong * v00 + hzShort / hzLong * v01
                                  ) * ((hyLong1 - hyShort) / hyLong1 * w00 + hyShort / hyLong1 * w10)
                                  - ((hzLong - hzShort) / hzLong * vM10 + hzShort / hzLong * vM11
                                    ) * ((hyLong0 - hyShort) / hyLong0 * w00 + hyShort / hyLong0 * wM10))
                                 / (2.0 * hyShort);

    const RealType kr = (hzLong - hzShort) / hzLong * v00 + hzShort / hzLong * v01;
    const RealType kl = (hzLong - hzShort) / hzLong * vM10 + hzShort / hzLong * vM11;

    const RealType firstOrder
      = 1.0 / (4.0 * hyShort)
        * (kr * (w00 + w10) - kl * (wM10 + w00) + fabs(kr) * (w00 - w10) - fabs(kl) * (wM10 - w00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in dvwdy");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. z for v*w at location of v-component. For details on implementation, see duvdx.
  inline RealType dvwdz(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 1, 0, 2)]) *
        (lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 0, 1, 1)])) -
        ((lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 1, -1, 2)]) *
            (lv[mapd(0, 0, -1, 1)] + lv[mapd(0, 0, 0, 1)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 1, 0, 2)]) *
            (lv[mapd(0, 0, 0, 1)] - lv[mapd(0, 0, 1, 1)])) -
            (fabs(lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 1, -1, 2)]) *
                (lv[mapd(0, 0, -1, 1)] - lv[mapd(0, 0, 0, 1)])))) /
        lm[mapd(0, 0, 0, 2)];
#endif

    const RealType hzShort = 0.5 * lm[mapd(0, 0, 0, 2)]; // Distance of corner points in x-direction from center v-value
    const RealType hzLong0 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, -1, 2)]); // Distance between center and west
                                                                                   // v-value
    const RealType hzLong1 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hyShort = 0.5 * lm[mapd(0, 0, 0, 1)]; // Distance of center u-value from upper edge of cell
    const RealType hyLong = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]); // Distance of north and center u-value

    const RealType w00 = lv[mapd(0, 0, 0, 2)];
    const RealType w10 = lv[mapd(0, 1, 0, 2)];
    const RealType v00 = lv[mapd(0, 0, 0, 1)];
    const RealType v01 = lv[mapd(0, 0, 1, 1)];

    const RealType w0M1 = lv[mapd(0, 0, -1, 2)];
    const RealType w1M1 = lv[mapd(0, 1, -1, 2)];
    const RealType v0M1 = lv[mapd(0, 0, -1, 1)];

    const RealType secondOrder = (((hyLong - hyShort) / hyLong * w00 + hyShort / hyLong * w10
                                  ) * ((hzLong1 - hzShort) / hzLong1 * v00 + hzShort / hzLong1 * v01)
                                  - ((hyLong - hyShort) / hyLong * w0M1 + hyShort / hyLong * w1M1
                                    ) * ((hzLong0 - hzShort) / hzLong0 * v00 + hzShort / hzLong0 * v0M1))
                                 / (2.0 * hzShort);

    const RealType kr = (hyLong - hyShort) / hyLong * w00 + hyShort / hyLong * w10;
    const RealType kl = (hyLong - hyShort) / hyLong * w0M1 + hyShort / hyLong * w1M1;

    const RealType firstOrder
      = 1.0 / (4.0 * hzShort)
        * (kr * (v00 + v01) - kl * (v0M1 + v00) + fabs(kr) * (v00 - v01) - fabs(kl) * (v0M1 - v00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      spdlog::error("{}, {}", tmp1, tmp2);
      throw std::runtime_error("Error in dvwdz");
    }
#endif

    return tmp2;
  }

  // First derivative of u*u w.r.t. x, evaluated at location of u-component.
  inline RealType du2dx(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 0)] + lv[mapd(1, 0, 0, 0)]) *
        (lv[mapd(0, 0, 0, 0)] + lv[mapd(1, 0, 0, 0)])) -
        ((lv[mapd(-1, 0, 0, 0)] + lv[mapd(0, 0, 0, 0)]) *
            (lv[mapd(-1, 0, 0, 0)] + lv[mapd(0, 0, 0, 0)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 0)] + lv[mapd(1, 0, 0, 0)]) *
            (lv[mapd(0, 0, 0, 0)] - lv[mapd(1, 0, 0, 0)])) -
            (fabs(lv[mapd(-1, 0, 0, 0)] + lv[mapd(0, 0, 0, 0)]) *
                (lv[mapd(-1, 0, 0, 0)] - lv[mapd(0, 0, 0, 0)])))) /
        lm[mapd(0, 0, 0, 0)];
#endif

    const RealType dxShort = 0.5 * lm[mapd(0, 0, 0, 0)];
    // const RealType dxLong0 = 0.5*(lm[mapd(-1,0,0,0)] + lm[mapd(0,0,0,0)]);
    const RealType dxLong1 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);

    const RealType u0  = lv[mapd(0, 0, 0, 0)];
    const RealType uM1 = lv[mapd(-1, 0, 0, 0)];
    const RealType u1  = lv[mapd(1, 0, 0, 0)];

    // const RealType kr = (dxLong1 - dxShort) / dxLong1 * u0 + dxShort / dxLong1 * u1;
    // const RealType kl = (dxLong0 - dxShort) / dxLong0 * u0 + dxShort / dxLong0 * uM1;
    const RealType kr = (u0 + u1) / 2;
    const RealType kl = (u0 + uM1) / 2;

    // Central difference expression which is second-order accurate for uniform meshes. We interpolate u half-way
    // between neighboured u-component values and afterwards build the central difference for u*u.

    /*const RealType secondOrder = (((dxLong1 - dxShort) / dxLong1 * u0 + dxShort / dxLong1 * u1) * ((dxLong1 - dxShort)
       / dxLong1 * u0 + dxShort / dxLong1 * u1)
        - ((dxLong0 - dxShort) / dxLong0 * u0 + dxShort / dxLong0 * uM1) * ((dxLong0 - dxShort) / dxLong0 * u0 + dxShort
       / dxLong0 * uM1) ) / (2.0 * dxShort);*/

    const RealType secondOrder = ((u0 + u1) * (u0 + u1) - (u0 + uM1) * (u0 + uM1)) / (4 * dxLong1);

    // Donor-cell like derivative expression. We evaluate u half-way between neighboured u-components and use this as a
    // prediction of the transport direction.
    const RealType firstOrder = 1.0 / (4.0 * dxShort)
                                * (kr * (u0 + u1) - kl * (uM1 + u0) + fabs(kr) * (u0 - u1) - fabs(kl) * (uM1 - u0));

    // Return linear combination of central- and upwind difference
    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      spdlog::error("{}, {}", tmp1, tmp2);
      throw std::runtime_error("Error in du2dx");
    }
#endif

    return tmp2;
  }

  // First derivative of v*v w.r.t. y, evaluated at location of v-component. For details, see du2dx.
  inline RealType dv2dy(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 1, 0, 1)]) *
        (lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 1, 0, 1)])) -
        ((lv[mapd(0, -1, 0, 1)] + lv[mapd(0, 0, 0, 1)]) *
            (lv[mapd(0, -1, 0, 1)] + lv[mapd(0, 0, 0, 1)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 1, 0, 1)]) *
            (lv[mapd(0, 0, 0, 1)] - lv[mapd(0, 1, 0, 1)])) -
            (fabs(lv[mapd(0, -1, 0, 1)] + lv[mapd(0, 0, 0, 1)]) *
                (lv[mapd(0, -1, 0, 1)] - lv[mapd(0, 0, 0, 1)])))) /
        lm[mapd(0, 0, 0, 1)];
#endif

    const RealType dyShort = 0.5 * lm[mapd(0, 0, 0, 1)];
    // const RealType dyLong0 = 0.5*(lm[mapd(0,-1,0,1)] + lm[mapd(0,0,0,1)]);
    const RealType dyLong1 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);

    const RealType v0  = lv[mapd(0, 0, 0, 1)];
    const RealType vM1 = lv[mapd(0, -1, 0, 1)];
    const RealType v1  = lv[mapd(0, 1, 0, 1)];

    // const RealType kr = (dyLong1 - dyShort) / dyLong1 * v0 + dyShort / dyLong1 * v1;
    // const RealType kl = (dyLong0 - dyShort) / dyLong0 * v0 + dyShort / dyLong0 * vM1;
    const RealType kr = (v0 + v1) / 2;
    const RealType kl = (v0 + vM1) / 2;

    /*const RealType secondOrder = (((dyLong1 - dyShort) / dyLong1 * v0 + dyShort / dyLong1 * v1) * ((dyLong1 - dyShort)
       / dyLong1 * v0 + dyShort / dyLong1 * v1)
        - ((dyLong0 - dyShort) / dyLong0 * v0 + dyShort / dyLong0 * vM1) * ((dyLong0 - dyShort) / dyLong0 * v0 + dyShort
       / dyLong0 * vM1) ) / (2.0 * dyShort);*/

    const RealType secondOrder = ((v0 + v1) * (v0 + v1) - (v0 + vM1) * (v0 + vM1)) / (4 * dyLong1);

    const RealType firstOrder = 1.0 / (4.0 * dyShort)
                                * (kr * (v0 + v1) - kl * (vM1 + v0) + fabs(kr) * (v0 - v1) - fabs(kl) * (vM1 - v0));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      spdlog::error("{}, {}", tmp1, tmp2);
      throw std::runtime_error("Error in dv2dy");
    }
#endif

    return tmp2;
  }

  // First derivative of w*w w.r.t. z, evaluated at location of w-component. For details, see du2dx.
  inline RealType dw2dz(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 0, 1, 2)]) *
        (lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 0, 1, 2)])) -
        ((lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 0, 0, 2)]) *
            (lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 0, 0, 2)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 0, 1, 2)]) *
            (lv[mapd(0, 0, 0, 2)] - lv[mapd(0, 0, 1, 2)])) -
            (fabs(lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 0, 0, 2)]) *
                (lv[mapd(0, 0, -1, 2)] - lv[mapd(0, 0, 0, 2)])))) /
        lm[mapd(0, 0, 0, 2)];
#endif

    const RealType dzShort = 0.5 * lm[mapd(0, 0, 0, 2)];
    // const RealType dzLong0 = 0.5 * (lm[mapd(0, 0, -1, 2)] + lm[mapd(0, 0, 0, 2)]);
    const RealType dzLong1 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);

    const RealType w0  = lv[mapd(0, 0, 0, 2)];
    const RealType wM1 = lv[mapd(0, 0, -1, 2)];
    const RealType w1  = lv[mapd(0, 0, 1, 2)];

    // const RealType kr = (dzLong1 - dzShort) / dzLong1 * w0 + dzShort / dzLong1 * w1;
    // const RealType kl = (dzLong0 - dzShort) / dzLong0 * w0 + dzShort / dzLong0 * wM1;
    const RealType kr = (w0 + w1) / 2;
    const RealType kl = (w0 + wM1) / 2;

    /*const RealType secondOrder = (((dzLong1 - dzShort) / dzLong1 * w0 + dzShort / dzLong1 * w1) * ((dzLong1 - dzShort)
       / dzLong1 * w0 + dzShort / dzLong1 * w1)
        - ((dzLong0 - dzShort) / dzLong0 * w0 + dzShort / dzLong0 * wM1) * ((dzLong0 - dzShort) / dzLong0 * w0 + dzShort
       / dzLong0 * wM1) ) / (2.0 * dzShort);*/

    const RealType secondOrder = ((w0 + w1) * (w0 + w1) - (w0 + wM1) * (w0 + wM1)) / (4 * dzLong1);

    const RealType firstOrder = 1.0 / (4.0 * dzShort)
                                * (kr * (w0 + w1) - kl * (wM1 + w0) + fabs(kr) * (w0 - w1) - fabs(kl) * (wM1 - w0));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      spdlog::error("{}, {}", tmp1, tmp2);
      throw std::runtime_error("Error in dw2dz");
    }
#endif

    return tmp2;
  }

  //************************************************************************************************************************//
  inline RealType fu(const Parameters& parameters, TurbulentFlowFieldKE& flowField, int i, int j) {

    RealType Rd = sqrt(flowField.getk().getScalar(i, j)) * flowField.getheight().getScalar(i, j) * parameters.flow.Re;

    RealType Rt = (flowField.getk().getScalar(i, j) * flowField.getk().getScalar(i, j)) * parameters.flow.Re
                  / (flowField.geteps().getScalar(i, j));

    return (1 - exp(-0.0165 * Rd)) * (1 - exp(-0.0165 * Rd)) * (1 + 20.5 / Rt);
  }

  inline RealType f1(const Parameters& parameters, TurbulentFlowFieldKE& flowField, int i, int j) {
    return 1
           + (0.05 / fu(parameters, flowField, i, j)) * (0.05 / fu(parameters, flowField, i, j))
               * (0.05 / fu(parameters, flowField, i, j));
  }

  inline RealType f2(const Parameters& parameters, TurbulentFlowFieldKE& flowField, int i, int j) {
    RealType Rt = (flowField.getk().getScalar(i, j) * flowField.getk().getScalar(i, j)) * parameters.flow.Re
                  / (flowField.geteps().getScalar(i, j));

    return 1 - exp(-Rt * Rt);
  }

  inline RealType dukdx(
    const RealType* const lv, const Parameters& parameters, const RealType* const lk, const RealType* const lm
  ) {
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType kR = (0.5 * dx_P1 * lk[mapd(0, 0, 0, 0)] + 0.5 * dx_0 * lk[mapd(1, 0, 0, 0)]) / dx1;
    const RealType kL = (0.5 * dx_M1 * lk[mapd(0, 0, 0, 0)] + 0.5 * dx_0 * lk[mapd(-1, 0, 0, 0)]) / dx0;

    return (kR * lv[mapd(0, 0, 0, 0)] - kL * lv[mapd(-1, 0, 0, 0)]) / dx_0
           + (parameters.solver.gamma / dx_0) * (kR * fabs(lv[mapd(0, 0, 0, 0)]) - kL * fabs(lv[mapd(-1, 0, 0, 0)]));
  }

  inline RealType dvkdy(
    const RealType* const lv, const Parameters& parameters, const RealType* const lk, const RealType* const lm
  ) {
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    const RealType kT = (0.5 * dy_P1 * lk[mapd(0, 0, 0, 0)] + 0.5 * dy_0 * lk[mapd(0, 1, 0, 0)]) / dy1;
    const RealType kB = (0.5 * dy_M1 * lk[mapd(0, 0, 0, 0)] + 0.5 * dy_0 * lk[mapd(0, -1, 0, 0)]) / dy0;

    return (kT * lv[mapd(0, 0, 0, 1)] - kB * lv[mapd(0, -1, 0, 1)]) / dy_0
           + (parameters.solver.gamma / dy_0) * (kT * fabs(lv[mapd(0, 0, 0, 1)]) - kB * fabs(lv[mapd(0, -1, 0, 1)]));
  }

  inline RealType dnuTkd2x(
    const RealType* const lv,
    const Parameters&     parameters,
    const RealType* const lvis,
    const RealType* const lk,
    const RealType* const lm
  ) {
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType dkdxR = (lk[mapd(1, 0, 0, 0)] - lk[mapd(0, 0, 0, 0)]) / dx1;
    const RealType dkdxL = (lk[mapd(0, 0, 0, 0)] - lk[mapd(-1, 0, 0, 0)]) / dx0;

    const RealType viscR = (0.5 * dx_P1 * lvis[mapd(0, 0, 0, 0)] + 0.5 * dx_0 * lvis[mapd(1, 0, 0, 0)]) / dx1;
    const RealType viscL = (0.5 * dx_M1 * lvis[mapd(0, 0, 0, 0)] + 0.5 * dx_0 * lvis[mapd(-1, 0, 0, 0)]) / dx0;

    return (viscR * dkdxR - viscL * dkdxL) / dx_0;
  }

  inline RealType dnuTkd2y(
    const RealType* const lv,
    const Parameters&     parameters,
    const RealType* const lvis,
    const RealType* const lk,
    const RealType* const lm
  ) {
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    const RealType dkdyT = (lk[mapd(0, 1, 0, 0)] - lk[mapd(0, 0, 0, 0)]) / dy1;
    const RealType dkdyB = (lk[mapd(0, 0, 0, 0)] - lk[mapd(0, -1, 0, 0)]) / dy0;

    const RealType viscT = (0.5 * dy_P1 * lvis[mapd(0, 0, 0, 0)] + 0.5 * dy_0 * lvis[mapd(0, 1, 0, 0)]) / dy1;
    const RealType viscB = (0.5 * dy_M1 * lvis[mapd(0, 0, 0, 0)] + 0.5 * dy_0 * lvis[mapd(0, -1, 0, 0)]) / dy0;

    return (viscT * dkdyT - viscB * dkdyB) / dy_0;
  }

  inline RealType duepsdx(
    const RealType* const lv,
    const Parameters&     parameters,
    const RealType* const leps,
    const RealType* const lk,
    const RealType* const lm
  ) {
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType lepsR = (0.5 * dx_P1 * lk[mapd(0, 0, 0, 0)] + 0.5 * dx_0 * lk[mapd(1, 0, 0, 0)]) / dx1;
    const RealType lepsL = (0.5 * dx_M1 * lk[mapd(0, 0, 0, 0)] + 0.5 * dx_0 * lk[mapd(-1, 0, 0, 0)]) / dx0;

    return (lepsR * lv[mapd(0, 0, 0, 0)] - lepsL * lv[mapd(-1, 0, 0, 0)]) / dx_0
           + (parameters.solver.gamma / dx_0
             ) * (lepsR * fabs(lv[mapd(0, 0, 0, 0)]) - lepsL * fabs(lv[mapd(-1, 0, 0, 0)]));
  }

  inline RealType dvepsdy(
    const RealType* const lv,
    const Parameters&     parameters,
    const RealType* const leps,
    const RealType* const lk,
    const RealType* const lm
  ) {
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    const RealType lepsT = (0.5 * dy_P1 * lk[mapd(0, 0, 0, 0)] + 0.5 * dy_0 * lk[mapd(0, 1, 0, 0)]) / dy1;
    const RealType lepsB = (0.5 * dy_M1 * lk[mapd(0, 0, 0, 0)] + 0.5 * dy_0 * lk[mapd(0, -1, 0, 0)]) / dy0;

    return (lepsT * lv[mapd(0, 0, 0, 1)] - lepsB * lv[mapd(0, -1, 0, 1)]) / dy_0
           + (parameters.solver.gamma / dy_0
             ) * (lepsT * fabs(lv[mapd(0, 0, 0, 1)]) - lepsB * fabs(lv[mapd(0, -1, 0, 1)]));
  }

  inline RealType dfunuTepsd2x(
    TurbulentFlowFieldKE& flowField,
    const RealType* const lv,
    const Parameters&     parameters,
    const RealType* const lvis,
    const RealType* const leps,
    const RealType* const lm,
    int                   i,
    int                   j
  ) {
    const RealType dx_0  = lm[mapd(0, 0, 0, 0)];
    const RealType dx_P1 = lm[mapd(1, 0, 0, 0)];
    const RealType dx_M1 = lm[mapd(-1, 0, 0, 0)];
    const RealType dx0   = 0.5 * (dx_0 + dx_M1);
    const RealType dx1   = 0.5 * (dx_0 + dx_P1);

    const RealType depsdxR = (leps[mapd(1, 0, 0, 0)] - leps[mapd(0, 0, 0, 0)]) / dx1;
    const RealType depsdxL = (leps[mapd(0, 0, 0, 0)] - leps[mapd(-1, 0, 0, 0)]) / dx0;

    const RealType viscR = (0.5 * dx_P1 * lvis[mapd(0, 0, 0, 0)] + 0.5 * dx_0 * lvis[mapd(1, 0, 0, 0)]) / dx1;
    const RealType viscL = (0.5 * dx_M1 * lvis[mapd(0, 0, 0, 0)] + 0.5 * dx_0 * lvis[mapd(-1, 0, 0, 0)]) / dx0;

    const RealType fuR = fu(parameters, flowField, i + 1, j);
    const RealType fuL = fu(parameters, flowField, i - 1, j);

    return (fuR * viscR * depsdxR - fuL * viscL * depsdxL) / dx_0;
  }

  inline RealType dfunuTepsd2y(
    TurbulentFlowFieldKE& flowField,
    const RealType* const lv,
    const Parameters&     parameters,
    const RealType* const lvis,
    const RealType* const leps,
    const RealType* const lm,
    int                   i,
    int                   j
  ) {
    const RealType dy_0  = lm[mapd(0, 0, 0, 1)];
    const RealType dy_P1 = lm[mapd(0, 1, 0, 1)];
    const RealType dy_M1 = lm[mapd(0, -1, 0, 1)];
    const RealType dy0   = 0.5 * (dy_0 + dy_M1);
    const RealType dy1   = 0.5 * (dy_0 + dy_P1);

    const RealType depsdyT = (leps[mapd(0, 1, 0, 0)] - leps[mapd(0, 0, 0, 0)]) / dy1;
    const RealType depsdyB = (leps[mapd(0, 0, 0, 0)] - leps[mapd(0, -1, 0, 0)]) / dy0;

    const RealType viscT = (0.5 * dy_P1 * lvis[mapd(0, 0, 0, 0)] + 0.5 * dy_0 * lvis[mapd(0, 1, 0, 0)]) / dy1;
    const RealType viscB = (0.5 * dy_M1 * lvis[mapd(0, 0, 0, 0)] + 0.5 * dy_0 * lvis[mapd(0, -1, 0, 0)]) / dy0;

    const RealType fuT = fu(parameters, flowField, i, j + 1);
    const RealType fuB = fu(parameters, flowField, i, j - 1);

    return (fuT * viscT * depsdyT - fuB * viscB * depsdyB) / dy_0;
  }

  //************************************************************************************************************************//

  // make if tree for FGH turbulent stencil with return values.
  inline RealType computeF2D(
    const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    return localVelocity[mapd(0, 0, 0, 0)]
        + dt * (1 / parameters.flow.Re * (d2udx2(localVelocity, localMeshsize)
            + d2udy2(localVelocity, localMeshsize)) - du2dx(localVelocity, parameters, localMeshsize)
            - duvdy(localVelocity, parameters, localMeshsize) + parameters.environment.gx);
  }

  inline RealType computeF2D_turbulent(
    const RealType* const localVelocity,
    const RealType* const localViscosity,
    const RealType* const localMeshsize,
    const Parameters&     parameters,
    RealType              dt
  ) {
    return localVelocity[mapd(0, 0, 0, 0)]
        + dt * (2*d2udx2(localVelocity, localViscosity, parameters, localMeshsize)
            + d2udy2(localVelocity, localViscosity, parameters, localMeshsize) + d2vdydx(localVelocity, localViscosity, parameters, localMeshsize) - du2dx(localVelocity, parameters, localMeshsize)
            - duvdy(localVelocity, parameters, localMeshsize) + parameters.environment.gx);
  }

  inline RealType computeF2D_turbulent_KE(
    const RealType* const localVelocity,
    const RealType* const localViscosity,
    const RealType* const localKineticEnergy,
    const RealType* const localMeshsize,
    const Parameters&     parameters,
    RealType              dt
  ) {
    return localVelocity[mapd(0, 0, 0, 0)]
        + dt * (2*d2udx2(localVelocity, localViscosity, parameters, localMeshsize)
            + d2udy2(localVelocity, localViscosity, parameters, localMeshsize) + d2vdydx(localVelocity, localViscosity, parameters, localMeshsize) - du2dx(localVelocity, parameters, localMeshsize)
            - duvdy(localVelocity, parameters, localMeshsize) - 2*dkdx(localKineticEnergy, localMeshsize)/3 + parameters.environment.gx);
  }

  inline RealType computeG2D_turbulent_KE(
    const RealType* const localVelocity,
    const RealType* const localViscosity,
    const RealType* const localKineticEnergy,
    const RealType* const localMeshsize,
    const Parameters&     parameters,
    RealType              dt
  ) {
    return localVelocity[mapd(0, 0, 0, 1)]
        + dt * ( d2vdx2(localVelocity, localViscosity, parameters, localMeshsize) + d2udxdy(localVelocity, localViscosity, parameters, localMeshsize)
            + 2*d2vdy2(localVelocity, localViscosity, parameters, localMeshsize) - duvdx(localVelocity, parameters, localMeshsize)
            - dv2dy(localVelocity, parameters, localMeshsize) - 2*dkdy(localKineticEnergy, localMeshsize)/3 +  + parameters.environment.gy);
  }

  inline RealType computeG2D(
    const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    return localVelocity[mapd(0, 0, 0, 1)]
        + dt * (1 / parameters.flow.Re * (d2vdx2(localVelocity, localMeshsize)
            + d2vdy2(localVelocity, localMeshsize)) - duvdx(localVelocity, parameters, localMeshsize)
            - dv2dy(localVelocity, parameters, localMeshsize) + parameters.environment.gy);
  }

  inline RealType computeG2D_turbulent(
    const RealType* const localVelocity,
    const RealType* const localViscosity,
    const RealType* const localMeshsize,
    const Parameters&     parameters,
    RealType              dt
  ) {
    return localVelocity[mapd(0, 0, 0, 1)]
        + dt * ( d2vdx2(localVelocity, localViscosity, parameters, localMeshsize) + d2udxdy(localVelocity, localViscosity, parameters, localMeshsize)
            + 2*d2vdy2(localVelocity, localViscosity, parameters, localMeshsize) - duvdx(localVelocity, parameters, localMeshsize)
            - dv2dy(localVelocity, parameters, localMeshsize) + parameters.environment.gy);
  }

  inline RealType computek2D(
    TurbulentFlowFieldKE& flowField,
    const RealType* const localVelocity,
    const RealType* const localViscosity,
    const RealType* const localk,
    const RealType* const localMeshsize,
    const Parameters&     parameters,
    RealType              dt,
    int                   i,
    int                   j
  ) {
    // std::cout
    //   << "Hello from "
    //      "k2d**************************************************************************************************"
    //      "**************\n";
    // std::cout << "dnuTkd2x = " << dnuTkd2x(localVelocity, parameters, localViscosity, localk, localMeshsize) << "\n";
    // std::cout << "dnuTkd2y = " << dnuTkd2y(localVelocity, parameters, localViscosity, localk, localMeshsize) << "\n";
    // std::cout << "dukdx = " << dukdx(localVelocity, parameters, localk, localMeshsize) << "\n";
    // std::cout << "dvkdy = " << dvkdy(localVelocity, parameters, localk, localMeshsize) << "\n";
    return localk[mapd(0, 0, 0, 0)]
           + dt
               * (dnuTkd2x(localVelocity, parameters, localViscosity, localk, localMeshsize) 
                + dnuTkd2y(localVelocity, parameters, localViscosity, localk, localMeshsize) 
                - dukdx(localVelocity, parameters, localk, localMeshsize) 
                - dvkdy(localVelocity, parameters, localk, localMeshsize) 
                + 0.5 * flowField.getnuT().getScalar(i,j) * (4 * pow(dudx(localVelocity, localMeshsize),2) + 2 * pow(dudy(localVelocity, localMeshsize) + dvdx(localVelocity, localMeshsize),2) + 4 * pow(dvdy(localVelocity, localMeshsize),2)) 
                - flowField.geteps().getScalar(i, j));
  }

  inline RealType computeEpsilon2D(
    TurbulentFlowFieldKE& flowField,
    const RealType* const localVelocity,
    const RealType* const localViscosity,
    const RealType* const localk,
    const RealType* const localEpsilon,
    const RealType* const localMeshsize,
    const Parameters&     parameters,
    RealType              dt,
    int                   i,
    int                   j
  ) {
    return localEpsilon[mapd(0, 0, 0, 0)]
           + dt
               * ((parameters.turbulent.ce / parameters.turbulent.cmu) 
               * (dfunuTepsd2x(flowField, localVelocity, parameters, localViscosity, localEpsilon, localMeshsize, i, j) 
                + dfunuTepsd2y(flowField, localVelocity, parameters, localViscosity, localEpsilon, localMeshsize,i, j)) 
                - duepsdx(localVelocity, parameters, localEpsilon, localk, localMeshsize) 
                - dvepsdy(localVelocity, parameters, localEpsilon, localk, localMeshsize)
                + 0.5*parameters.turbulent.c1*flowField.getk().getScalar(i,j)*f1(parameters, flowField,i,j)*(4 * pow(dudx(localVelocity, localMeshsize),2) + 2 * pow(dudy(localVelocity, localMeshsize) + dvdx(localVelocity, localMeshsize),2) + 4 * pow(dvdy(localVelocity, localMeshsize) ,2))
                - parameters.turbulent.c2*f2(parameters, flowField,i,j)*pow(flowField.geteps().getScalar(i,j),2)/flowField.getk().getScalar(i,j));
  }

  inline RealType computeF3D(
    const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    return localVelocity[mapd(0, 0, 0, 0)]
        + dt * (1 / parameters.flow.Re * (d2udx2(localVelocity, localMeshsize)
            + d2udy2(localVelocity, localMeshsize) + d2udz2(localVelocity, localMeshsize))
            - du2dx(localVelocity, parameters, localMeshsize) - duvdy(localVelocity, parameters, localMeshsize)
            - duwdz(localVelocity, parameters, localMeshsize) + parameters.environment.gx);
  }

  inline RealType computeF3D_turbulent(
    const RealType* const localVelocity,
    const RealType* const localViscosity,
    const RealType* const localMeshsize,
    const Parameters&     parameters,
    RealType              dt
  ) {
    return localVelocity[mapd(0, 0, 0, 0)]
        + dt * (2*d2udx2(localVelocity, localViscosity, parameters, localMeshsize)
            + d2udy2(localVelocity, localViscosity, parameters, localMeshsize) 
            + d2udz2(localVelocity, localViscosity, parameters, localMeshsize)  
            + d2vdydx(localVelocity, localViscosity, parameters, localMeshsize) 
            + d2wdxdz(localVelocity, localViscosity, parameters, localMeshsize)
            - du2dx(localVelocity, parameters, localMeshsize) - duvdy(localVelocity, parameters, localMeshsize)
            - duwdz(localVelocity, parameters, localMeshsize) + parameters.environment.gx);
  }

  inline RealType computeG3D(
    const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    return localVelocity[mapd(0, 0, 0, 1)]
        + dt * (1 / parameters.flow.Re * (d2vdx2(localVelocity, localMeshsize)
            + d2vdy2(localVelocity, localMeshsize) + d2vdz2(localVelocity, localMeshsize))
            - dv2dy(localVelocity, parameters, localMeshsize) - duvdx(localVelocity, parameters, localMeshsize)
            - dvwdz(localVelocity, parameters, localMeshsize) + parameters.environment.gy);
  }

  inline RealType computeG3D_turbulent(
    const RealType* const localVelocity,
    const RealType* const localViscosity,
    const RealType* const localMeshsize,
    const Parameters&     parameters,
    RealType              dt
  ) {
    return localVelocity[mapd(0, 0, 0, 1)]
        + dt * (d2vdx2(localVelocity, localViscosity, parameters, localMeshsize)            
              + 2*d2vdy2(localVelocity, localViscosity, parameters, localMeshsize)
              + d2vdz2(localVelocity, localViscosity, parameters, localMeshsize)
              + d2udxdy(localVelocity, localViscosity, parameters, localMeshsize)
              + d2wdydz(localVelocity, localViscosity, parameters, localMeshsize)
              - dv2dy(localVelocity, parameters, localMeshsize) - duvdx(localVelocity, parameters, localMeshsize)
              - dvwdz(localVelocity, parameters, localMeshsize) + parameters.environment.gy);
  }

  inline RealType computeH3D(
    const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    return localVelocity[mapd(0, 0, 0, 2)]
        + dt * (1 / parameters.flow.Re * (d2wdx2(localVelocity, localMeshsize)
            + d2wdy2(localVelocity, localMeshsize) + d2wdz2(localVelocity, localMeshsize))
            - dw2dz(localVelocity, parameters, localMeshsize) - duwdx(localVelocity, parameters, localMeshsize)
            - dvwdy(localVelocity, parameters, localMeshsize) + parameters.environment.gz);
  }

  inline RealType computeH3D_turbulent(
    const RealType* const localVelocity,
    const RealType* const localViscosity,
    const RealType* const localMeshsize,
    const Parameters&     parameters,
    RealType              dt
  ) {
    return localVelocity[mapd(0, 0, 0, 1)]
        + dt * (d2wdx2(localVelocity, localViscosity, parameters, localMeshsize)            
              + d2wdy2(localVelocity, localViscosity, parameters, localMeshsize)
              + 2*d2wdz2(localVelocity, localViscosity, parameters, localMeshsize)
              + d2udzdx(localVelocity, localViscosity, parameters, localMeshsize)
              + d2vdzdy(localVelocity, localViscosity, parameters, localMeshsize)
              - dw2dz(localVelocity, parameters, localMeshsize) - duwdx(localVelocity, parameters, localMeshsize)
              - dvwdy(localVelocity, parameters, localMeshsize) + parameters.environment.gz);
  }

} // namespace Stencils
