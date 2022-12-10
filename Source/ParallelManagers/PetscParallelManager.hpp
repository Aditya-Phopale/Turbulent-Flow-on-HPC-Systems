#pragma once
// #include "Parameters.hpp"

// #include "Stencils/BoundaryStencil.hpp"
#include "Iterators.hpp"

#include "Stencils/PressureBufferFillStencil.cpph"
#include "Stencils/PressureBufferFillStencil.hpp"
#include "Stencils/PressureBufferReadStencil.cpph"
#include "Stencils/PressureBufferReadStencil.hpp"
#include "Stencils/VelocityBufferFillStencil.cpph"
#include "Stencils/VelocityBufferFillStencil.hpp"
#include "Stencils/VelocityBufferReadStencil.cpph"
#include "Stencils/VelocityBufferReadStencil.hpp"
#include "Stencils/ViscosityBufferFillStencil.hpp"
#include "Stencils/ViscosityBufferReadStencil.hpp"

// #include <mpi.h>

namespace ParallelManagers {

  template <class FlowFieldType>
  class PetscParallelManager {
  private:
    Parameters&    parameters_;
    FlowFieldType& flowfield_;

  public:
    PetscParallelManager(Parameters& parameters, FlowFieldType& flowfield);

    ~PetscParallelManager() = default;
    void communicatePressure();
    void communicateVelocities();
    void communicateViscosity();
  };
} // namespace ParallelManagers
