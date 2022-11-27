#pragma once
// #include "Parameters.hpp"

// #include "Stencils/BoundaryStencil.hpp"
#include "Iterators.hpp"

#include "Stencils/PressureBufferFillStencil.hpp"
#include "Stencils/PressureBufferReadStencil.hpp"
#include "Stencils/VelocityBufferFillStencil.hpp"
#include "Stencils/VelocityBufferReadStencil.hpp"
// #include <mpi.h>

namespace ParallelManagers {
  class PetscParallelManager {
  private:
    Parameters& parameters_;
    FlowField&  flowfield_;

  public:
    PetscParallelManager(Parameters& parameters, FlowField& flowfield);

    ~PetscParallelManager() = default;
    void communicatePressure();
    void communicateVelocities();
  };
} // namespace ParallelManagers
