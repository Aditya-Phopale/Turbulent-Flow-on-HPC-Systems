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
    Stencils::PressureBufferFillStencil pressureFillStencil_;
    Stencils::PressureBufferReadStencil pressureReadStencil_;
    Stencils::VelocityBufferFillStencil velocityFillStencil_;
    Stencils::VelocityBufferReadStencil VeclocityBufferStencil_;
    ParallelBoundaryIterator<FlowField> pressureFillIterator_;
    ParallelBoundaryIterator<FlowField> pressureReadIterator_;
    ParallelBoundaryIterator<FlowField> VelocityFillIterator_;
    ParallelBoundaryIterator<FlowField> VelocityReadIterator_;

  public:
    PetscParallelManager(Parameters& parameters, FlowField& flowfield);

    ~PetscParallelManager() = default;
    void communicatePressure();
    void communicateVelocities();
  };
} // namespace ParallelManagers


// I am still at pressurefill stencil need to update all the files