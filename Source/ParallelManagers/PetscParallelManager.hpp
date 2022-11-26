#include "Parameters.hpp"

#include "Stencils/BoundaryStencil.hpp"
#include "Stencils/PressureBufferFillStencil.hpp"
#include "Stencils/PressureBufferReadStencil.hpp"
#include "Stencils/VelocityBufferFillStencil.hpp"
#include "Stencils/VelocityBufferReadStencil.hpp"
#include "mpi.h"

namespace ParallelManagers {
  class PetscParallelManager {
  private:
    Parameters parameters_;
    FlowField  flowfield_;
    // Stencils::PressureBufferFillStencil pfill_;
    // Stencils::PressureBufferReadStencil pread_;
    // Stencils::VelocityBufferFillStencil vfill_;
    // Stencils::VelocityBufferReadStencil vread_;

    // ParallelBoundaryIterator<FlowField> pressureFillIterator;
    // ParallelBoundaryIterator<FlowField> pressureReadIterator;
    // ParallelBoundaryIterator<FlowField> velocityFillIterator;
    // ParallelBoundaryIterator<FlowField> velocityReadIterator;

  public:
    PetscParallelManager(const Parameters& parameters, FlowField& flowfield):
      parameters_(parameters),
      flowfield_(flowfield) {}
    ~PetscParallelManager() = default;
    void communicatePressure();
    void communicateVelocities();
  };
} // namespace ParallelManagers
