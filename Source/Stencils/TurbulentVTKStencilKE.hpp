#include <fstream>
#include <memory>

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowFieldKE.hpp"

namespace Stencils {

  /** Stencil for writting VTK files
   *
   * When iterated with, creates a VTK file.
   */
  class TurbulentVTKStencilKE: public FieldStencil<TurbulentFlowFieldKE> {
  private:
    bool          written_; //! Whether the file has already been written
    std::string   prefix_;  //! Prefix to be attached to the vtk files
    std::ofstream ofile_;   //! Output file stream

    std::stringstream pressureStream_; //! Stream for the pressure data
    std::stringstream velocityStream_; //! Stream for the velocity data
    std::stringstream nuTStream_;
    std::stringstream kStream_;
    std::stringstream epsilonStream_;

    void writeVTKHeader(std::ostream& file) const;
    void writePoints(std::ostream& file, RealType simulationTime) const;

    /** Open output file
     * Opens the output file and prepares for writing.
     * It also writes the header of the file and the points of the grid.
     * @param timeStep Current time step of the simulation
     * @param simulationTime Current simulation time in double format
     */
    void openFile(int timeStep, RealType simulationTime);

    /** Finish writing. Must be called once the file has been written.
     *
     * Stores all the streams and closes the file.
     */
    void closeFile();

  public:
    TurbulentVTKStencilKE(const Parameters& parameters);
    ~TurbulentVTKStencilKE() override = default;

    void apply(TurbulentFlowFieldKE& flowField, int i, int j) override;
    void apply(TurbulentFlowFieldKE& flowField, int i, int j, int k) override;

    void write(TurbulentFlowFieldKE& flowField, int timeStep, RealType simulationTime);
  };
} // namespace Stencils
