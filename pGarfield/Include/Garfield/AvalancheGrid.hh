#ifndef G_AVALANCHE_GRID_H
#define G_AVALANCHE_GRID_H

#include <iostream>
#include <string>
#include <vector>

#include "AvalancheMicroscopic.hh"
#include "Garfield/ComponentParallelPlate.hh"
#include "GarfieldConstants.hh"
#include "Sensor.hh"

namespace Garfield {

/// Calculate avalanches in a uniform electric field using avalanche
/// statistics.

class AvalancheGrid {
public:
  /// Constructor
  AvalancheGrid() {}
  /// Destructor
  ~AvalancheGrid() {}
  /// Set the sensor.
  void SetSensor(Sensor *sensor) { m_sensor = sensor; }
  /// Set the AvalancheMicroscopic.
  void SetAvalancheMicroscopic(AvalancheMicroscopic *avmc) { m_avmc = avmc; }

  /** Start grid based avalanche simulation.
   *
   * \param zmin,zmax z-coordinate range of grid [cm].
   * \param zsteps amount of z-coordinate points in grid.
   * \param xmin,xmax x-coordinate range of grid [cm].
   * \param xsteps amount of x-coordinate points in grid.
   */
  void StartGridAvalanche();
  /// Set the electron drift velocity (in cm / ns).
  void SetElectronVelocity(const double vx, const double vy, const double vz) {
    double vel = sqrt(vx * vx + vy * vy + vz * vz);
    if (vel != std::abs(vx) && vel != std::abs(vy) && vel != std::abs(vz))
      return;
    int nx = (int)vx / vel;
    int ny = (int)vy / vel;
    int nz = (int)vz / vel;
    m_velNormal = {nx, ny, nz};
    m_Velocity = -std::abs(vel);
  }
  /// Set the electron Townsend coefficient (in 1 / cm).
  void SetElectronTownsend(const double town) { m_Townsend = town; }
  /// Set the electron attachment coefficient (in 1 / cm).
  void SetElectronAttachment(const double att) { m_Attachment = att; }
  /// Set the maximum avalanche size (1e7 by default).
  void SetMaxAvalancheSize(const double size) { m_MaxSize = size; }
  /// Enable transverse diffusion of electrons with transverse diffusion
  /// coefficients (in √cm).
  void EnableDiffusion(const double diffSigma) {
    m_diffusion = true;
    m_DiffSigma = diffSigma;
  }
  /** Setting the starting point of an electron that .
   *
   * \param z z-coordinate of initial electron.
   * \param x x-coordinate of initial electron.
   * \param v speed of initial electron.
   * \param t starting time of avalanche.
   */
  void CreateAvalanche(const double x, const double y, const double z,
                       const double t = 0, const int n = 1);
  /// Import electron data from AvalancheMicroscopic class
  void GetElectronsFromAvalancheMicroscopic();

  /// Import electron data from AvalancheMicroscopic class
  void SetGrid(const double xmin, const double xmax, const int xsteps,
               const double ymin, const double ymax, const int ysteps,
               const double zmin, const double zmax, const int zsteps);

  /// Returns the initial number of electrons in the avalanche.
  int GetAmountOfStartingElectrons() { return m_nestart; }
  /// Returns the final number of electrons in the avalanche.
  int GetAvalancheSize() { return m_avgrid.N; }

  /// Asigning layer index to all Avalanche nodes.
  void AsignLayerIndex(ComponentParallelPlate *RPC);

  void EnableDebugging() { m_debug = true; }

  void Reset();

private:
  bool m_debug = false;

  double m_Townsend = -1; // [1/cm]

  double m_Attachment = -1; // [1/cm]

  double m_Velocity = 0.; // [cm/ns]

  std::vector<int> m_velNormal = {0, 0, 0};

  double m_MaxSize = 1.6e7; // Saturations size

  bool m_Saturated = false; // Check if avalanche has reached maximum size

  double m_SaturationTime =
      -1.; // Time when the avalanche has reached maximum size

  bool m_diffusion = false; // Check if transverse diffusion is enabled.

  double m_DiffSigma = 0.; // Transverse diffusion coefficients (in √cm).

  int m_nestart = 0.;

  bool m_driftAvalanche = false;
  bool m_importAvalanche = false;

  bool m_layerIndix = false;
  std::vector<double> m_NLayer;

  std::string m_className = "AvalancheGrid";

  Sensor *m_sensor = nullptr;

  AvalancheMicroscopic *m_avmc = nullptr;

  bool m_printPar = false;

  struct Grid {
    std::vector<double> zgrid; ///< Grid points of z-coordinate.
    int zsteps = 0;            ///< Amount of grid points.
    double zStepSize =
        0.; ///< Distance between the grid points of z-coordinate.

    std::vector<double> ygrid; ///< Grid points of y-coordinate.
    double yStepSize = 0.;     ///< Amount of grid points.
    int ysteps = 0.; ///< Distance between the grid points of y-coordinate.

    std::vector<double> xgrid; ///< Grid points of x-coordinate.
    double xStepSize = 0.;     ///< Amount of grid points.
    int xsteps = 0.; ///< Distance between the grid points of x-coordinate.

    bool gridset = false; ///< Keeps track if the grid has been defined.
    int N = 0;            ///< Total amount of charge.

    double time = 0; ///< Clock.

    bool run = true; ///< Tracking if the charges are still in the drift gap.
  };

  struct AvalancheNode {
    double ix = 0;
    double iy = 0;
    double iz = 0;

    int n = 1;

    int layer = 0;

    double townsend = 0;
    double attachment = 0;
    double velocity = 0;

    double stepSize = 0;
    std::vector<int> velNormal = {0, 0, 0};

    double time = 0.; ///< Clock.
    double dt = -1.;  ///< time step.

    bool active = true;
    double dSigmaL = 0;
    double dSigmaT = 0;
  };

  std::vector<AvalancheNode> m_activeNodes = {};

  Grid m_avgrid;
  // Setting z-coordinate grid.
  void SetZGrid(Grid &av, const double top, const double bottom,
                const int steps);
  // Setting y-coordinate grid.
  void SetYGrid(Grid &av, const double top, const double bottom,
                const int steps);
  // Setting x-coordinate grid.
  void SetXGrid(Grid &av, const double top, const double bottom,
                const int steps);
  // Get size of avalanche when going from z to z-dz.
  int GetAvalancheSize(double dz, const int nsize, const double alpha,
                       const double eta);
  // Assign electron to the closest grid point.
  bool SnapToGrid(Grid &av, const double x, const double y, const double z,
                  const double v, const int n = 1);
  // Go to next time step.
  void NextAvalancheGridPoint(Grid &av);
  // Obtain the Townsend coef., Attachment coef. and velocity vector from
  // sensor class.
  bool GetParameters(AvalancheNode &newNode);

  void DeactivateNode(AvalancheNode &node);
};
} // namespace Garfield

#endif
