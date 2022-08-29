#ifndef G_COMPONENT_PP_H
#define G_COMPONENT_PP_H

#include <string>

#include "Component.hh"
#include "ComponentGrid.hh"
#include "Medium.hh"

#include <TF1.h>
#include <TF2.h>

namespace Garfield {

/// Component for parallel-plate geometries.

class ComponentParallelPlate : public Component {
public:
  /// Constructor
  ComponentParallelPlate();
  /// Destructor
  ~ComponentParallelPlate() {}

  /** Define the geometry.
   * \param N amount of layers in the geometry, this includes the gas gaps
   * \f$y\f$. \param d thickness of the layers starting from the bottom to the
   * top lauer along \f$y\f$. \param eps relative permittivities of the layers
   * starting from the bottom to the top lauer along \f$y\f$ . Here, the  gas
   * gaps having a value of 1. \param sigmaIndex Indices of the resistive
   * layers. \param V applied potential difference between the parallel plates.
   */
  void Setup(const int N, std::vector<double> eps, std::vector<double> d,
             const double V, std::vector<int> sigmaIndex = {});

  void ElectricField(const double x, const double y, const double z, double &ex,
                     double &ey, double &ez, Medium *&m, int &status) override;
  void ElectricField(const double x, const double y, const double z, double &ex,
                     double &ey, double &ez, double &v, Medium *&m,
                     int &status) override;

  double WeightingPotential(const double x, const double y, const double z,
                            const std::string &label) override;

  bool GetVoltageRange(double &vmin, double &vmax) override;

  /** Add a pixel electrode.
   * \param x,z position of the center of the electrode in the xz-plane.
   * \param lx width in the along \f$x\f$.
   * \param lz width in the along \f$z\f$.
   * \param label give name using a string.
   * \param fromAnode is \f$true\f$ is the electrode is the andode and
   * \f$false\f$ if it is the cathode.
   */
  void AddPixel(double x, double z, double lx, double lz,
                const std::string &label, bool fromAnode = true);
  /// Add strip electrode.
  void AddStrip(double z, double lz, const std::string &label,
                bool fromAnode = true);

  /// Add plane electrode, if you want to read the signal from the cathode set
  /// the second argument to false.
  void AddPlane(const std::string &label, bool fromAnode = true);

  /// Setting the medium
  void SetMedium(Medium *medium) { m_medium = medium; }

  /** Calculate time-dependent weighting potential on a grid.
   * \param xmin,ymin,zmin minimum value of the interval in the \f$x\f$-,
   * \f$y\f$- and \f$z\f$-direction. \param xmax,ymax,zmax maximum value of the
   * interval in the \f$x\f$-,\f$y\f$- and \f$z\f$-direction. \param
   * xsteps,ysteps,zsteps mumber of grid nodes in the \f$x\f$-,\f$y\f$- and
   * \f$z\f$-direction. \param label give name using a string.
   */
  void SetWeightingPotentialGrid(const double xmin, const double xmax,
                                 const double xsteps, const double ymin,
                                 const double ymax, const double ysteps,
                                 const double zmin, const double zmax,
                                 const double zsteps, const std::string &label);

  /// This will calculate all electrodes time-dependent weighting potential on
  /// the specified grid.
  void SetWeightingPotentialGrids(const double xmin, const double xmax,
                                  const double xsteps, const double ymin,
                                  const double ymax, const double ysteps,
                                  const double zmin, const double zmax,
                                  const double zsteps);

  /// This will load a previously calculated grid of time-dependant weighting
  /// potential values.
  void LoadWeightingPotentialGrid(const std::string &label) {
    for (auto &electrode : m_readout_p) {
      if (electrode.label == label) {
        if (electrode.grid.LoadWeightingField(label + "map", "xyz", true)) {
          std::cerr
              << m_className
              << "::LoadWeightingPotentialGrid: Weighting potential set for "
              << label << ".\n";
          electrode.m_usegrid = true;
          return;
        }
      }
    }
    std::cerr << m_className
              << "::LoadWeightingPotentialGrid: Could not find file for "
              << label << ".\n";
  }

  Medium *GetMedium(const double x, const double y, const double z) override;

  bool GetBoundingBox(double &xmin, double &ymin, double &zmin, double &xmax,
                      double &ymax, double &zmax) override;

  // Obtain the index and permitivity of of the layer at hight z.
  bool getLayer(const double y, int &m, double &epsM) {

    int mholer = -1;

    for (int i = 1; i < m_N; i++) {
      if (y <= m_z[i]) {
        mholer = i;
        break;
      }
    }
    if (mholer == -1)
      return false;

    m = mholer;
    epsM = m_epsHolder[m - 1];
    return true;
  }
  int NumberOfLayers() { return m_N - 1; }

private:
  static constexpr double m_precision = 1.e-30;
  static constexpr double m_Vw = 1.;
  double m_eps0 = 8.85418782e-3;
  // Voltage difference between the parallel plates.
  double m_V = 0.;

  double m_dt = 0.;

  int m_N = 0; ///< amount of layers

  double m_upperBoundIntegration = 30;

  std::vector<double> m_eps; ///< relative permittivity of each layer
  std::vector<double> m_epsHolder;
  std::vector<double> m_d; ///< thickness of each layer
  std::vector<double> m_dHolder;
  std::vector<double> m_z; ///< list of indices of conducting layers

  std::vector<int> m_sigmaIndex; ///< list of indices of conducting layers

  TF2 m_hIntegrand;

  TF1 m_wpStripIntegral; ///< Weighting potential integrant for strips
  TF2 m_wpPixelIntegral; ///< Weighting potential integrant for pixels

  std::vector<std::vector<std::vector<int>>>
      m_sigmaMatrix; // sigma_{i,j}^n, where n goes from 1 to N;
  std::vector<std::vector<std::vector<int>>>
      m_thetaMatrix; // theta_{i,j}^n, where n goes from 1 to N;

  std::vector<std::vector<double>> m_cMatrix; ///< c-matrixl.
  std::vector<std::vector<double>> m_vMatrix; ///< v-matrixl.
  std::vector<std::vector<double>> m_gMatrix; ///< g-matrixl.
  std::vector<std::vector<double>> m_wMatrix; ///< w-matrixl.

  int m_currentLayer = 0; ///< Index of the current layer.

  Medium *m_medium = nullptr;

  /// Structure that captures the information of the electrodes under study
  struct Electrode {
    std::string label;                    ///< Label.
    int ind = structureelectrode::NotSet; ///< Readout group.
    double xpos, ypos;                    ///< Coordinates in x/y.
    double lx, ly;                        ///< Dimensions in the x-y plane.
    bool formAnode = true;                ///< Dimensions in the x-y plane.

    bool m_usegrid = false; ///< Enabeling grid based calculations.
    ComponentGrid grid;     ///< grid object.
  };

  enum fieldcomponent { xcomp = 0, ycomp, zcomp };

  /// Possible readout groups
  enum structureelectrode { NotSet = -1, Plane, Strip, Pixel };

  // Vectors storing the readout electrodes.
  std::vector<std::string> m_readout;
  std::vector<Electrode> m_readout_p;

  // Functions that calculate the weighting potential

  double IntegratePromptPotential(const Electrode &el, const double x,
                                  const double y, const double z);

  void CalculateDynamicalWeightingPotential(const Electrode &el);

  double FindWeightingPotentialInGrid(Electrode &el, const double x,
                                      const double y, const double z);

  // function returning 0 if layer with specific index is conductive.
  double kroneckerDelta(const int index) {
    if (std::find(m_sigmaIndex.begin(), m_sigmaIndex.end(), index) !=
        m_sigmaIndex.end())
      return 0;
    else
      return 1;
  }

  // function construct the sigma matrix needed to calculate the w, v, c and g
  // matrices
  bool Nsigma(int N, std::vector<std::vector<int>> &sigmaMatrix);

  // function construct the theta matrix needed to calculate the w, v, c and g
  // matrices
  bool Ntheta(int N, std::vector<std::vector<int>> &thetaMatrix,
              std::vector<std::vector<int>> &sigmaMatrix);

  // function constructing the sigma an theta matrices.
  void constructGeometryMatrices(const int N);

  // function connstructing the w, v, c and g matrices needed for constructing
  // the weighting potentials equations.
  void constructGeometryFunction(const int N);

  // build function h needed for the integrant of the weighting potential of a
  // stip and pixel
  void setHIntegrand();

  // build integrant of weighting potential of a strip
  void setwpPixelIntegrand();

  // build integrant of weighting potential of a pixel
  void setwpStripIntegrand();

  // weighting field of a plane in layer with index "indexLayer"
  double constWEFieldLayer(const int indexLayer) {
    double invEz = 0;
    for (int i = 1; i <= m_N - 1; i++) {
      invEz += (m_z[i] - m_z[i - 1]) / m_epsHolder[i - 1];
    }
    return 1 / (m_epsHolder[indexLayer - 1] * invEz);
  }

  // weighting potential of a plane
  double wpPlane(const double z) {
    int im = -1;
    double epsM = -1;
    if (!getLayer(z, im, epsM))
      return 0.;
    double v = 1 - (z - m_z[im - 1]) * constWEFieldLayer(im);
    for (int i = 1; i <= im - 1; i++) {
      v -= (m_z[i] - m_z[i - 1]) * constWEFieldLayer(i);
    }

    return v;
  }

  // electric field in layer with index "indexLayer"
  double constEFieldLayer(const int indexLayer) {
    if (kroneckerDelta(indexLayer) == 0)
      return 0.;
    double invEz = 0;
    for (int i = 1; i <= m_N - 1; i++) {
      invEz += -(m_z[i] - m_z[i - 1]) * kroneckerDelta(i) / m_epsHolder[i - 1];
    }
    return m_V / (m_epsHolder[indexLayer - 1] * invEz);
  }

  // function to convert decimal to binary expressed in n digits.
  bool decToBinary(int n, std::vector<int> &binaryNum);

  // Rebuilds c, v, g and w matrix.
  void LayerUpdate(const double z, const int im, const double epsM) {

    if (im != m_currentLayer) {
      m_currentLayer = im;
      for (int i = 0; i < im - 1; i++)
        m_eps[i] = m_epsHolder[i];
      m_eps[im - 1] = epsM;
      m_eps[im] = epsM;
      for (int i = im + 1; i < m_N; i++)
        m_eps[i] = m_epsHolder[i - 1];
    }

    double diff1 = m_z[im] - z;
    double diff2 = z - m_z[im - 1];

    for (int i = 0; i < im - 1; i++)
      m_d[i] = m_dHolder[i];
    m_d[im - 1] = diff2;
    m_d[im] = diff1;
    for (int i = im + 1; i < m_N; i++)
      m_d[i] = m_dHolder[i - 1];
    // TODO::Construct c and g matrices only for im != m_currentLayer.
    constructGeometryFunction(m_N);
  };

  void UpdatePeriodicity() override;
  void Reset() override;
};
} // namespace Garfield
#endif
