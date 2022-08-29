#ifndef G_COMPONENT_GRID_H
#define G_COMPONENT_GRID_H

#include <vector>

#include "Component.hh"

namespace Garfield {

/// Component for interpolating field maps on a regular mesh.

class ComponentGrid : public Component {
 public:
  /// Constructor
  ComponentGrid();
  /// Destructor
  ~ComponentGrid() {}

  /** Define the grid.
   * \param nx,ny,nz number of nodes along \f$x, y, z\f$.
   * \param xmin,xmax range along \f$x\f$.
   * \param ymin,ymax range along \f$y\f$.
   * \param zmin,zmax range along \f$z\f$.
   */
  bool SetMesh(const unsigned int nx, const unsigned int ny,
               const unsigned int nz, const double xmin, const double xmax,
               const double ymin, const double ymax, const double zmin,
               const double zmax);
  /// Retrieve the parameters of the grid.
  bool GetMesh(unsigned int& nx, unsigned int& ny, unsigned int& nz,
               double& xmin, double& xmax, double& ymin, double& ymax,
               double& zmin, double& zmax) const;
  /// Use Cartesian coordinates (default).
  void SetCartesianCoordinates() { m_coordinates = Coordinates::Cartesian; }
  /// Use cylindrical coordinates.
  void SetCylindricalCoordinates();

  /** Import electric field and potential values from a file.
   * The file is supposed to contain one line for each grid point starting with
   *   - either two or three floating point numbers,
   *     specifying the coordinates (in cm) of the grid node or
   *   - two or three integers specifying the index of the node,
   *
   * followed by
   *   - two or three floating point numbers for the electric field (in V/cm),
   * and (depending on the value of withPotential and withFlag),
   *   - a floating point number specifying the potential (in V), and
   *   - an integer flag indicating whether the point is in an active region (1)
   *     or not (0).
   *
   * Format types are:
   *  - "xy", "xz", "xyz": nodes are specified by their coordinates
   *  - "ij", "ik", "ijk": nodes are specified by their indices
   * 
   * If cylindrical coordinates are used, the first coordinate (x)
   * corresponds to the radial distance and the second coordinate (y)
   * corresponds to the azimuth (in radian).
   */
  bool LoadElectricField(const std::string& filename, const std::string& format,
                         const bool withPotential, const bool withFlag,
                         const double scaleX = 1., const double scaleE = 1.,
                         const double scaleP = 1.);

  /// Import (prompt) weighting field from file.
  bool LoadWeightingField(const std::string& filename,
                          const std::string& format, const bool withPotential,
                          const double scaleX = 1., const double scaleE = 1.,
                          const double scaleP = 1.);
  /// Import delayed weighting field from file.
  bool LoadWeightingField(const std::string& filename,
                          const std::string& format, const double time,
                          const bool withPotential, const double scaleX = 1.,
                          const double scaleE = 1., const double scaleP = 1.);
  /// Offset coordinates in the weighting field, such that the
  /// same numerical weighting field map can be used for electrodes at
  /// different positions.
  void SetWeightingFieldOffset(const double x, const double y, const double z);

  /// Import magnetic field values from a file.
  bool LoadMagneticField(const std::string& filename, const std::string& format,
                         const double scaleX = 1., const double scaleB = 1.);

  /** Export the electric field and potential of a component to a text file.
   * \param cmp Component object for which to export the field/potential
   * \param filename name of the text file
   * \param fmt format string, see @ref LoadElectricField
   */
  bool SaveElectricField(Component* cmp, const std::string& filename,
                         const std::string& fmt);
  /** Export the weighting field and potential of a component to a text file.
   * \param cmp Component object for which to export the field/potential
   * \param id identifier of the weighting field
   * \param filename name of the text file
   * \param fmt format string, see @ref LoadElectricField
   */
  bool SaveWeightingField(Component* cmp, const std::string& id,
                          const std::string& filename,
                          const std::string& fmt);

  /// Return the field at a given node.
  bool GetElectricField(const unsigned int i, const unsigned int j,
                        const unsigned int k, double& v, double& ex, double& ey,
                        double& ez) const;

  /// Set the medium.
  void SetMedium(Medium* m);
  /// Get the medium.
  Medium* GetMedium() const { return m_medium; }

  /// Print information about the mesh and the available data.
  void Print();

  /** Import electron attachment coefficients from a file.
   * \param fname name of the text file.
   * \param fmt format string, see @ref LoadElectricField.
   * \param col column in the file which has the attachment coefficient.
   * \param scaleX scaling factor to be applied to the coordinates.
   */ 
  bool LoadElectronAttachment(const std::string& fname, 
                              const std::string& fmt,
                              const unsigned int col, 
                              const double scaleX = 1.);
  /// Import hole attachment coefficients from a file.
  bool LoadHoleAttachment(const std::string& fname, 
                          const std::string& fmt,
                          const unsigned int col,
                          const double scaleX = 1.);

  /** Import a map of electron drift velocities from a file.
   * \param fname name of the text file.
   * \param fmt format string, see @ref LoadElectricField
   * \param scaleX scaling factor to be applied to the coordinates.
   * \param scaleV scaling factor to be applied to the velocity components.
   */ 
  bool LoadElectronVelocity(const std::string& fname, 
                            const std::string& fmt,
                            const double scaleX = 1.,
                            const double scaleV = 1.e-9);
  /// Import a map of hole drift velocities from a file.
  bool LoadHoleVelocity(const std::string& fname, 
                        const std::string& fmt,
                        const double scaleX = 1.,
                        const double scaleV = 1.e-9);

  void Clear() override { Reset(); }
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;
  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;
  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;
  void DelayedWeightingField(const double x, const double y, const double z,
                             const double t, double& wx, double& wy, double& wz,
                             const std::string& label) override;
  double DelayedWeightingPotential(const double x, const double y,
                                   const double z, const double t,
                                   const std::string& label) override;
  void MagneticField(const double x, const double y, const double z, double& bx,
                     double& by, double& bz, int& status) override;

  Medium* GetMedium(const double x, const double y, const double z) override;

  bool GetVoltageRange(double& vmin, double& vmax) override;
  bool GetElectricFieldRange(double& exmin, double& exmax, double& eymin,
                             double& eymax, double& ezmin, double& ezmax);
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, 
                      double& xmax, double& ymax, double& zmax) override;
  bool GetElementaryCell(double& xmin, double& ymin, double& zmin, 
                         double& xmax, double& ymax, double& zmax) override;

  bool HasMagneticField() const override;

  bool HasAttachmentMap() const override {
    return !(m_eAttachment.empty() && m_hAttachment.empty());
  } 
  bool ElectronAttachment(const double x, const double y, const double z,
                          double& att) override;
  bool HoleAttachment(const double x, const double y, const double z,
                      double& att) override;

  bool HasVelocityMap() const override {
    return !(m_eVelocity.empty() && m_hVelocity.empty());
  } 
  bool ElectronVelocity(const double x, const double y, const double z,
                        double& vx, double& vy, double& vz) override;
  bool HoleVelocity(const double x, const double y, const double z,
                    double& vx, double& vy, double& vz) override;
 private:
  enum class Format {
    Unknown,
    XY,
    XZ,
    XYZ,
    IJ,
    IK,
    IJK,
    YXZ
  };
  enum class Coordinates {
    Cartesian,
    Cylindrical
  };
  Coordinates m_coordinates = Coordinates::Cartesian;
 
  Medium* m_medium = nullptr;
  struct Node {
    double fx, fy, fz;  ///< Field
    double v;           ///< Potential
  };

  /// Electric field values and potentials.
  std::vector<std::vector<std::vector<Node> > > m_efields;
  /// Magnetic field values.
  std::vector<std::vector<std::vector<Node> > > m_bfields;
  /// Prompt weighting field values and potentials.
  std::vector<std::vector<std::vector<Node> > > m_wfields;
  /// Delayed weighting field values and potentials.
  std::vector<std::vector<std::vector<std::vector<Node> > > > m_wdfields;
  std::vector<double> m_wdtimes;
  /// Attachment maps for electrons and holes.
  std::vector<std::vector<std::vector<double> > > m_eAttachment;
  std::vector<std::vector<std::vector<double> > > m_hAttachment;
  /// Velocity maps for electrons and holes.
  std::vector<std::vector<std::vector<Node> > > m_eVelocity;
  std::vector<std::vector<std::vector<Node> > > m_hVelocity;
  /// Active medium flag.
  std::vector<std::vector<std::vector<bool> > > m_active;

  // Dimensions of the mesh
  std::array<unsigned int, 3> m_nX = {{1, 1, 1}};
  std::array<double, 3> m_xMin = {{0., 0., 0.}};
  std::array<double, 3> m_xMax = {{0., 0., 0.}};
  std::array<double, 3> m_sX = {{0., 0., 0.}};

  bool m_hasMesh = false;
  bool m_hasPotential = false;

  // Offset for weighting field
  std::array<double, 3> m_wFieldOffset = {{0., 0., 0.}};

  // Voltage range
  double m_pMin = 0., m_pMax = 0.;

  /// Read/determine mesh parameters from file.
  bool LoadMesh(const std::string& filename, std::string format,
                const double scaleX);

  /// Read electric field and potential from file.
  bool LoadData(const std::string& filename, std::string format,
                const bool withPotential, const bool withFlag,
                const double scaleX, const double scaleF, const double scaleP,
                std::vector<std::vector<std::vector<Node> > >& field);
  /// Load scalar data (e. g. attachment coefficients) from file. 
  bool LoadData(const std::string& filename, std::string format,
                const double scaleX,
                std::vector<std::vector<std::vector<double> > >& tab,
                const unsigned int col);

  void Reset() override;
  void UpdatePeriodicity() override;

  /// Interpolation of the field and potential at a given point.
  bool GetField(const double x, const double y, const double z,
                const std::vector<std::vector<std::vector<Node> > >& field,
                double& fx, double& fy, double& fz, double& p, bool& active);
  /// Interpolation in a table of scalars.
  bool GetData(const double x, const double y, const double z,
      const std::vector<std::vector<std::vector<double> > >& table,
      double& value);

  /// Reduce a coordinate to the basic cell (in case of periodicity).
  double Reduce(const double xin, const double xmin, const double xmax,
                const bool simplePeriodic, const bool mirrorPeriodic,
                bool& isMirrored) const;
  /// Set the dimensions of a table according to the mesh.
  void Initialise(std::vector<std::vector<std::vector<Node> > >& fields);
  /// Decode a format string.
  Format GetFormat(std::string fmt);
};
}  // namespace Garfield
#endif
